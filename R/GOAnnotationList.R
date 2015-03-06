#' @include AnnotationList.R
NULL

setClass("GOAnnotationList", 
         slots=c(type="character"),
         contains="AnnotationList"  
         )

#' GOAnnotationList
#' 
#' \code{GOAnnotationList} retrieves and manages dense GO annotations
#' for a list of proteins. By dense it means that it doesn't only shows the tip/leaf
#' of the gene ontology to which a protein is associated, but includes all the is_a
#' and part_of ancestors for that term being assigned to the protein.
#' 
#' The object is a sub-class of AnnotationList.
#' 
#' @param proteins The list of proteins for which the annotation should be produced.
#' @param type The branch of the Gene Ontology to use, either "MF" for molecular function,
#' "BP" for biological process, or "CC" for cellular container
#'  
#' @export GOAnnotationList   
GOAnnotationList<-function(proteins,type, ...) {
  new("GOAnnotationList",type=type,proteins=proteins)
}

#' @export
setMethod("retrieveAnnotation",signature(object="GOAnnotationList"), function(object) {
  # first, retrieve all direct annotations from either biomart or a GO annotation dbi
  branch<-NA_character_
  root<-NA_character_
  allAncestorsDB<-list()
  directChildrenDB<-list()
  if(object@type == "MF") {
    branch<-"db2go_f__dm_primary_id"
    root<-"GO:0003674"    
    allAncestorsDB<-AnnotationDbi::as.list(GOMFANCESTOR)
    directChildrenDB<-AnnotationDbi::as.list(GOMFCHILDREN)
  } else if(object@type == "BP") {
    branch<-"db2go_p__dm_primary_id"
    root<-"GO:0008150"
    allAncestorsDB<-AnnotationDbi::as.list(GOBPANCESTOR)
    directChildrenDB<-AnnotationDbi::as.list(GOBPCHILDREN)
  } else if(object@type == "CC") {
    branch<-"db2go_c__dm_primary_id"
    root<-"GO:0005575"
    allAncestorsDB<-AnnotationDbi::as.list(GOCCANCESTOR)
    directChildrenDB<-AnnotationDbi::as.list(GOCCCHILDREN)
  }
  
  mart<-useMart(biomart = "unimart",dataset = "uniprot")
  annotRes<-data.table(getBM(attributes=c("accession",branch),
                   filters=c("accession"),
                   values=object@proteins,mart=mart,uniqueRows=T))#,key=key)
  setnames(annotRes,old = c(branch), new= c("category"))
  unique(annotRes$category)->uniqueCategory
  #unique(c(unlist(allAncestorsDB[names(allAncestorsDB) %in% uniqueCategory]),uniqueCategory))->allAncestors
  
  # fill Cij matrix: can a protein be annotated by a GO. Proteins are rows, GOs are columns
  setkey(annotRes,category)
  annotRes[,
      ancestors:=paste(category,paste(unlist(unique(allAncestorsDB[names(allAncestorsDB) %in% category])),collapse=";"),sep=";"),
      by=category]
  setkey(annotRes,accession)
  annotRes[,list(ancestor=unique(unlist(strsplit(ancestors,split = ";")))),by=accession]->cij.dt
  setkey(cij.dt)
  cij.dt[ancestor!="all" & ancestor!=root,][,,by=c("accession","ancestor")]->cij.dt
  
  table(cij.dt)->object@annotation # Cij
  
  # fill in Dj, the vector of depths for each term. Make a breadth first recursion from the root using GO.db,
  # only following those elements that are present in object@allAncestors, recording their minimal distance from the root 
  
  data.table(term=colnames(object@annotation),depth=rep(0,times=ncol(object@annotation)),jindex=1:ncol(object@annotation))->Dj
  
  fillDj<-function(parentTerm,iterationDepth=1) {
    directChildrenDB[names(directChildrenDB) %in% parentTerm]->children
    unlist(children[[1]])->children
    children[names(children)=="is_a"]->children
    # we select only those of interest
    children[children %in% Dj$term]->children
    # and mark their depth, using as depth the maximum depth found for a node.
    # This condition implies that it doesn't matter if we do this depth first or breadth first.
    Dj[term %in% children & depth<iterationDepth,depth:=iterationDepth,by=term]
    for(childAsParent in children) {
      fillDj(parentTerm=childAsParent,iterationDepth=iterationDepth+1)
    }
  }
  
  fillDj(root)
  
  object@annotationWeight<-Dj
  object@categories<-as.factor(Dj$term)
  
  return(object)
})

#'@export
setMethod("annotationWeights",signature(object="GOAnnotationList"), function(object) {
  return(object@annotationWeight)
})

#'@export
setMethod("annotationMatrix",signature(object="GOAnnotationList"), function(object) {
  return(object@annotation)
})

