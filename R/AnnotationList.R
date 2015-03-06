setOldClass("table")
#' The annotation data.table should have columns "accession" and "category"
AnnotationList<-setClass("AnnotationList",
                         slots=c(
                           proteins="character",
                           categories="factor",
                           annotationWeight="data.table",
                           annotation="table"
                           ),
                         prototype=prototype(
                           annotationWeight=structure(list(), class="data.table"),
                           annotation=structure(list(), class="table")
                         )
                         )

#' Retrieve Annotation
#' 
#' \code{retrieveAnnotation} Uses a data source (depending on 
#' implementation) to retrieve annotations for the
#' defined set of proteins. The result is stored internally.
#' 
#' For the case of deep DAGs like the gene ontology, this method
#' is expected to provide the entire list of GOs from root to tip
#' that can be assigned to the proteins.
#' 
#' @param object An object which inherits from AnnotationList
#' 
#' @return A new object of the same class with the annotation filled.
#' @export
setGeneric("retrieveAnnotation",function(object) standardGeneric("retrieveAnnotation"))

#' Annotation weights
#' 
#' \code{annotationWeights} provides the weights for each annotation
#' category. In the case of GO, this could be the distance from
#' the root object. Given this, all weights should be integers > 0. 
#' In terms of the LP problem, the returned vector is equivalent to Dj.
#' The implementation needs to make sure that the order returned is consistent
#' with other outputs.
#' 
#' @param object An object which inherits from AnnotationList
#' 
#' @return The data.table Dj holding the annotation weights.
#' @export
setGeneric("annotationWeights", function(object) standardGeneric("annotationWeights"))

#' Annotation Matrix
#' 
#' \code{annotationMatrix} produces the Cij matrix, where position i,j is 1 if gene i
#' can be annotated with category j.
#' 
#' @param object An object which inherits from AnnotationList
#' 
#' @return matrix A matrix of n=number of proteins rows and m=number of categories 
#' columns
setGeneric("annotationMatrix", function(object) standardGeneric("annotationMatrix"))