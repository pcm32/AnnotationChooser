#' Produce Gurobi Model
#' 
#' Produces a formatted gurobi model based on the annotation list object given and 
#' the additional parameters that control the problem definition.
#' 
#' @param annotationList An AnnotationList object, from where the annotation matrix and
#' dimensions are obtained.
#' @param maxAnnotationTerms An integer with the maximum numbers of terms that should be used
#' for the overall annotation. Defaults to 20.
#' @param minDistanceFromRoot An integer with the minimum distance that terms need to have from
#' the root of the ontology to be eligible. Defaults to 3.
#' 
#' @return the gurobi model, ready to be obtimized using gurobi(). It is a list containing
#' the A matrix, the objective function, the b vector, the senses vector, the right hand side,
#' the nature of the variables and the type (maximization.)
produceGurobiModel<-function(annotationList,maxAnnotationTerms=20,minDistanceFromRoot=3) {
  annotMatrix<-annotationMatrix(annotationList)
  geneAmount<-nrow(annotMatrix)
  ontSize<-ncol(annotMatrix)
  c_ij<-as.vector(annotMatrix) #sample(0:1,size = geneAmount*ontSize, replace = TRUE)
  
  d_j<-annotationWeights(annotationList)$depth #sample(3:8,size = ontSize, replace=TRUE)
  
  maxTerms<-maxAnnotationTerms
  minDistance<-minDistanceFromRoot
  
  # All the following parts should be re-written to use the sparse matrix packages.
  # package slam 
  
  ## A1
  
  rowsA1<-geneAmount*ontSize
  colsA1<-rowsA1+ontSize
  
  #A1<-matrix(0,rowsA1,colsA1)
  #A1<-simple_triplet_zero_matrix(nrow = rowsA1, ncol = colsA1, mode = "integer")
  A1<-Matrix(0,rowsA1,colsA1,sparse = TRUE)
  
  i<-0
  j<-1
  tj<-0
  # we don't need to look at all columns, only those below rowsA1
  for(col in 1:min(colsA1,rowsA1))  {
    # For this part is not being used.
    #if(i==geneAmount) {
    #  i<-1
    #  j<-j+1
    #} else {
    #  i<-i+1
    #}
    
    # Since we know which row we want, we don't need to go through all of them!!!
    #for(row in 1:rowsA1) {
    #  if( col <= rowsA1 ) {
    #    if( row == col ) {
    #      A1[row, col]<-1
    #    }
    #  }
    #}
    # Instead we can do:
    # And we can avoid the if if we only visit the required columns
    #if(col <= rowsA1) {
      A1[col, col]<-1
    #}
  }
  
  ## Senses1
  
  senses1<-c(rep('<=',times = rowsA1))
  
  ## b1
  
  b1<-c(c_ij)
  
  ## A2
  
  rowsA2<-geneAmount
  colsA2<-colsA1
  #A2<-matrix(0,rowsA2,colsA2)
  #A2<-simple_triplet_zero_matrix(nrow = rowsA2, ncol = colsA2)
  A2<-Matrix(0,rowsA2,colsA2,sparse = TRUE)
  
  i<-0
  j<-1
  tj<-0
  
  # we only need to visit cols that are smaller than geneAmount*ontSize
  for(col in 1:min(colsA2,geneAmount*ontSize) ) {
    if(i==geneAmount) {
      i<-1
      #j<-j+1 # not used here
    } else {
      i<-i+1
    }
    
    # Since we know which row we want, we don't need to go through all of them!!!
    #for(row in 1:rowsA2) {
    #  if( col <= geneAmount*ontSize ) {
    #    if( row == i ) {
    #      A2[row, col]<-1
    #    }
    #  }
    #}
    # Instead we can do:
    # if we visit the correct columns, we don't need the if
    #if(col <= geneAmount*ontSize) {
      A2[i,col]<-1
    #}
  }
  
  ## senses2
  
  senses2<-c(rep('<=',times = rowsA2))
  
  ## b2
  
  b2<-c(rep(1,times = geneAmount))
  
  ## A3
  
  rowsA3<-ontSize
  colsA3<-colsA1
  #A3<-matrix(0,rowsA3,colsA3)
  #A3<-simple_triplet_zero_matrix(nrow = rowsA3, ncol = colsA3)
  A3<-Matrix(0,rowsA3,colsA3,sparse = TRUE)
  
  i<-0
  j<-1
  tj<-0
  
  # we should split this for into two parts, 1 to geneAmount*ontSize and geneAmount*ontSize+1 to colsA3
#   for(col in 1:colsA3) {
#     if(col <= geneAmount*ontSize) {
#       if(i==geneAmount) {
#         i<-1
#         j<-j+1
#       } else {
#         i<-i+1
#       }
#     } else {
#       tj<-col-(geneAmount*ontSize)
#     }
#     
#     # Since we know which row we want, we don't need to go through all of them!!!
#     #for(row in 1:rowsA3) {
#     #  if( col <= geneAmount*ontSize ) {
#     #    if( row == j ) {
#     #      A3[row, col]<-1
#     #    }
#     #  } else {
#     #    if( row == tj) {
#     #      A3[row, col]<-geneAmount*-1
#     #    }
#     #  }
#     #}
#     # Instead we can do
#     if( col <= geneAmount*ontSize ) {
#         A3[j, col]<-1
#     } else {
#         A3[tj, col]<-geneAmount*-1
#     }
#   }
  
  # so we can re-write like this:
  for(col in 1:(geneAmount*ontSize) ) {
    if(i==geneAmount) {
      i<-1
      j<-j+1
    } else {
      i<-i+1
    }
    
    A3[j, col]<-1
  } 
  for(col in (geneAmount*ontSize+1):colsA3  ) {
    tj<-col-(geneAmount*ontSize)
    A3[tj, col]<-geneAmount*-1
  }
  
  ## senses3
  
  senses3<-c(rep('<=',times = rowsA3))
  
  ## b3
  
  b3<-c(rep(0,times = ontSize))
  
  ## A5
  
  rowsA5<-ontSize
  colsA5<-colsA1
  #A5<-matrix(0,rowsA5,colsA5)
  #A5<-simple_triplet_zero_matrix(nrow = rowsA5, ncol = colsA5)
  A5<-Matrix(0,rowsA5,colsA5,sparse = TRUE)
  
  i<-0
  tj<-0
  
  for(col in (geneAmount*ontSize+1):colsA5) {
    tj<-col-(geneAmount*ontSize)
    
    # Since we know which row we want, we don't need to go through all of them!!!
    #for(row in 1:rowsA5) {
    #  if( row == tj) {
    #    A5[row, col]<-(minDistance - d_j[tj])
    #  }
    #}
    # Instead we can do:
    A5[tj,col]<-(minDistance - d_j[tj])
  }
  
  ## senses5
  
  senses5<-c(rep('<=',times = rowsA5))
  
  ## b5
  
  b5<-c(rep(0,times = ontSize))
  
  ## A4
  
  rowsA4<-1
  colsA4<-colsA1
  #A4<-matrix(0,rowsA4,colsA4)
  #A4<-simple_triplet_zero_matrix(nrow = rowsA4, ncol = colsA4)
  A4<-Matrix(0,rowsA4,colsA4,sparse=TRUE)
  
  for(col in ((geneAmount*ontSize)+1):colsA4) {
    A4[1,col]<-1
  }
  
  ## senses4
  
  senses4<-c(rep('<=',times = rowsA4))
  
  ## b4
  
  b4<-c(maxTerms)
  
  ## Put everything together
  
  objFun<-c(rep(d_j,times = geneAmount),rep(0,times = ontSize))
  #rbind(A1, A2, A3, A4, A5)->Amatrix
  rBind(A1, A2, A3, A4, A5)->Amatrix
  c(senses1,senses2,senses3,senses4,senses5)->senses
  c(b1,b2,b3,b4,b5)->bVector
  
  model<-list()
  model$A<-Amatrix
  model$obj<-objFun
  model$sense<-senses
  model$rhs<-bVector
  model$vtypes<-c(rep('B',times = ncol(Amatrix)))
  model$modelsense<-'max'
  
  return(model)
}

#' Obtain annotated solution
#' 
#' Produces an annotated version of the optimal within the gurobi result given (the result of calling gurobi()).
#' The annotation is extracted from the annotation list object.
#' 
#' @param gurobiRes The result of invoking gurobi() with the model produced by \code{produceGurobiModel}.
#' @param annotationList The same AnnotationList object provided to \code{produceGurobiModel}.
#' 
#' @return a list containing \code{resMat}, which contains the result matrix (Terms as rows, Proteins as Columns, 
#' a 1 in i,j means Term i is annotation Protein j), and \code{coverage}, which contains a summary per Term of 
#' the maximum depth of the term and the number of proteins it is annotating (the coverage).
#' 
obtainAnnotatedSolution<-function(gurobiRes,annotationList) {
  annotMatrix<-annotationMatrix(annotationList)
  geneAmount<-nrow(annotMatrix)
  ontSize<-ncol(annotMatrix)
  resMat<-matrix(gurobiRes$x[1:(geneAmount*ontSize)],nrow=ontSize, byrow=TRUE)
  rownames(resMat)<-colnames(annotMatrix)
  colnames(resMat)<-rownames(annotMatrix)
  
  annotationWeights(annotationList)->annotWeights
  
  solution<-list()
  solution$resMat<-resMat[rowSums(resMat)>0,]
  solution$coverage<-annotWeights[term %in% rownames(solution$resMat),
                                  list(term,depth,coverage=rowSums(solution$resMat[rownames(solution$resMat)==term,])),
                                  ]
  
  return(solution)
}