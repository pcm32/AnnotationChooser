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
  
  A1<-matrix(0,rowsA1,colsA1)
  
  i<-0
  j<-1
  tj<-0
  for(col in 1:colsA1) {
    if(i==geneAmount) {
      i<-1
      j<-j+1
    } else {
      i<-i+1
    }
    
    for(row in 1:rowsA1) {
      if( col <= rowsA1 ) {
        if( row == col ) {
          A1[row, col]<-1
        }
      }
    }
  }
  
  ## Senses1
  
  senses1<-c(rep('<=',times = rowsA1))
  
  ## b1
  
  b1<-c(c_ij)
  
  ## A2
  
  rowsA2<-geneAmount
  colsA2<-colsA1
  A2<-matrix(0,rowsA2,colsA2)
  
  i<-0
  j<-1
  tj<-0
  
  for(col in 1:colsA2) {
    if(i==geneAmount) {
      i<-1
      j<-j+1
    } else {
      i<-i+1
    }
    
    for(row in 1:rowsA2) {
      if( col <= geneAmount*ontSize ) {
        if( row == i ) {
          A2[row, col]<-1
        }
      }
    }
  }
  
  ## senses2
  
  senses2<-c(rep('<=',times = rowsA2))
  
  ## b2
  
  b2<-c(rep(1,times = geneAmount))
  
  ## A3
  
  rowsA3<-ontSize
  colsA3<-colsA1
  A3<-matrix(0,rowsA3,colsA3)
  
  i<-0
  j<-1
  tj<-0
  
  for(col in 1:colsA3) {
    if(col <= geneAmount*ontSize) {
      if(i==geneAmount) {
        i<-1
        j<-j+1
      } else {
        i<-i+1
      }
    } else {
      tj<-col-(geneAmount*ontSize)
    }
    
    for(row in 1:rowsA3) {
      if( col <= geneAmount*ontSize ) {
        if( row == j ) {
          A3[row, col]<-1
        }
      } else {
        if( row == tj) {
          A3[row, col]<-geneAmount*-1
        }
      }
    }
  }
  
  ## senses3
  
  senses3<-c(rep('<=',times = rowsA3))
  
  ## b3
  
  b3<-c(rep(0,times = ontSize))
  
  ## A5
  
  rowsA5<-ontSize
  colsA5<-colsA1
  A5<-matrix(0,rowsA5,colsA5)
  
  i<-0
  tj<-0
  
  for(col in (geneAmount*ontSize+1):colsA5) {
    tj<-col-(geneAmount*ontSize)
    
    for(row in 1:rowsA5) {
      if( row == tj) {
        A5[row, col]<-(minDistance - d_j[tj])
      }
    }
  }
  
  ## senses5
  
  senses5<-c(rep('<=',times = rowsA5))
  
  ## b5
  
  b5<-c(rep(0,times = ontSize))
  
  ## A4
  
  rowsA4<-1
  colsA4<-colsA1
  A4<-matrix(0,rowsA4,colsA4)
  
  for(col in ((geneAmount*ontSize)+1):colsA4) {
    A4[1,col]<-1
  }
  
  ## senses4
  
  senses4<-c(rep('<=',times = rowsA4))
  
  ## b4
  
  b4<-c(maxTerms)
  
  ## Put everything together
  
  objFun<-c(rep(d_j,times = geneAmount),rep(0,times = ontSize))
  rbind(A1, A2, A3, A4, A5)->Amatrix
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