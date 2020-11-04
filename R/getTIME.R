#'@title  Get the TIME phenotype and TIMEscore
#'
#'@description  Use gene expression data to get the corresponding TIME phenotype and TIMEscore of each patient
#'
#'@details  Make sure that the data you enter is a matrix, and rows are samples and columns are genes.
#'
#'@param data The matrix of gene expression data, rows are samples and columns are genes.
#'
#'@return return a dataframe, the first column is the ID of each sample, the second column is the
#'TIME phenotype of each sample, and the third column is the TIME index of each sample.
#'
#'@export
getTIME <- function(data=data){
  ssgsea <- t(scale(t((gsva(as.matrix(data),immune,method = 'ssgsea')))))
  x <- intersect(rownames(ssgsea),rownames(TIME))
  ssgsea <- ssgsea[x,]
  Centroids <- TIME[x,]
  tt <- data.frame()
  for (i in 1:ncol(ssgsea)) {
    c1 <- cor(ssgsea[,i],Centroids[,1])
    c2 <- cor(ssgsea[,i],Centroids[,2])
    c3 <- cor(ssgsea[,i],Centroids[,3])
    if(c1>c2&c1>c3)
      tt <- rbind(tt,data.frame(ID=colnames(ssgsea)[i],Cluster='TIME-1'))
    if(c2>c1&c2>c3)
      tt <- rbind(tt,data.frame(ID=colnames(ssgsea)[i],Cluster='TIME-2'))
    if(c3>c1&c3>c2)
      tt <- rbind(tt,data.frame(ID=colnames(ssgsea)[i],Cluster='TIME-3'))
  }
  TI <- as.data.frame(t(gsva(as.matrix(data),TI,method='ssgsea')))
  tt <- merge(tt,TI,by.x = 1,by.y = 0)
  colnames(tt)[3] <- 'TI'
  return(tt)
}
