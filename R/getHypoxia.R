#'@title  Get the Hypoxia phenotype and HARS
#'
#'@description  Use gene expression data to get the corresponding Hypoxia phenotype and Hypoxia associated risk score of each patient
#'
#'@details  Make sure that the data you enter is a matrix, and rows are samples and columns are genes.
#'
#'@param data The matrix of gene expression data, rows are samples and columns are genes.
#'
#'@return return a dataframe, the first column is the ID of each sample, the second column is the
#'Hypoxia phenotype of each sample, and the third column is the HARS of each sample.
#'
#'@export
getHypoxia <- function(data=data){
  centroid <- Hypoxia
  x <- intersect(rownames(data),rownames(centroid))
  data2 <- data[x,]
  centroid <- centroid[x,]
  tt <- data.frame()
  for (i in 1:ncol(data2)) {
    c1 <- cor(data2[,i],centroid[,1])
    c2 <- cor(data2[,i],centroid[,2])
    if(c1>c2)
      tt <- rbind(tt,data.frame(ID=colnames(data2)[i],Cluster='C1'))
    if(c2>c1)
      tt <- rbind(tt,data.frame(ID=colnames(data2)[i],Cluster='C2'))
  }
  gg <- c("KIF15","KIF23","MCM10","CEP55","CENPA","NEIL3","GTSE1","KIF18B","RAD54L","KIF18A")
  if(sum(rownames(data)%in%gg)<10){
    stop(paste0("Absent ",setdiff(gg,rownames(data))))
  }else{
    data.pca <- prcomp(t(data[gg,]), scale. = TRUE)
    pcaPredict=as.data.frame(predict(data.pca))
    tt$HARS <- pcaPredict$PC1
  }
  return(tt)
}
