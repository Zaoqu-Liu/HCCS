#'@title  Get the Ferroptosis phenotype and FRRS
#'
#'@description  Use gene expression data to get the corresponding Ferroptosis phenotype and Ferroptosis related risk score of each patient
#'
#'@details  Make sure that the data you enter is a matrix, and rows are samples and columns are genes.
#'
#'@param data The matrix of gene expression data, rows are samples and columns are genes.
#'
#'@return return a dataframe, the first column is the ID of each sample, the second column is the
#'Ferroptosis phenotype of each sample, and the third column is the FRRS of each sample.
#'
#'@export
getFerroptosis <- function(data=data){
  centroid <- Ferroptosis
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
  if(sum(rownames(data)%in%c('SLC16A3','CPS1'))<2){
    stop(paste0("Absent ",setdiff(c('SLC16A3','CPS1'),rownames(data))))
  }else{
    tt$FRRS <- rowSums(t(data[c('SLC16A3','CPS1'),])*c(0.348,-0.151))
  }
  return(tt)
}



