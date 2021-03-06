% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getHypoxia.R
\name{getHypoxia}
\alias{getHypoxia}
\title{Get the Hypoxia phenotype and HARS}
\usage{
getHypoxia(data = data)
}
\arguments{
\item{data}{The matrix of gene expression data, rows are samples and columns are genes.}
}
\value{
return a dataframe, the first column is the ID of each sample, the second column is the
Hypoxia phenotype of each sample, and the third column is the HARS of each sample.
}
\description{
Use gene expression data to get the corresponding Hypoxia phenotype and Hypoxia associated risk score of each patient
}
\details{
Make sure that the data you enter is a matrix, and rows are samples and columns are genes.
}
