# HCCS
Each sample would be partitioned into the corresponding subtype based on the nearest centroid method and Pearson correlation.


devtools::install_github('Zaoqu-Liu/HCCS')

library(HCCS)

library(GSVA)

time_phenotype <- getTIME(data)

ferroptosis_phenotype <- getFerroptosis(data)

hypoxia_phenotype <- getHypoxia(data)
