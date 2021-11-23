setwd("D:\\1.BMAP\\1.data")
mydata <- read.table("BRCA_gene_tpm.txt",header=T,sep='\t',row.names=1,check.names=F)
library(pRRophetic)
library(ggplot2)

setwd("D:\\1.BMAP\\35.preDrug")
predictedPtype_ccle <- pRRopheticPredict(as.matrix(mydata), "Tamoxifen",tissueType="breast",selection=1)
write.table(predictedPtype_ccle,"predictedPtype_Methotrexate.txt",sep='\t',quote=F,col.names=F)