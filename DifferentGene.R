library('gplots')
library('limma')
setwd("F:\\8.AgeModel\\1.data\\Exp\\fpkm")

mydata <- read.table("BRCAUniGeneSample.txt",header=T,sep='\t',row.names=1)
mygroup <- read.table("F:\\8.AgeModel\\GSE108213Test\\1.data\\CancerSampleGroup.txt",header=T,sep='\t',row.names=1)
mygroup <- mygroup[colnames(mydata),]
group <- matrix(,dim(mygroup)[1],1)
group[which(mygroup[,4]=="Low"),1]="Low"
group[which(mygroup[,4]!="Low"),1]="Control"
rownames(group) <- rownames(mygroup)
design <- model.matrix(~0+group)
colnames(design)=levels(factor(group))
rownames(design)=rownames(group)
contrast.matrix<-makeContrasts(paste(c("Low","Control"),collapse = "-"),levels = design)

##step1
fit <- lmFit(mydata,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
length(which(2^nrDEG[,1]>1.5 & nrDEG[,5]<0.05))

setwd("F:\\8.AgeModel\\GSE108213Test\\2.DiffExp")
write.table(nrDEG,"LowFPKMlimmaResult.txt",sep='\t',quote=F)