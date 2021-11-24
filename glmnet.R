setwd("F:\\8.AgeModel\\GSE108213Test\\1.data")
library(glmnet)
mydata <- read.table("TCGABRCAinterCG.txt",sep='\t',header=T,row.names=1)

set.seed(1234)
data <- t(mydata[,c(-45,-60,-80,-96)]) #remove samples
RandomSam<-sample(rownames(data),46)
TestSam<- setdiff(rownames(data), RandomSam)
Train_x <- as.matrix(data[RandomSam,2:length(data[1,])])
Train_y <- as.numeric(data[RandomSam,1])
Test_x <- as.matrix(data[TestSam,2:length(data[1,])])
Test_y <- as.numeric(data[TestSam,1])

cv.fit = cv.glmnet(Train_x,Train_y, family = "gaussian",nfolds=10,type.measure="mae",alpha=1)
fit <- glmnet(Train_x,Train_y, family = "gaussian",alpha=1)
coefficients<-coef(fit,s=cv.fit$lambda.min)
Active.Index<-which(as.numeric(coefficients)!=0)
Active.coefficients<-coefficients[Active.Index]
length(Active.Index)

predict_Train <- as.numeric(predict(fit,newx=Train_x,s=cv.fit$lambda.min))
TrainMAE=mean(abs(as.numeric(predict_Train)-Train_y))
TrainMAE
predict_Test <- as.numeric(predict(fit,newx=Test_x,s=cv.fit$lambda.min))
TestMAE=mean(abs(as.numeric(predict_Test)-Test_y))
TestMAE

CGmatrix <- as.matrix(coef(fit,s=cv.fit$lambda.min))
CGindex <- CGmatrix[which(CGmatrix[,1]!=0),1]
write.table(CGindex,"CGindex.txt",sep='\t',quote=F,col.names=F)


data1 <- Train_x[,names(CGindex)[2:length(CGindex)]]
data2 <- matrix(,dim(data1)[1],dim(data1)[2])
for (i in 1:length(CGindex)-1){
	data2[,i] <- data1[,i]*as.numeric(CGindex[i+1])
}
DNAm <- rowSums(data2)+as.numeric(CGindex[1])
mean(abs(as.numeric(Train_y)-DNAm))
AgeData <- cbind(Train_y,DNAm)
rownames(AgeData) <- rownames(Train_x)
write.table(AgeData,"TrainAge.txt",quote=F,sep='\t')

data3 <- Test_x[,names(CGindex)[2:length(CGindex)]]
data4 <- matrix(,dim(data3)[1],dim(data3)[2])
for (i in 1:length(CGindex)-1){
	data4[,i] <- data3[,i]*as.numeric(CGindex[i+1])
}
TestDNAm <- rowSums(data4)+as.numeric(CGindex[1])
mean(abs(as.numeric(Test_y)-TestDNAm))
AgeDataTest <- cbind(Test_y,TestDNAm)
rownames(AgeDataTest) <- rownames(Test_x)
write.table(AgeDataTest,"TestAge.txt",quote=F,sep='\t')

#GEO data
GEOdata <- read.table("GSE108213Data.txt",header=T,sep='\t',row.names=1)

data5 <- t(GEOdata)[,names(CGindex)[2:length(CGindex)]]
data6 <- matrix(,dim(data5)[1],dim(data5)[2])
for (i in 1:length(CGindex)-1){
	data6[,i] <- data5[,i]*as.numeric(CGindex[i+1])
}
GEODNAm <- rowSums(data6)+as.numeric(CGindex[1])
AgeDataGEO <- cbind(as.numeric(GEOdata[1,]),GEODNAm)
mean(abs(as.numeric(GEOdata[1,])-GEODNAm))
rownames(AgeDataGEO) <- colnames(GEOdata)
write.table(AgeDataGEO,"GEOAge.txt",quote=F,sep='\t')
