###enrichment
library(GOplot)
library(stringr)
library(clusterProfiler)
setwd("D:\\11.BMAP\\2.DiffExp")
enrich=read.table("medianDiffGene.txt",header=T,sep="\t",row.names=1)
gene<-rownames(enrich)

eg <- bitr(gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
gene<- eg$ENTREZID
ego<-enrichGO(gene=gene,OrgDb=org.Hs.eg.db,ont="ALL",
pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.05,readable=TRUE) #MF,CC,BP
ego2=data.frame(ego)
library(stringr)
GO=ego2[1:5,c(2,3,9,7)]
GO$geneID=str_replace_all(GO$geneID,'/',',')
names(GO)=c("ID","Term","Genes","adj_pval")
GO$Category="BP"

#read logFC
genelogFC=read.table("LowFPKMlimmaResult.txt",header=T,sep="\t",check.names=F)
names(genelogFC)[1]="ID"
genelist<-as.data.frame(genelogFC)
genelist <- genelist[,c("ID","logFC")]

plot_data=list(DA=GO,GE=genelist)
circ2=data.frame()
circ2=circle_dat(plot_data$DA,plot_data$GE)
chord<-chord_dat(data=circ2,genes=genelist,process=unique(circ2$term))
pdf("low_GO1.pdf",width=18,height=7)
pdf("low_GO.pdf",width=6,height=7)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 2,lfc.min=0.6,lfc.max=1.3,
ribbon.col=c("#ED5565","#FC6E51","#FFCE54","#48CFAD","#4FC1E9"), #GO term颜色
lfc.col=c("firebrick3", "white", "royalblue3"))
dev.off()
######KEGG
kegg <- enrichKEGG(gene, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 2,maxGSSize = 500,qvalueCutoff = 0.05,use_internal_data = FALSE)
kegg1 <- setReadable(kegg,OrgDb=org.Hs.eg.db,keyType="ENTREZID")

setwd("D:\\11.BMAP\\7.cnv")
mydata <- read.table("cnvClassCount.txt",header=T,sep='\t',check.names=F,row.names=1)

###CNV
library(ggpubr)
pdf("del.pdf",width=7,height=7)
mydata$group <- factor(mydata$group,level=c("High","Median","Low"))
compar<-list(c("High","Median"),c("High","Low"),c("Median","Low"))
ggviolin(mydata, "group", "del",
   fill = "group", palette =c("#E13220", "#DDAC2A", "#519DCC"),add = c("jitter"))+
stat_compare_means(comparisons = compar,method="wilcox.test")+stat_compare_means(label.y = 15)+
ylab("log2(deletion)")+xlab("")+theme(legend.position = "none")+
scale_y_continuous(limits = c(4,15))

dev.off()

###Stage VS TMB plot
setwd("D:\\1.BMAP\\29.PAM50")
mydata <- read.table("pam50_group.txt",header=T,sep='\t',row.names=1,check.names=F)

setwd("D:\\1.BMAP\\8.stage")
data1 <- read.table("BRCAstageClass.txt",header=F,sep='\t',row.names=1,check.names=F)

setwd("D:\\1.BMAP\\mutation")
mygroup <- read.table("BRCA_TMB.txt",header=T,sep='\t',row.names=1,check.names=F)

name1 <- intersect(rownames(data1),rownames(mygroup))
data2 <- data.frame(data1[name1,],tmb=mygroup[name1,])

data3$logTMB <- log2(data3$tmb+1)
data3 <- data2[which(data2[,4]!= "not reported" & data2[,4]!= "stage x"),]
data3 <- data2[which(data2[,2]!= "NX"),]
data3$logTMB <- log2(data3$tmb+1)
data3$V2 <- factor(data3$V2)
pdf("N_tmb.pdf",height=5,width=5)
p <- ggviolin(data3, "V3", "logTMB", 
fill = "V3", palette =c("#ED5565", "#FFCE54", "#48CFAD","#33CCFF"),
add = c("jitter"))+stat_compare_means()
print(p)
dev.off()

###immune cell
library(ggplot2)
library(ggpubr)
setwd("D:\\11.BMAP\\12.CIBERSORT")
mydata <- read.table("BRCACancerCellProportions.txt",header=T,sep='\t',row.names=1,check.names=F)
mygroup <- read.table("BRCA_DNAm_group.txt",header=T,sep='\t',row.names=1,check.names=F)

data1 <- mydata[which(mydata[,23]<0.05),1:22]
data2 <- data.frame(data1,group=mygroup[rownames(data1),4])
write.table(data2,"cell_group.txt",quote=F,sep='\t')

data2$group <- factor(data2$group,level=c("High","Median","Low"))

setwd("D:\\11.BMAP\\12.CIBERSORT\\singleCell")
pdf("Neutrophils.pdf",width=7,height=7)
p <- ggplot(data2, aes(x=group, y=Neutrophils,fill=group)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill="white")+ #绘制箱线图
  scale_fill_manual(values = c("#E13220", "#DDAC2A", "#519DCC"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
	  legend.position="none")+xlab("")+stat_compare_means(label.x.npc="middle")
print(p)
dev.off()

###
library(reshape2)
mydata <- melt(data2,id.vars="group",variable.name="cell",value.name="value")
pdf("Cellboxplot.pdf",height=7,width=15)
ggplot(mydata,aes(x=cell,y=value,fill=group,col=group))+ 
  coord_cartesian(ylim = c(-0.05,0.8))+
  stat_compare_means(aes(group=group), label = "p.signif",method = "kruskal.test",
                     label.y = 0.7,hide.ns = TRUE)+#增加p值
  scale_color_manual(values=c("black","black","black"))+
  scale_fill_manual(values=c("#E13220", "#DDAC2A", "#519DCC"))+
  geom_boxplot(outlier.size=0.8,outlier.alpha=0.8)+
  theme_classic() +
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size = 16, 
                                   color = "black", face = "plain", 
                                   vjust = 1, hjust = 1,
                                   angle = 60))+
  theme(legend.spacing.x=unit(1.2,"cm"),
        axis.line = element_line(colour = "black",size=1))+ 
  theme(axis.ticks.x = element_blank())
dev.off()