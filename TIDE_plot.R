setwd("D:\\11.BMAP\\16.TIDE")
mydata <- read.table("TIDE_groupResult.txt",header=T,sep='\t',row.names=1,check.names=F)

library(ggplot2)
library(ggpubr)
mydata$group <- factor(mydata$group,level=c("High","Median","Low"))
compar<-list(c("High","Median"),c("High","Low"),c("Median","Low"))
pdf("TIDEviolin_group.pdf",height=7,width=7)

ggplot(data = mydata,aes(x = group, #分组列名
                      y = TIDE, #连续变量列名
                      fill = group))+ #按分组填充颜色
  scale_fill_manual(values = c("#E13220", "#DDAC2A", "#519DCC")) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  geom_point(shape = 21, size=1.5, # 点的性状和大小
             position = position_jitterdodge(jitter.width=1), # 让点散开
             color="black", alpha = 1) +
  theme_classic() + 
  ylab("TIDE") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) +
stat_compare_means(comparisons = compar)+stat_compare_means()
dev.off()