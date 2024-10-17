#setwd("~/Desktop/imple/microarray_integration_new")
load("GBM10_150.RData")
auc=as.data.frame(l[[2]])
gbm=c()
for (i in seq(1,30,2))
{
  gbm=rbind(gbm,auc[4,i])
}
gbm=as.data.frame(gbm)
gbm$gene_size=c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)
#Connected scatter plot

library(tidyr)
library(ggplot2)

ggplot(gbm, aes(x=gene_size, y=V1))+ geom_point()+ labs(y="AUC", x = "Gene size")

g=gbm %>% 
  gather(key, value, -gene_size) %>% 
  ggplot(aes(gene_size, value)) + geom_point(aes(color = key))  +
  scale_colour_manual(values=c("red", "blue", "black","yellow"))

g+labs(y="Performance",x="Bag size",color="Classifier" )

