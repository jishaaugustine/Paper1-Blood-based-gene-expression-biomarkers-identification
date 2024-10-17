library(org.Hs.eg.db)
library(AnnotationHub)
library(clusterProfiler)
library(AnnotationDbi)
setwd("~/Desktop/imple/microarray_integration_new")
load("Regulation.RData")
###############################
library(data.table)
library(stringr)
data1 <- row.names(sigPositive)

gene.df <- bitr(data1, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

ALL_positive <- enrichGO(gene         = gene.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod  = "BH",
                        readable = TRUE
)

ego_all <- setReadable(ALL_positive, OrgDb = org.Hs.eg.db)
#dotplot(ego2, showCategory=5)
library(DOSE)
p1=dotplot(ego_all, showCategory=15)

#BP,MF,CC
BP_positive <- enrichGO(gene         = gene.df$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP" ,
                 pAdjustMethod  = "BH",
                 readable = TRUE
                )

MF_positive <- enrichGO(gene         = gene.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF" ,
                        pAdjustMethod  = "BH",
                        readable = TRUE
)

CC_positive <- enrichGO(gene         = gene.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC" ,
                        pAdjustMethod  = "BH",
                        readable = TRUE
)
dim(GO_positive)


ego1 <- setReadable(BP_positive, OrgDb = org.Hs.eg.db)
ego2 <- setReadable(MF_positive, OrgDb = org.Hs.eg.db)

#dotplot(ego2, showCategory=5)
p1=dotplot(ego2, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))

# save(GO_positive,file="GO_positive.RData")
# write.table(GO_positive[1:1353],file="GO_positive.csv")
#############################################################
data1 <- row.names(sigNegative)

gene.df <- bitr(data1, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

BP_negative <- enrichGO(gene         = gene.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP" ,
                        pAdjustMethod  = "BH",
                        readable = TRUE
)

dim(GO_negative)

MF_negative <- enrichGO(gene         = gene.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF" ,
                        pAdjustMethod  = "BH",
                        readable = TRUE
)

CC_negative <- enrichGO(gene         = gene.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC" ,
                        pAdjustMethod  = "BH",
                        readable = TRUE
)
# save(GO_negative,file="GO_negative.RData")
# write.table(GO_negative[1:29],file="GO_negative.csv")


ego3 <- setReadable(BP_negative, OrgDb = org.Hs.eg.db)
p1 = dotplot(ego3, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
#dotplot(ego3, showCategory=5,label_format=10)
p1
library(enrichplot)

########################################################
#KEGG
data1 <- row.names(sigPositive)

gene.df <- bitr(data1, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

KEGG_positive<-enrichKEGG(gene.df$ENTREZID, 
                organism = "hsa", 
                keyType = "kegg"
)

ego3 <- setReadable(KEGG_positive, OrgDb = org.Hs.eg.db)
library(enrichplot)
dotplot(KEGG_positive, showCategory=15)
p1 = barplot(KEGG_positive, showCategory=15) + scale_y_discrete(labels=function(x) str_wrap(x, width=50))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))

dim(KEGG_positive)
save(KEGG_positive,file="KEGG_positive.RData")
write.table(KEGG_positive[1:101],file="KEGG_positive.csv")

###########################################################################
data1 <- row.names(sigNegative)

gene.df <- bitr(data1, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

KEGG_negative<-enrichKEGG(gene.df$ENTREZID, 
                          organism = "hsa", 
                          keyType = "kegg"
)

dim(KEGG_negative)
save(KEGG_negative,file="KEGG_negative.RData")
write.table(KEGG_negative[1:101],file="KEGG_negative.csv")

########################
#PLOT
library(png)
library(eps)
library(grid)
library(gridExtra)

plot1 <- readPNG('BP_positive.png')
plot2 <- readPNG('BP_negative.png')
plot3 <- readPNG('MF_positive.png')
plot4 <- readPNG('MF_negative.png')
plot5 <- readPNG('CC_positive.png')
plot6 <- readPNG('CC_negative.png')
plot7 <- readPNG('KEGG_positive1.png')

library(ggpubr)
figure <- ggarrange(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot5),rasterGrob(plot6),rasterGrob(plot7),
                    labels = c("a", "b", "c","d","e","f"),
                    ncol = 2, nrow = 4)


library("cowplot")
p1=plot_grid(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot4),rasterGrob(plot7),
          labels=c("a", "b", "c","d"),ncol = 2, nrow = 2)
#ggsave(file="test.eps",limitsize=TRUE)
ggsave("Fig2.tiff", plot = p1, width = 17, height = 20, units = "cm")
library(ggplot2)
png("test.png",res=300)
plot_grid(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),rasterGrob(plot4),rasterGrob(plot7),
          labels="AUTO",label_size = 5,ncol = 2, nrow = 3)
dev.off()


ego1 <- setReadable(BP_positive, OrgDb = org.Hs.eg.db)
ego2 <- setReadable(BP_negative, OrgDb = org.Hs.eg.db)
ego3 <- setReadable(MF_positive, OrgDb = org.Hs.eg.db)
ego4 <- setReadable(MF_negative, OrgDb = org.Hs.eg.db)
ego5 <- setReadable(CC_positive, OrgDb = org.Hs.eg.db)
ego6 <- setReadable(CC_negative, OrgDb = org.Hs.eg.db)


#dotplot(ego2, showCategory=5)
png("BP_positive.png",res=100)
dotplot(ego1, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
dev.off()
png("BP_negative.png",res=100)
dotplot(ego2, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
dev.off()
png("MF_positive.png",res=100)
dotplot(ego3, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
dev.off()
png("MF_negative.png",res=100)
dotplot(ego4, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
dev.off()
p5=dotplot(ego5, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
p6=dotplot(ego6, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
png("KEGG_positve.png",res=100)
barplot(KEGG_positive, showCategory=15) + scale_y_discrete(labels=function(x) str_wrap(x, width=50))+theme(plot.margin=unit(c(.25,0,.25,.25), "cm"))
dev.off()
png("test.png",res=100)
plot_grid(p1,p2,p3,p4,p7,labels="AUTO",label_size = 5,ncol = 2, nrow = 3)
dev.off()
