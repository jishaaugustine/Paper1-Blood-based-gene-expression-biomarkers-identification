library(GEOquery)
library(affy)
library(limma)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(hgu133a.db)
library(hgu133acdf)
library(hgu133a2.db)
library(hgu133a2cdf)
library(hgu133plus2.db)
library(Biobase)
library(limma)
library(illuminaio)
library(illuminaHumanv3.db)
library(illuminaHumanv4.db)
library(metaMA)

setwd("~/Desktop/imple/microarray_integration_new")

#Read normalized data from GEO
library(preprocessCore)
##Read data
load("GSE72267.RData")
#boxplot(data.frame(exprs(GSE72267),main="GSE72267",outline=FALSE))
load("GSE6613.RData")
#boxplot(data.frame(exprs(GSE6613),main="GSE6613",outline=FALSE))
load("GSE99039.RData")
#boxplot(data.frame(exprs(GSE99039),main="GSE99039",outline=FALSE))
load("GSE57475.RData")
#boxplot(data.frame(exprs(GSE57475),main="GSE57475",outline=FALSE))
#Feature selection Algorithm
library(limma)
library(gtools)
feature_limma=function(data,label)
{
  design=model.matrix(~0+label)
  colnames(design)[1] <- 'control'
  colnames(design)[2] <- 'diseased'
  fit<-lmFit(as.data.frame(data),design)
  cont.matrix <- makeContrasts(diseased-control, levels=design)
  fit1 <- contrasts.fit(fit, cont.matrix)
  fit2 <-eBayes(fit1)
  top=topTable(fit2 , number=Inf, adjust.method = "BH", coef=1,sort.by = "p")
  top <- top[top$P.Val <= 0.05,]
  top$FC <- logratio2foldchange(top$logFC, base = 2)
  #top <- top[abs(top$FC) >= 1.2,]
  return(top)
}

# orglist=list(GSE72267,GSE99039,GSE6613,GSE57475)
# dataList=list(exprs(GSE72267),exprs(GSE99039),exprs(GSE6613),exprs(GSE57475))
# x <- list(hgu133a2SYMBOL,hgu133plus2SYMBOL, hgu133aSYMBOL, illuminaHumanv3SYMBOL)
orglist=list(GSE99039,GSE6613,GSE57475)
dataList=list(exprs(GSE99039),exprs(GSE6613),exprs(GSE57475))
x <- list(hgu133plus2SYMBOL, hgu133aSYMBOL,illuminaHumanv3SYMBOL)
# orglist=list(GSE99039,GSE6613)
# dataList=list(exprs(GSE99039),exprs(GSE6613))
# x <- list(hgu133plus2SYMBOL, hgu133aSYMBOL)


mapped_genes <- lapply(x,mappedkeys)
sub <- function(x, m) as.list(x[m])
linklist <- mapply(sub,x,mapped_genes)

probe2unigene<-function(expset,dataset,link)
{
  #construction of the map probe->unigene
  probes=rownames(expset)
  unigene=link[probes]
  names(unigene)<-probes
  probe_unigene=unigene
}

unigene2probe<-function(map)
{
  unigene_probe=reverseSplit(map)
}
convert2GeneID<-function(listStudies,datalist,linklist)
{
  listStudies = dataList
  datalist = orglist
  linklist = linklist
  if (!(class(listStudies) %in% c("list")))
  {
    stop("listStudies must be a list")
  }
  foo<-function(x,y,z) unigene2probe(probe2unigene(x,y,z))
  a=unigene2probe(probe2unigene(listStudies[[1]],datalist[[1]],linklist[[1]]))
  b=unigene2probe(probe2unigene(listStudies[[2]],datalist[[2]],linklist[[2]]))
  c=unigene2probe(probe2unigene(listStudies[[3]],datalist[[3]],linklist[[3]]))
  conv_unigene=mapply(foo,listStudies,datalist,linklist)
  id=lapply(conv_unigene,names)
  inter=Reduce(intersect,id)
  if(length(inter)<=0){stop("no common genes")}
  print(paste(length(inter),"genes in common"))
  lr=list()
  esets=lapply(1:length(listStudies),FUN=function(i)
  {
    l=lapply(conv_unigene[[i]][inter],FUN=function(x) listStudies[[i]][x,,drop=FALSE])
    
    esetsgr=t(sapply(l,FUN=function(y)
      if(is.null(dim(y))){y}
      else{y[which.max(rowMeans(y,na.rm=TRUE)),]}))
    esetsgr
  })
  
  probs=lapply(1:length(listStudies),FUN=function(i)
  {
    l=lapply(conv_unigene[[i]][inter],FUN=function(x) listStudies[[i]][x,,drop=FALSE])
    
    prob=lapply(l,FUN=function(y)
      if(is.null(dim(y))){y}
      else{row.names(y)[which.max(rowMeans(y, na.rm=TRUE))]})
    prob
  })
  
  return(list(esets=esets,conv.unigene=probs))
}

conv=convert2GeneID(dataList,orglist,linklist)
esets=conv$esets
probeList=conv$conv.unigene
###################################################################################
#Filter batches in each dataset which has less than 10 samples in each batch
phenoSub=list()
pheno=lapply(orglist,FUN=function(x) pData(x))
c=list()

# phenoSub[[1]]=pheno[[1]][pheno[[1]]$`diagnosis:ch1` %in% list("Parkinson's disease","Healthy"), c("geo_accession","diagnosis:ch1")]
# c[[1]]=as.factor(as.numeric(phenoSub[[1]]["diagnosis:ch1"]=="Parkinson's disease"))
# 
# phenoSub[[2]]=pheno[[2]][pheno[[2]]$`disease label:ch1` %in% list("IPD","CONTROL"),c("geo_accession","disease label:ch1","batch:ch1") ]
# c[[2]]=as.factor(as.numeric(phenoSub[[2]]["disease label:ch1"]=="IPD"))
# 
# phenoSub[[3]]=pheno[[3]][pheno[[3]]$`characteristics_ch1` %in% list("Parkinson's disease","healthy control"), c("geo_accession","characteristics_ch1")]
# c[[3]]=as.factor(as.numeric(phenoSub[[3]]["characteristics_ch1"]=="Parkinson's disease"))
# 
# phenoSub[[4]]=pheno[[4]][pheno[[4]]$`description` %in% list("PD","HC"),c("geo_accession","description") ]
# c[[4]]=as.factor(as.numeric(phenoSub[[4]]["description"]=="PD"))
phenoSub[[1]]=pheno[[1]][pheno[[1]]$`disease label:ch1` %in% list("IPD","CONTROL"),c("geo_accession","disease label:ch1","batch:ch1") ]
c[[1]]=as.factor(as.numeric(phenoSub[[1]]["disease label:ch1"]=="IPD"))
#
phenoSub[[2]]=pheno[[2]][pheno[[2]]$`characteristics_ch1` %in% list("Parkinson's disease","healthy control"), c("geo_accession","characteristics_ch1")]
c[[2]]=as.factor(as.numeric(phenoSub[[2]]["characteristics_ch1"]=="Parkinson's disease"))

phenoSub[[3]]=pheno[[3]][pheno[[3]]$`description` %in% list("PD","HC"),c("geo_accession","description") ]
c[[3]]=as.factor(as.numeric(phenoSub[[3]]["description"]=="PD"))

subfun <- function(x,y) x[,match(row.names(y),colnames(x))]
esets_new <- mapply(subfun,esets,phenoSub)

############################################
esetsT=list()
esetsT=lapply(esets_new,FUN=function(x) as.data.frame(t(x)))
#d1=feature_limma(as.data.frame(t(esetsT[[2]])),c[[2]])
# esetsT[[1]]$batch=as.factor(1)
# esetsT[[2]]$batch=as.factor(2)
# esetsT[[3]]$batch=as.factor(3)
# esetsT[[4]]$batch=as.factor(4)
# esetsT[[1]]$label=c[[1]]
# esetsT[[2]]$label=c[[2]]
# esetsT[[3]]$label=c[[3]]
# esetsT[[4]]$label=c[[4]]
# 
# mergeEsets=rbind(esetsT[[1]],esetsT[[2]],esetsT[[3]],esetsT[[4]])
esetsT[[1]]$batch=as.factor(1)
esetsT[[2]]$batch=as.factor(2)
esetsT[[3]]$batch=as.factor(3)
esetsT[[1]]$label=c[[1]]
esetsT[[2]]$label=c[[2]]
esetsT[[3]]$label=c[[3]]

#mergeEsets=rbind(esetsT[[1]],esetsT[[2]],esetsT[[3]],esetsT[[4]])
mergeEsets=rbind(esetsT[[1]],esetsT[[2]],esetsT[[3]])
#mergeEsets=rbind(esetsT[[1]],esetsT[[2]])

mergeSub=mergeEsets[mergeEsets$label==1,]
set.seed(123)
index <- sample(1:nrow(mergeSub), 30)
mergeSub=mergeSub[index,]
mergeSub1=mergeSub[,-which(names(mergeSub) %in% c("label","batch"))]

#Clustering to check if batch normalization is required
dist_mat <- dist(mergeSub1, method = 'euclidean')
#dist_mat <- as.dist(1-cor(t(mergeSub1), method="pearson"))
hclust_avg <- hclust(dist_mat, method = 'average')
dend=as.dendrogram(hclust_avg)

library("dendextend")
labels(dend)=mergeSub$batch[hclust_avg$order]
plot(dend)

#Batch normalization using COMBAT

library(sva)
#mergeData=cbind(esets_new[[1]],esets_new[[2]],esets_new[[3]],esets_new[[4]])
# #Z transform
#mergeDataZ=rbind(scale(t(esets_new[[1]])),scale(t(esets_new[[2]])),scale(t(esets_new[[3]])),scale(t(esets_new[[4]])))
##log2normalize data
mergeData=cbind(esets_new[[1]],esets_new[[2]],esets_new[[3]])
# #Z transform
mergeDataZ=rbind(scale(t(esets_new[[1]])),scale(t(esets_new[[2]])),scale(t(esets_new[[3]])))
#mergeData=cbind(esets_new[[1]],esets_new[[2]])
#Z transform
#mergeDataZ=rbind(scale(t(esets_new[[1]])),scale(t(esets_new[[2]])))

phenoDta=mergeEsets[,c("batch","label")] 
# batch=phenoDta$batch
# modcombat = model.matrix(~1, data=phenoDta)
# mergeDataC = ComBat(dat=mergeData, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# 


#Significance analysis directly on the adjusted data 
# mod = model.matrix(~as.factor(label), data=phenoDta)
# mod0=model.matrix(~1,data=phenoDta)
# pValuesComBat = f.pvalue(combat_edata,mod,mod0)
# qValuesComBat = p.adjust(pValuesComBat,method="BH")

#Clustering to check if batch normalization is required
#dist_mat <- dist(mergeSub1, method = 'euclidean')
combat_edataT=as.data.frame(mergeDataZ)
#combat_edataT=as.data.frame(t(mergeDataC))
combat_edataT$batch=phenoDta$batch
combat_edataT$label=phenoDta$label

mergeSub=combat_edataT[combat_edataT$label==1,]
set.seed(123)
index <- sample(1:nrow(mergeSub), 30)
mergeSub=mergeSub[index,]
mergeSub1=mergeSub[,-which(names(mergeSub) %in% c("label","batch"))]

#Clustering after batch normalization
dist_mat <- dist(mergeSub1, method = 'euclidean')
#dist_mat <- as.dist(1-cor(t(mergeSub1), method="pearson"))
hclust_avg <- hclust(dist_mat, method = 'average')
dend=as.dendrogram(hclust_avg)

library("dendextend")
labels(dend)=mergeSub$batch[hclust_avg$order]
plot(dend)

mergeEsetsZ=t(mergeDataZ)
#mergeEsetsC=mergeDataC
#save(mergeEsetsC,phenoDta,probeList,file="CombindCNolimma_result.RData")
#save(mergeEsetsZ,phenoDta,probeList,file="CombindZNolimma2data_result.RData")
#save(mergeEsetsC,file="dataC_nolimma.csv")
#save(mergeEsetsZ,file="dataZ_nolimma.csv")
label1=phenoDta$label
sigZ=feature_limma(mergeEsetsZ,label1)

mergeEsetsZT=as.data.frame(t(mergeEsetsZ))
#combat_edata=mergeEsetsZT[,colnames(mergeEsetsZT) %in% row.names(sigZ)]
combat_edata=mergeEsetsZT
# #################################################################
save(sigZ,file="limmaZ_3data.RData")
# 
save(combat_edata,phenoDta,probeList,file="CombindZ_4datafull.RData")
# 
#######################################################
#Validation GSE72267
load("GSE72267.RData")
#Filter low intensity probes in each dataset
phenoData=pData(GSE72267)
#phenoSub = phenoData[phenoData$`characteristics_ch1` %in% list("status: AD","status: CTL"),c("geo_accession","characteristics_ch1")]
phenoSub=phenoData[phenoData$`diagnosis:ch1` %in% list("Parkinson's disease","Healthy"), c("geo_accession","diagnosis:ch1")]

esets=exprs(GSE72267)[,match(row.names(phenoSub),colnames(exprs(GSE72267)))]
c1=as.factor(as.numeric(phenoSub["diagnosis:ch1"]=="Parkinson's disease"))

#c1=as.factor(as.numeric(phenoSub["characteristics_ch1"]=="status: AD"))
##Skipped batch normalization###

#GENE SYMBOL
#x <- illuminaHumanv3SYMBOL
x<-hgu133a2SYMBOL
mapped_genes <- mappedkeys(x)
link <- as.list(x[mapped_genes])
probe2unigene<-function(expset)
{
  probes=rownames(expset)
  unigene=link[probes]
  names(unigene)<-probes
  probe_unigene=unigene
}
unigene2probe<-function(map)
{
  unigene_probe=reverseSplit(map)
}


lstSig=unigene2probe(probe2unigene(esets))
filtered_data=esets

l=lapply(lstSig,FUN=function(x) filtered_data[x,,drop=FALSE])

filtered_symbol=t(sapply(l,FUN=function(y)
  if(is.null(dim(y))){y}
  else{y[which.max(rowMeans(y,na.rm=TRUE)),]}))


filtered_probe=lapply(l,FUN=function(y)
  if(is.null(dim(y))){y}
  else{row.names(y)[which.max(rowMeans(y, na.rm=TRUE))]})

filtered_probe1=stack(filtered_probe)

validationset=scale(t(filtered_symbol))
validationset_label=c1
save(validationset,validationset_label,file="GSE72267_validation.RData")
#################################################################################
#GSE51799 validataion
library(biomaRt)
library("org.Hs.eg.db")

#Read normalized data
GSE51799<-read.csv("GSE51799.txt", header=T, row.names=1,sep="\t")

mart = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
ann <- getBM(c("hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id"), "ensembl_gene_id", rownames(GSE51799), mart)
ann.fixed <- ann[match(rownames(GSE51799), ann[,8]),c("hgnc_symbol","ensembl_gene_id")]
ann.fixed=ann.fixed[!(is.na(ann.fixed$hgnc_symbol) | ann.fixed$hgnc_symbol==""), ]
GSE51799=GSE51799[ann.fixed$ensembl_gene_id,]
rownames(GSE51799)=ann.fixed$hgnc_symbol

#Read pData
eset_GSE51799 <- getGEO("GSE51799", GSEMatrix=TRUE, getGPL=TRUE)
if(length(eset_GSE51799) > 1)
  idx <- grep("GPL10999", attr(eset_GSE51799, "names")) else idx <- 1
eset_GSE51799_1 <- eset_GSE51799[[idx]]
if(length(eset_GSE51799) > 1)
  idx <- grep("GPL11154", attr(eset_GSE51799, "names")) else idx <- 1
eset_GSE51799_2 <- eset_GSE51799[[idx]]

GSE51799_pheno=rbind(pData(eset_GSE51799_1),pData(eset_GSE51799_2))
GSE51799_phenoData=GSE51799_pheno[match(colnames(GSE51799),GSE51799_pheno$description),c("geo_accession","source_name_ch1","description")]
colnames(GSE51799)=GSE51799_phenoData$geo_accession
c1=as.factor(as.numeric(GSE51799_phenoData["source_name_ch1"]=="Whole blood, carrier"))
validationset=as.data.frame(t(GSE51799))
validationset_label=c1

save(validationset,validationset_label,file="GSE51799_validation.RData")

#################################################################################################
#GSE85426 validation
#Read normalized data
GSE85426_data<- read.csv("GSE85426.txt", header=T, row.names=1,sep="\t")
#Read pData
eset_GSE85426 <- getGEO("GSE85426", GSEMatrix=TRUE, getGPL=TRUE)
if(length(eset_GSE85426) > 1)
  idx <- grep("GPL14550", attr(eset_GSE85426, "names")) else idx <- 1
eset_GSE85426 <- eset_GSE85426[[idx]]
sample=as.data.frame(colnames(GSE85426_data))
sampleSplit <- data.frame(do.call('rbind', strsplit(as.character(sample$`colnames(GSE85426_data)`),'.',fixed=TRUE)))
sampleName=sampleSplit$X1
colnames(GSE85426_data)=sampleName
eset_GSE85426_df=as.data.frame(pData(eset_GSE85426))
GSE85426_phenoData=eset_GSE85426_df[match(sampleName,eset_GSE85426_df$description.1),c("geo_accession","description.1","characteristics_ch1")]

library(tidyr)
library(reshape)
gene=as.data.frame(rownames(GSE85426_data))
#gene %>% separate(gene$`rownames(GSE85426)`,c("ID", "name"),"\\|",extra='drop')
foo <- data.frame(do.call('rbind', strsplit(as.character(gene[,1]),'|',fixed=TRUE)))
foo1 <- data.frame(do.call('rbind', strsplit(as.character(foo$X1),',',fixed=TRUE)))
# # df2=with(gene, colsplit(gene$`rownames(GSE85426)`, split = "\\|",names=c("col1","col2")))
# # df3=with(df2, colsplit(df2$col1, split = "\\,",names=c("col11","col12","col13","col14","col15","col16","col17","col18","col19","col20","col21","col22","col23","col24","col25")))
geneName=foo1$X1
row.names(GSE85426_data)=geneName
# annotLookup2.fixed <- annotLookup2[match(geneName, annotLookup2$agilent_sureprint_g3_ge_8x60k_v2),c("agilent_sureprint_g3_ge_8x60k_v2","external_gene_name")]
# 
# 
# library(GEOquery)
# anno=getGEO("GPL13607")
# annodata=Table(anno)
# annodata.fix=annodata[match(df3$NA..15, annodata$ProbeName),c("ProbeName","GeneName")]

library(hgug4110b.db)
x<-hgug4110bSYMBOL
mapped_genes <- mappedkeys(x)
linklist <- as.list(x[mapped_genes])

probe2unigene<-function(expset)
{
  probes=rownames(expset)
  unigene=linklist[probes]
  names(unigene)<-probes
  probe_unigene=unigene
}
unigene2probe<-function(map)
{
  unigene_probe=reverseSplit(map)
}
lstSig=unigene2probe(probe2unigene(GSE85426_data))
filtered_data=GSE85426_data

l=lapply(lstSig,FUN=function(x) filtered_data[x,,drop=FALSE])

filtered_symbol=t(sapply(l,FUN=function(y)
  if(is.null(dim(y))){y}
  else{y[which.max(rowMeans(y,na.rm=TRUE)),]}))


filtered_probe=lapply(l,FUN=function(y)
  if(is.null(dim(y))){y}
  else{row.names(y)[which.max(rowMeans(y, na.rm=TRUE))]})

filtered_probe1=stack(filtered_probe)

filtered_symbol=filtered_symbol[,match(GSE85426_phenoData$description.1,colnames(filtered_symbol))]
colnames(filtered_symbol)=GSE85426_phenoData$geo_accession
c1=as.factor(as.numeric(GSE85426_phenoData["characteristics_ch1"]=="diagnosis: Probable Alzheimer's Disease"))

validationset=as.data.frame(t(filtered_symbol))
validationset_label=c1
save(validationset,validationset_label,file="GSE85426_validation.RData")


######################################################################################
library(illuminaio)
library(illuminaHumanv3.db)
library(limma)
data1 = getGEO('GSE63060',GSEMatrix =TRUE, getGPL=TRUE, AnnotGPL=TRUE)
phenoData=pData(data1[[1]])

phenoSub = phenoData[phenoData$`status:ch1` %in% list("AD","CTL"),c("geo_accession","status:ch1")]
esets=data1[[1]][,row.names(phenoSub)]
c1=as.factor(as.numeric(phenoSub["status:ch1"]=="AD"))

#GENE SYMBOL
x <- illuminaHumanv3SYMBOL
#x<-org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
link <- as.list(x[mapped_genes])
probe2unigene<-function(expset)
{
  #construction of the map probe->unigene
  probes=row.names(expset)
  #gene_id=fData(expset)[probes,"ENTREZ_GENE_ID"]
  unigene=link[probes] #Illumina
  #unigene=link[gene_id]
  names(unigene)<-probes
  probe_unigene=unigene
}
unigene2probe<-function(map)
{
  unigene_probe=reverseSplit(map)
}


lstSig=unigene2probe(probe2unigene(esets))
filtered_data=exprs(esets)

l=lapply(lstSig,FUN=function(x) filtered_data[x,,drop=FALSE])

filtered_symbol=t(sapply(l,FUN=function(y)
  if(is.null(dim(y))){y}
  else{y[which.max(rowMeans(y,na.rm=TRUE)),]}))


filtered_probe=lapply(l,FUN=function(y)
  if(is.null(dim(y))){y}
  else{row.names(y)[which.max(rowMeans(y, na.rm=TRUE))]})

filtered_probe1=stack(filtered_probe)

validationset=scale(t(filtered_symbol))
validationset_label=c1
save(validationset,validationset_label,file="GSE63060_validation.RData")

