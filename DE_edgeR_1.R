# Created on 09/12/2016 to process RNA-seq data 
# Author :JX

# input: 
# test:data frame with all genes count, m*n, m is the number of genes, n is the total sample size
# class: group assignment for each sample
# CB is fetus
# install.packages("edgeR")
# source("https://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# biocLite("limma")
# install.packages("gplots")
rm(list =ls())
library(edgeR)
workdir = "C:/Users/jackie/Google Drive/projects/RNA-seq/"
datadir = paste(workdir,"data/",sep="")
figdir = paste(workdir,"figure/",sep="")
outdir = paste(workdir,"output/",sep="")
refdir = paste(workdir,"to_deliver/",sep="")
#refdir ="/rna-seq/raw/Project_GNA_11920_B01_EXS_RNA.2016-08-24/Results_Aug23_2016/to_deliver/"
#### gene annotation file ######
fn.gene="gene_annotation_info_1.csv"
gene.anno=read.csv(file=paste(refdir,fn.gene,sep=""))
colnames(gene.anno) = c("Geneid","Genename")
setwd(workdir)
fn ="fullTable"
raw =read.csv(file = paste(datadir,fn,".csv",sep=""))[,-1]

gene.inf = raw[,c(1,2)]
gene.inf= merge(x=gene.inf,y= gene.anno,by="Geneid", all=TRUE)

raw.1 =raw[,c(-1,-2)]
rownames(raw.1)= raw[,1]
sample.name = colnames(raw.1)

sampleDict = c("10-CB","11-N","13-N","14-N","15-N","16-N","17-N","18-N",
               "19-N","1-CB","20-N","21-CB-R","22-CB","23-CB","2-CB","3-CB",
               "4-CB","5-CB","7-CB","9-CB")
ind.cb= grep("CB",sampleDict)
ind.n =grep("N",sampleDict)
class =rep(0,length(sampleDict))
class[ind.cb] = 0
class[ind.n] = 1
################## DE ####################
## filtering
keep = rowSums(cpm(raw.1)>1) >=2
test = as.data.frame(raw.1[keep,])
gene.inf.1 = gene.inf[keep,]
#rownames(test) = gene.inf.1[,"Genename"]
# 12941 unique genenames, how can it be possible?
# from  57445 to only 12944
# create DGEList
edgeR.dgelist = DGEList(counts = test, group = factor(class))

######### filtering#########
# edgeR.dgelist$samples
edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
#edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM", logratioTrim = 0.6, sumTrim = 0.5)

edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
 tmp= edgeR.dgelist$samples
 rownames(tmp) =sampleDict
 # write.csv(tmp,file = paste(outdir,"library_size.csv",sep=""))
edgeR.test = exactTest(edgeR.dgelist)
edgeR.pvalues = edgeR.test$table$PValue
edgeR.test$table$adjP = p.adjust(edgeR.pvalues, method = "BH")
summary(edgeR.test$table$adjP)



# 10257
sum(edgeR.test$table$adjP < 0.05)
edgeR.test$table$Genename = gene.inf.1[,"Genename"]
head(edgeR.test$table,3)

edgeR.test$table$metric =-log10(edgeR.test$table$adjP)/sign(edgeR.test$table$logFC)
out = edgeR.test$table
out1=out[order(out$metric,decreasing = T),]
out1$ID =rownames(out1)
write.table(out1[,c("Genename","metric")],file =paste(outdir,"DE_edgeR_all",Sys.Date(),".rnk",sep=""),quote=F,sep="\t",row.names=F)
tt =rbind(head(out1,3),tail(out1,3))
write.csv(tt,file =paste(outdir,"DE_edgeR_rank_top_tail",Sys.Date(),".csv",sep=""),row.names=F)

############## heatmap of RPKM for top 30 DE genes#################
# transform the data into RPKM ,we also need to know the gene length
logcpm=cpm(edgeR.dgelist,prior.count = 2,log=TRUE)
edgeR.dgelist$genes$Length = gene.inf.1[,2]
gene.rpkm= rpkm(edgeR.dgelist)
rownames(gene.rpkm) = gene.inf.1[,"Genename"]

colnames(gene.rpkm) = sampleDict 
topn =30
o= order(edgeR.test$table$adjP,decreasing = FALSE)


####### #1 updated on 09/26/2016 to visualize top expressed genes in both groups##########
# order(gene.rpkm)
# gene.rpkm[,ind.cb]
ave.cb= rowMeans(gene.rpkm[,ind.cb])
ave.n= rowMeans(gene.rpkm[,ind.n])
summary(ave.cb)
summary(ave.n)
ord.cb=rank(-ave.cb)
ord.n = rank(-ave.n)
adult.fetus = ord.n - ord.cb 

### add DE analysis result there
gene.rpkm.ord= as.data.frame(cbind(gene.rpkm,ave.n, ave.cb,ord.n,ord.cb,adult.fetus))
gene.rpkm.ord$ID =gene.inf.1[,"Geneid"]
gene.rpkm.ord$Genename=rownames(gene.rpkm.ord)
head(out1,2)
head(gene.rpkm.ord,2)
gene.tot= merge(gene.rpkm.ord,out1,by="ID")

topn =100

t1=(gene.tot$ord.n %in% 1:topn)
t2=(gene.tot$ord.cb %in% 1:topn)
t3 = cbind(t1,t2)
t3.1 =t3[rowSums(t3)>0,]
library(limma)

t4 = vennCounts(t3.1)
t4
fname =paste("Venn",topn,sep="")
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 600, height = 600, units = "px", pointsize = 12) 
vennDiagram(t4, include = "both", 
            names = c("High in N", "High in CB"), 
            cex = 1,  circle.col=c("green","blue"),counts.col = "black")
title(paste("Top", topn,"genes" ,sep=""))
dev.off()

gene.both.high = gene.tot[rowSums(t3)>1,]
sum(gene.both.high$adjP <0.05)
# 400 are also DE
head(gene.both.high)
write.csv(gene.both.high,file =paste(outdir,"both_high_top_",topn,"_",Sys.Date(),".csv",sep=""))

## see if those genes are DE
gene.both.high



############ #2 highly different expressed #######

# length(unique(gene.rpkm.ord$Genename))
## N bigger than CB 6461

tmp.1 = gene.tot[gene.tot$adult.fetus <0,] 

ord.adjP =rank(tmp.1$adjP)
ord.logFC = rank(-tmp.1$logFC)
ord.rpkm =  rank(tmp.1$adult.fetus)
ord.all =(ord.adjP+ord.logFC+ord.rpkm)/3
tmp.2=cbind(tmp.1, ord.adjP,ord.logFC,ord.rpkm,ord.all)
tmp.3 = tmp.2[(ord.adjP <=3000) & (ord.logFC <= 3000) & (ord.rpkm <=3000),]

tmp.4=tmp.3[order(tmp.3$ord.all),][1:100,]

### cluster samples by using all genes######
gene.rpkm.0 = t(scale(t(gene.rpkm)))
sampleDists = dist( t( gene.rpkm.0 ) )
# sampleDists
#install.packages("pheatmap")
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )

fname ="cluster_samples"
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 950, height = 500, units = "px", pointsize = 12) 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
########### heatmap for top 30 genes######
gene.rpkm.1 = gene.rpkm[o[1:topn],]
### standardize 
gene.rpkm.1 = t(scale(t(gene.rpkm.1)))

library(gplots)
col.pan = colorpanel(100,"blue","white","red")
fname ="heatmap"
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
     width = 1200, height = 600, units = "px", pointsize = 12) 
heatmap.2(gene.rpkm.1, col=col.pan, Rowv=TRUE, scale="none",
           trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

dev.off()
################### compare to DESeq2 result
res.1=edgeR.test$table
res.1 =res.1[o,]
de.edgeR =res.1[which(res.1$adjP < 0.05),]
write.csv(de.edgeR,file =paste(outdir,"DE_gene_list_edgeR.csv",sep=""))
res.2=res.1[1:topn,]["Genename"]

deseq.gene="DESeq2_results_CB_vs_N.txt"
res.0=read.delim(file=paste(refdir,deseq.gene,sep=""),sep="/t",header = T)
# 15342 DE
# 28658
deseq.0= res.0[which(!is.na(res.0$padj)),]

o.1 =order(deseq.0$padj,decreasing = FALSE)
deseq.0 =deseq.0[o.1,]
de.deseq =deseq.0[which(deseq.0$padj < 0.05),]
nrow(de.deseq)
# 15342
deseq.1=de.deseq[1:topn,]["Gene_name"]
sum(!is.na(res.0$padj))
# 28658 genes enter the study
#########
sum(de.edgeR$Genename %in% de.deseq$Gene_name)
## 10151 all detected #######
topn=1000
sum(de.edgeR$Genename[1:topn] %in% de.deseq$Gene_name[1:topn])
# 849 genes 