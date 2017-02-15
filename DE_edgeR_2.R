# Created on 09/12/2016 to process RNA-seq data 
# Author :JX
# updated on 10/12/2016 to remove the genes annotated by chrM

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
# workdir ="/Users/jianjinxu/Google Drive/projects/RNA-seq/"
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
#write.table(out1[,c("Genename","metric")],file =paste(outdir,"DE_edgeR_all",Sys.Date(),".rnk",sep=""),quote=F,sep="\t",row.names=F)
tt =rbind(head(out1,3),tail(out1,3))
#write.csv(tt,file =paste(outdir,"DE_edgeR_rank_top_tail",Sys.Date(),".csv",sep=""),row.names=F)

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

## save top 100 genes 
n.top=gene.tot[t1,]
write.csv(n.top,file =paste(outdir,"n_top_",topn,"_",Sys.Date(),".csv",sep=""))

cb.top=gene.tot[t2,]
write.csv(cb.top,file =paste(outdir,"cb_top_",topn,"_",Sys.Date(),".csv",sep=""))


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
#write.csv(gene.both.high,file =paste(outdir,"both_high_top_",topn,"_",Sys.Date(),".csv",sep=""))

## see if those genes are DE



########### heatmap for the top common genes ############

gene.rpkm.1 = gene.both.high[,sampleDict]
### standardize 
gene.rpkm.1 = t(scale(t(gene.rpkm.1)))

library(gplots)
col.pan = colorpanel(100,"blue","white","red")
fname =paste("heatmap_top_",topn,"_common_",sep="")
print(fname)
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 1200, height = 600, units = "px", pointsize = 12) 
heatmap.2(gene.rpkm.1, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

dev.off()

##########################highly expressed vs. no expression##########################
# histogram of RPKM in one sample
tmp= gene.rpkm
dim(tmp)
sum.stat=apply(tmp,2,summary)

#write.csv(sum.stat,file =paste(outdir,"RPKM_summary_stat_",Sys.Date(),".csv",sep=""))
## no exp in fetus 0.1, high exp at adult 75% 5
## no exp in adult 0.04, high exp at fetus 75% 4
tt=apply(tmp,2,function(x) quantile(x,probs = c(0.05,0.95)))
print(tt)
aa=tmp[tmp[,1]<5,1]

tt=apply(tmp,2,function(x) quantile(x,probs = c(0.05,0.90)))
print(tt)
print(tt)
## define 5% as low cutoff
fname ="density_plot_RPKM_samples"
## 10-CB, 11-N
## define 50% as moderate
samp =c("10-CB","11-N")
cols =c ("red","green")
i=1
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 600, height = 500, units = "px", pointsize = 20) 
par(new = F)
xrange =c(0,20)
yrange =c(0,1)

for (i in  1:2){
  aa=tmp[tmp[,samp[i]]<20,samp[i]]
  # aa=tmp[,samp[i]]
  d=density(aa)
  plot(d, col=cols[i],main="",xlab = "", ylab = "", yaxt="n",xaxt="n",ylim=yrange,xlim=xrange,lwd=2)
  par(new = T)
}
title(xlab ="RPKM", ylab="Density")


axis(1)
axis(2)
legend("topright", lwd = 2,
       legend = samp,
       col = cols,
       lty = 1, box.col="white")

box()

dev.off()

### select genes###########
tt=apply(tmp,2,function(x) quantile(x,probs = c(0.3,0.5)))
print(tt)
cf =c(0.4,1)
d1= tmp[,ind.n]

low.0 = (rowSums((d1< cf[1])) == ncol(d1))
sum(low.0)
# 198 genes selected#####
high.0 = (rowSums((d1 >cf[2])) == ncol(d1))
sum(high.0)
# 2580
cf =c(0.6,1.3)
d1= tmp[,ind.cb]

low.1 = (rowSums((d1< cf[1])) == ncol(d1))
sum(low.1)
# 198 genes selected#####
high.1 = (rowSums((d1 >cf[2])) == ncol(d1))
sum(high.1)
## adult low and fetus high
(sum(low.0 & high.1))
out.0 =gene.tot[which(low.0 & high.1),]
(sum(low.1 & high.0))
out.1 =gene.tot[which(low.1 & high.0),]

out.2=rbind(out.0,out.1)
write.csv(out.2,file =paste(outdir,"high_low_list_better",Sys.Date(),".csv",sep=""))

####################### contamination analysis##################

# CD45,71 and 235 
# PTPRC, TFRC, ADCY4
# gene.tot[,]

markers =c("PTPRC", "TFRC", "GYPA")
check.0= gene.tot[which(gene.tot$Genename.x %in% markers),]
which(gene.inf[,"Genename"] =="ADCY4")


check.1=t(scale(t(check.0[,sampleDict])))
# check.1=t((t(check.0[,sampleDict])))
rownames(check.1) = check.0$Genename.x
library(gplots)
col.pan = colorpanel(100,"blue","white","red")
fname =paste("heatmap_contamination_",sep="")
print(fname)
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 700, height = 500, units = "px", pointsize = 20) 
heatmap.2(check.1, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

dev.off()


write.csv(check.0,file =paste(outdir,"check_marker_",Sys.Date(),".csv",sep=""))

write.csv(gene.tot,file =paste(outdir,"all_info_",Sys.Date(),".csv",sep=""))

