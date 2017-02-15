############## explore the  raw data #############
rm(list =ls())
workdir = "C:/Users/jianjin/Google Drive/paired_RNA_seq/"
# workdir ="/Users/jianjinxu/Google Drive/paired_RNA_seq/"
codedir = paste(workdir,"code/",sep="")
datadir = paste(workdir,"data/",sep="")
outdir = paste(workdir,"output/",sep="")
figdir = paste(outdir,"figures/",sep="")
library(SARTools)


## simulate the raw data
library(RColorBrewer)
col = c(brewer.pal(12, "Paired"))
# display.brewer.pal(12, "Paired")

workdir = "C:/Users/jianjin/Google Drive/projects/RNA-seq/"
workdir ="/Users/jianjinxu/Google Drive/projects/RNA-seq/"
datadir = paste(workdir,"data/",sep="")
figdir = paste(workdir,"figure/presentation/",sep="")
outdir = paste(workdir,"output/",sep="")
refdir = paste(workdir,"to_deliver/",sep="")
codedir = paste(workdir,"code/",sep="")


#### gene annotation file ######
fn.gene="gene_annotation_info.txt"
gene.anno=read.delim2(file=paste(refdir,fn.gene,sep=""),sep="\t")

# gene.anno[gene.anno$Chromosome=="chrM",]
setwd(workdir)
fn ="fullTable"
raw =read.csv(file = paste(datadir,fn,".csv",sep=""))[,c(-1)]
colnames(raw)[1] ="Ensembl_Gene_ID"
colSums(raw[,c(-1,-2)])


gene.inf= merge(x=raw,y= gene.anno,by="Ensembl_Gene_ID", all=TRUE)

sampleDict = c("10-CB","11-N","13-N","14-N","15-N","16-N","17-N","18-N",
               "19-N","1-CB","20-N","21-CB-R","22-CB","23-CB","2-CB","3-CB",
               "4-CB","5-CB","7-CB","9-CB")



############# chrM ######################


## normaliztion methods

test =raw[,-c(1,2)]
colnames(test) =sampleDict


sampleDict1 = c("20-N","11-N","13-N","14-N","15-N","16-N","17-N","18-N",
                "19-N","1-CB","21-CB-R","22-CB","23-CB","2-CB","3-CB",
                "4-CB","5-CB","7-CB","9-CB","10-CB")

X=test[,sampleDict1]

(condition = unlist(lapply(colnames(X), function(x){unlist(strsplit(x,split = "-"))[2]})))
# condition = ifelse(condition =="N", "pre","post")

# the level statement deteremined the reference level is CB
group = factor(condition, levels=c("CB","N"))
levels(group)

# group = factor(condition, levels=c("pre","post"))

## filtering

cri = 0.5
keep = rowSums(cpm(X)>= cri) >=2
table(keep)
#14391
counts = as.data.frame(X[keep,])

gene.inf.1 = gene.inf[keep,]

table(gene.inf.1$Gene_type)


##### ready for exploratoray analysis#######
#######  ready for exploratory analysis

nGenes <- dim(counts)[1]
dim(counts)
summary(counts)


rownames(counts)=gene.inf.1$Ensembl_Gene_ID

geneLength =gene.inf.1$Length

geneName =gene.inf.1$Gene_name

n.p=ncol(counts)/2

# distribution of counts per sample
densityPlot(counts=counts, group=group, col=col,outfile=F)
par(new=T)
plot(density(log2(counts[,1]+1)),col="black",lwd=4)
majSequences(counts, n=10, group,geneName=geneName, col=col, outfile=F)
topgenelist = topgene(counts=counts, group=group, n=10,geneName=geneName,col=col,outfile=F)
write.csv(topgenelist,file=paste0(outdir,"top10gene.csv"))


##################################normalization########################################

i=1
func.0 = function(m){
  df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
  names(df) =paste0(m,".",names(df))
  return(df)
}


# norm method  
m="none"
m="tmm"
m="qq"
fn=paste0(m,"_bar_density")


tmp.0= func.0(m)



png(filename=paste0("figures/",fn,".png"),width=1200,height=500,units = "px")
par(mfrow=c(1,2))
boxplot(log2(tmp.0 + 1), col = rep(col[1:length(norm)], each=2*n.p), las = 2,
        main = "Counts distribution", ylab = expression(log[2] ~ (norm ~ count + 1)),outline=F)

densityPlot(counts=tmp.0, group=group, col=col,outfile=F)
dev.off()


############################## DE analysis#################
m="tmm"
# y.test= runedgeR1(m,counts,group, geneLength)

scale = c("none","tmm", "rle","mrn","tc","ps","uq","tbt")
setCounts =c("qq","ruvg")
if (m%in%(scale)){
  f=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[1]]
  y <- DGEList(counts=counts, group=group, norm.factors = f)
}
if (m%in%(setCounts)){
  df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
  y <- DGEList(counts=df, group=group, norm.factors = rep(1, ncol(df)))
}



y <- estimateGLMCommonDisp(y, verbose = T)
y <- estimateGLMTagwiseDisp(y)

# y <- estimateCommonDisp(y, verbose = T)
# y <- estimateTagwiseDisp(y)
y.test = exactTest(y)




################ glm method#################
# design <- model.matrix(~group)
# y <- estimateGLMCommonDisp(y,design, verbose = T)
# y <- estimateGLMTagwiseDisp(y, design)
# fit <- glmFit(y, design)
# lrt <- glmLRT(fit)
# topTags(lrt)
# 
# y.test=lrt
### plot normalized result
# plotMDS(y, col=col[as.integer(group)])
# abline(h=0,v=0,lty=2,col="lightgray")
# clusterPlot(counts.trans=(tmp.0), group=group, outfile=F)


y.out=topTags(y.test,n=nrow(counts),adjust.method="BH",sort.by="none")$table
sum(y.out$FDR<0.05)
y.out$id =rownames(y.out)
# colnames(y.out)[c(4,5)] =c("pval","padj")

colnames(y.out)[c(3,4)] =c("pval","padj")
res=y.out
par(mfrow=c(1,1))
hist(-log10(res$pval),xlab="Raw -log10(p-value)", main="Histogram of raw p-value")

hist(-log10(res$padj),xlab=" -log10(adj p-value)", main="Histogram of adj p-value")

sum(gene.inf$Chromosome=="chrM")
sum(gene.inf.1$Chromosome=="chrM")
gene.inf.2 =gene.inf.1[,c(2,23:30)]

res1=cbind(res, gene.inf.1)
res2=res1[order(res1$padj,decreasing = F),]
crit=0.05
crit =0.0001
out=res2[res2$padj < crit,]

write.csv(out1, file=paste0(outdir,"most_sig_genes_TMM_N.csv"))
# 9890
cm =out[out$Chromosome=="chrM",]
out1=out[out$Chromosome!="chrM" & out$Gene_type %in% c("protein_coding","pseudogene"),]
out$Chromosome!="chrM" & out$Gene_type %in% c("protein_coding","pseudogene")


############### visualize the result###########

crit.p =0.05
crit.lfc =0
is.de = decideTestsDGE(y.test, adjust.method="BH", p.value=crit.p, lfc=crit.lfc)

table(is.de)

incd=gene.inf.1$Chromosome!="chrM" & gene.inf.1$Gene_type %in% c("protein_coding","pseudogene")
res2=cbind(res1,is.de, incd)
de.out=res2[is.de & incd,]
table(de.out$is.de)


chrM=res2[res2$Chromosome=="chrM",]
summary(is.de)
plotMD(y.test, status=is.de, values=c(1,-1), col=c("red","blue"),ylab="Log2 fold change", main="N vs. CB")



d =res2

with(d, plot(logFC, -log10(padj), pch=20, xlab="Log2 fold change", ylab=expression(-log[10] ~ adjusted ~ p-value)))

with(subset(d, (logFC)>0 & padj<.05), points(logFC, -log10(padj), pch=20, col="red"))
with(subset(d, padj<.05 & logFC < 0), points(logFC, -log10(padj), pch=20, col="blue"))
# install.packages("calibrate")

############################## qvalue cut-off#################
summary(d$padj)


out=NULL
norm=c("tmm", "tbt","qq")

i=1

m="tmm"

  
 qva.f= function(m){
    y.test= runedgeR1(m,counts,group, geneLength)
    y.out=topTags(y.test,n=nrow(counts),adjust.method="BH",sort.by="none")$table
    sum(y.out$FDR<0.05)
    y.out$id =rownames(y.out)
    colnames(y.out)[c(4,5)] =c("pval","padj")
    # colnames(y.out)[c(3,4)] =c("pval","padj")
    
    res=y.out
    
    qv =c(0.05,0.01,0.005,0.001,0.0001)
    
    
    decideDE=function(y.test,crit.p,crit.lfc){
      
      is.de = decideTestsDGE(y.test, adjust.method="BH", p.value=crit.p, lfc=crit.lfc)
      table(is.de)
      return(is.de)
      
    }
    
    f.o=NULL
    for (j in (1:length(qv))){
      crit.p= qv[j]
      print(crit.p)
      is.de=decideDE(y.test,crit.p = crit.p,crit.lfc =1)
      incd=gene.inf.1$Chromosome!="chrM" & gene.inf.1$Gene_type %in% c("protein_coding","pseudogene")
      
      res2=cbind(res,is.de, incd)
      de.out=res2[is.de & incd,]
      f= c(m,crit.p,sum(de.out$is.de==-1), sum(de.out$is.de==1))
      f.o=rbind(f.o, f)
    }
    
    print(f.o)
    
  }
 print(m)
 
 y.o=NULL
 for (i in 1:length(norm)){
   m=norm[i]
   print(m)
   y.0=qva.f(m)
   y.o=cbind(y.o, y.0)
   
 }

write.csv(y.o, file=paste0(outdir,"qvalue_lfc1.csv"))

##############

# original counts
ct=rpkm(counts,gene.length = geneLength)
incd=gene.inf.1$Chromosome!="chrM" & gene.inf.1$Gene_type %in% c("protein_coding","pseudogene")

ct1=cbind(gene.inf.2,res, ct)[incd,]

ave.cb= rowMeans(ct1[,grep("CB",colnames(ct1))])
ave.n= rowMeans(ct1[,grep("N",colnames(ct1))])
summary(ave.cb)
summary(ave.n)
ord.cb=rank(-ave.cb)
ord.n = rank(-ave.n)
adult.fetus = ord.n - ord.cb 

gene.rpkm.ord= as.data.frame(cbind(ct1,ave.n, ave.cb,ord.n,ord.cb,adult.fetus))


# gene.tot[which(gene.tot$ID=="ENSG00000151789.5"),]

gene.tot =gene.rpkm.ord

topn =100

t1=(gene.tot$ord.n %in% 1:topn)
t2=(gene.tot$ord.cb %in% 1:topn)

## save top 100 genes 
n.top=gene.tot[order(gene.tot$ord.n, decreasing = F),][1:topn,]
head(n.top)
# write.csv(n.top,file =paste(outdir,"n_top_",topn,"_",Sys.Date(),".csv",sep=""))

cb.top=gene.tot[order(gene.tot$ord.cb, decreasing = F),][1:topn,]
# write.csv(cb.top,file =paste(outdir,"cb_top_",topn,"_",Sys.Date(),".csv",sep=""))


t3 = cbind(t1,t2)
t3.1 =t3[rowSums(t3)>0,]
library(limma)

t4 = vennCounts(t3.1)
t4
fname =paste("Venn",topn,sep="")
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 600, height = 600, units = "px", pointsize = 20) 
vennDiagram(t4, include = "both", 
            names = c("Abundant in N", "Abundant in CB"), 
            cex = 1,  circle.col=c("green","blue"),counts.col = "black")
# title(paste("Top ", topn," genes" ,sep=""))
title("")
dev.off()
sum(gene.both.high$logFC>2)


gene.both.high = gene.tot[rowSums(t3)>1,]

both.high.sorted=gene.both.high[order(gene.both.high$padj,decreasing = F),]
sum(gene.both.high$padj<0.05) 

write.csv(both.high.sorted,file=paste0(outdir,"both_high_sorted.csv"))


