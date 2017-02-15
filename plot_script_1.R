## explore the difference
# gene.rpkm.raw = rpkm(edgeR.dgelist,normalized.lib.sizes=F)

gene.rpkm.0 = t(scale(t(gene.rpkm.raw)))
gene.rpkm.0 = t(scale(t(gene.rpkm)))
colnames(gene.rpkm.0) = sampleDict 
sampleDists = dist( t( gene.rpkm.0 ) )

sampleDists = dist( t( gene.rpkm.0 ) )

sampleDists=as.dist(sqrt(1-(cor(gene.rpkm.0, method="spearman"))^2))



# library('corrplot') #package corrplot
# install.packages("corrplot")
# corrplot(t, method = "circle") #plot matrix
# sampleDists
# install.packages("pheatmap")
library("pheatmap")
library("RColorBrewer")
sampleDists=as.dist(sqrt(1-(cor(X))^2))
sampleDistMatrix <- as.matrix( sampleDists )

plot(hclust(sampleDists), 
     main="Dissimilarity = 1 - Correlation", xlab="")
cor.1=cor(X)
plot(hclust(dist(abs(cor.1))))
colSums(X)

fname ="cluster_samples"
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 950, height = 500, units = "px", pointsize = 12) 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


x= sampleDistMatrix
#  d<-scale(gene.rpkm)
#  colnames(d) = sampleDict 
# x <- as.matrix( dist(t(d)))
fname ="MDS_plot"
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 550, height = 550, units = "px", pointsize = 16) 
plot(cmdscale(x), xlab="Coordinate 1", ylab="Coordinate 2", type = "n") ; 
text(cmdscale(x), labels=colnames(x)) 
dev.off()


tmp=read.delim(file=paste(outdir,"gsea_report_for_na_pos_1476992616496.xls",sep=""))

tmp=read.delim(file=paste(outdir,"gsea_report_for_na_neg_1476992616496.xls",sep=""))
topn=20
tmp.1= tmp[1:topn,c("NAME","NES")][topn:1,]


fname ="GSEA_adult"
fname ="GSEA_CB"
png(filename = paste(figdir,fname,Sys.Date(),".png",sep=""),
    width = 750, height = 550, units = "px", pointsize = 16) 
op <- par(mar = c(2,24,2,2) + 0.1)
barplot(tmp.1$NES, horiz=T,las=2,space=0)
axis(2, at =0.5:19.5 ,labels =tmp.1$NAME,las=2, cex.axis=0.7)

dev.off()
par(op) ## reset