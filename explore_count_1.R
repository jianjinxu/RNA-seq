
############ exploratory result##############
setwd(outdir)
outfile =T
if (outfile) png(filename="figures/barplotRaw.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
boxplot(log2(counts+ 1), col = rep(col[1:ncol(counts)], each=2*n.p), las = 2,
        main = "Counts distribution", ylab = expression(log[2] ~ (raw ~ count + 1)),outline=F)
if (outfile) dev.off()
par(mfrow=c(1,1))
# all the data visualization process
topgenelist= desVisu(counts, group,col=col,outfile=T)
# cluster and PCA,MDS for raw data
par(mfrow=c(1,1))
PCAPlot(counts.trans=counts, group=group, col=col, outfile = F)


MDSPlot(dge=count1, group=group, col=col, outfile=F)

MDSPlot(dge=data.frame(ct.cpm), group=group, col=col, outfile=F)
clusterPlot(counts.trans=(counts), group=group, outfile=F)
clusterPlot(counts.trans=rpkm(counts,gene.length = geneLength), group=group,fn="cluster_rpkm_remove", outfile=F)


clusterPlot(counts.trans=cpm(counts), group=group,fn="cluster_rpkm_remove", outfile=F)


ct.cpm=cpm(count1)
ct=rpkm(counts,gene.length = geneLength)
