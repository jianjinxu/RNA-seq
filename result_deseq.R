# read in result analyzed by the New York Genome center

workdir = "/Users/jianjinxu/Google Drive/projects/RNA-seq/"
datadir = paste(workdir,"data/",sep="")

refdir = paste(workdir,"to_deliver/",sep="")
deseq.gene="DESeq2_results_CB_vs_N.txt"
res.0=read.delim(file=paste(refdir,deseq.gene,sep=""),sep="\t",header = T)
# 15342 DE
out.deseq=res.0[which(res.0$padj<0.05),]
sum(!is.na(res.0$padj))
# 28658 genes enter the study
table(p1)
# 15342 DE, 13316 EE
sum(res.0)

######### order ########
####### top 100 DE genes######
