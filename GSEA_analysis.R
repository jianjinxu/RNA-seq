
# created on 09/26/2016 to do GSEA
rm(list =ls())
library(edgeR)
workdir = "C:/Users/jackie/Google Drive/projects/RNA-seq/"
datadir = paste(workdir,"data/",sep="")
figdir = paste(workdir,"figure/",sep="")
outdir = paste(workdir,"output/",sep="")
refdir = paste(workdir,"to_deliver/",sep="")
source("https://bioconductor.org/biocLite.R")
# biocLite("GSEABase")
# library("GSEABase")
biocLite(c("biomaRt", "topGO", "org.Mm.eg.db"))
library(biomaRt)
library(org.Mm.eg.db)
library(topGO)  

listEnsembl()
listEnsembl("GRCh=37")
ensembl = useEnsembl(biomart="ensembl")
head(listDatasets(ensembl))
grch37 = useEnsembl(biomart="ensembl",GRCh=37)
listDatasets(grch37)[31:35,]
listDatasets(grch37)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
listAttributes(ensembl) 

genes.with.id=getBM(attributes=c("ensembl_gene_id", "external_gene_id","entrezgene"),values=gene_names, mart= ensembl)



###############

biocLite("DESeq2")
library(DESeq2)

 coldat=DataFrame(grp=factor(class))
 dds <- DESeqDataSetFromMatrix(test, colData=coldat, design = ~ grp)
 dds <- DESeq(dds)
deseq2.res <- results(dds)
#direction of fc, depends on levels(coldat$grp), the first level
#taken as reference (or control) and the second one as experiment.
deseq2.fc=deseq2.res$log2FoldChange
 names(deseq2.fc)=rownames(deseq2.res)
 exp.fc=deseq2.fc
 out.suffix="deseq2"
 biocLite("gage")
 library(gage)
 data(kegg.gs)
 fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &  !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
 sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
 path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
 
 biocLite("pathview")
 library(pathview)
 #view first 3 pathways as demo
   pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview( gene.data = exp.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix))
   