
## compare the normalized counts

## distribution of counts
    
    # 13 in total
norm=c("none","tmm","tbt", "rle","mrn","ps","uq","qq","ruvg")
norm.1=c("rpkm","tpm")
   
    setwd(outdir)
#######fig1a boxplot of normalized counts    
    i=1
    func.0 = function(m){
      df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
      names(df) =paste0(m,".",names(df))
      return(df)
    }
    
     tmp.0= func.0("tpm")
# norm method  
    fn ="Fig1a"
    norm =norm
    tmp.0= do.call(cbind, lapply(norm, func.0))
    
    png(filename=paste0("figures/",fn,".png"),width=900,height=500,units = "px")
    boxplot(log2(tmp.0 + 1), col = rep(col[1:length(norm)], each=2*n.p), las = 2,
            main = "Counts distribution", ylab = expression(log[2] ~ (norm ~ count + 1)),outline=FALSE)
    dev.off()
    
#######fig1b boxplot of normalized counts    
    
    fn ="Fig1b"
    png(filename=paste0("figures/",fn,".png"),width=900,height=500,units = "px")
    par(mfrow =c(1,2))
    tmp.0= func.0("rpkm")
    boxplot(log2(tmp.0 + 1), col = rep(col[5], each=2*n.p), las = 2,
              main = "Counts distribution", ylab = expression(log[2] ~ (rpkm + 1)),outline=FALSE)
    tmp.0= func.0("tpm")
    boxplot(log2(tmp.0 + 1), col = rep(col[11], each=2*n.p), las = 2,
            main = "Counts distribution", ylab = expression(log[2] ~ (tpm ~ count + 1)),outline=FALSE)
    
    dev.off()

######################################## within group variataion#########################
  ############# 
    
    i=1
    m="tbt"
    func.1 = function(m){
      df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
      sigma = apply(df,2, sd)
      mu = apply(df,2, mean)
      (cv.0= sigma/mu)
      
      return(cv.0)
    }

   (norm =norm)
  
    # func.1("none")
    tmp.0= do.call(rbind, lapply(norm, func.1))
    rownames(tmp.0) = norm
    
    
    fn ="Fig2aCV"
    png(filename=paste0("figures/",fn,".png"),width=900,height=500,units = "px")
    par(mfrow =c(1,2))
    d1 = t(tmp.0[,c(grep("pre", colnames(tmp.0)))])
    boxplot(d1, col = (col[1:ncol(d1)]), las = 2,
            main = "Intra-group coef of variation (Pre)", ylab = "Coefficient of variance",outline=FALSE)
    d1 = t(tmp.0[,c(grep("post", colnames(tmp.0)))])
    boxplot(d1, col = (col[1:ncol(d1)]), las = 2,
            main = "Intra-group coef of variation (Post)", ylab = "Coefficient of variance",outline=FALSE)
      dev.off()
    
    

    