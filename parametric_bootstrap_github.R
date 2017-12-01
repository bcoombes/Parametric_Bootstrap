##################################################################################
## Author: Brandon Coombes
## Date: 11/30/2017
## Purpose: Code to be uploaded to GitHub for performing parametric bootstrap
##
## Requirements: GMMAT R package
##                  -GMMAT can be downloaded here: https://content.sph.harvard.edu/xlin/software.html#gmmat
##                      This package is being updated to run in GENESIS bioconductor package 
##                      but implementation in GENESIS may be different so best to download from above
##               parallel R package (optional if you want to parallelize bootstrap [highly recommended])  
##                
###############################################################################################

library(GMMAT)
library(parallel)

keep_PCA <- function(G,R2){ 
  #----------------------------------------
  # Creates PC matrix of genotype matrix 
  #
  # INPUTS:
  # G = SNP matrix with each column representing each SNP coded 0-2 (could be dosage value)
  # R2 = a numeric value giving percent of the variation to retain with the PCA output
  #
  # OUTPUT: a numeric matrix of PCA's that explain 'R2' percent of variance
  #----------------------------------------
  fit <- princomp(G,corr=T)
  var.PC <- as.numeric(fit$sdev)^2
  total.var <- sum(var.PC)
  greater.index <- which(cumsum(var.PC)>R2*total.var)
  keep.count <- min(greater.index)
  total.PCA <- fit$scores
  keep.PCAs <- total.PCA[,1:keep.count] # only keep PCA's that explain specified percent of variance for a given gene
  return(keep.PCAs)
}

fisher.p <- function(x){
  #---------------------------------------
  # Fisher's method of combining p-values
  #  
  # INPUTS:
  # x = set of p-values from genes
  #
  # OUTPUT: Overall combination p-value
  #---------------------------------------
  Ts <- -2*sum(log(x),na.rm=TRUE)
  pchisq(Ts,df=2*sum(!is.na(x)),lower.tail=FALSE)  
}

gamma.p <- function(x,STT){
  #---------------------------------------
  # Gamma method of combining p-values
  #  
  # INPUTS:
  # x = set of p-values from genes
  # STT = soft truncation threshold (Equivalent to Fisher's with STT = 1/exp(1))
  #
  # OUTPUT: Overall combination p-value
  #---------------------------------------
  sh <- STTtoShapeParameter(STT)
  Ts <- sum(qgamma(1-x,shape=sh),na.rm=TRUE) 
  pgamma(Ts,shape=sh*sum(!is.na(x)),lower.tail=FALSE)
}

STTtoShapeParameter <-function(STT) {
  # ----------------------------------------
  # Converts provided STT to a shape parameter in Gamma distribution
  #
  # INPUTS:
  # STT = soft truncation threshold 
  #
  # OUTPUT: shape parameter
  #----------------------------------------
  if(STT> 0.4) stop("only works for STT <= 0.4")
  f<-function(w,STT) abs(w-qgamma(1-STT,shape=w))
  out<-optimize(f,lower=0,upper=2,STT=STT)$minimum
  return(out)
}

lrtest.GE <- function(Y,G,E,X=NULL){
  #----------------------------------------------
  # Compute score test p-value for a GxE 
  #
  # Y = phenotype (binary or quantitative)
  # G = matrix of SNPs or PCs of SNPs
  # E = Matrix of environments (binary or quantitative)
  # X = Matrix of covariates 
  #---------------------------------------------
  if (length(unique(Y))==2){fam<-"binomial"}else{fam<-"gaussian"}
  if (is.null(X)){
    out <- anova(glm(Y~G*E,family=fam),test="Rao")$"Pr(>Chi)"[4]
  }else{out <- anova(glm(Y~G*E+X,family=fam),test="Rao")$"Pr(>Chi)"[5]}
  return(out)
}

Geno.matrix.to.list <- function(G,snp.map,R2=NULL){
  #----------------------------------------------
  # Take a user genotype matrix and set up into a list of genotypes for each gene
  # 
  # G = genotype matrix or dataframe containing SNPs on each column with rsIDs
  # snp.map = dataframe with columns "snp" and "gene" the snp belongs to
  # R2 = option to convert each gene genotypes to PCs with some R2 specified
  #----------------------------------------------
  genes <- unique(snp.map$gene)
  snpdat <- PCdat <- list()
  for (i in genes){
    #print(i)
    snpdat[[i]] <- as.matrix(G[,snp.map[snp.map$gene==i,"snp"]])
    if (!is.null(R2)){
      PCdat[[i]] <- keep_PCA(snpdat[[i]],R2)  
    }
  }
  return(list(Glist=snpdat,PClist=PCdat))
}


parboot.GE <- function(Y,Glist,PClist=NULL,E,X=NULL,B=1000,numcores=1){
  #----------------------------------------------
  # Perform parametric bootstrap to test for GE interaction in pathway
  # 
  # Y = phenotype (binary or quantitative)
  # Glist = list containing genotype matrices for each gene
  # PClist = list containing PC matrices for each gene (optional)
  # E = environment of interest
  # X = matrix of covariates 
  # B = number of bootstraps to be performed
  # numcores = number of cores to split the jobs into (highly recommended > 1)
  #----------------------------------------------
  if (length(unique(Y))==2){fam<-"binomial"}else{fam<-"gaussian"}
  if (is.null(PClist)){G <- Glist}else{G <- PClist}
  ## compute gene-level p-values (can insert your own GE interaction gene-level test if needed instead of lrtest.GE)
  if (numcores==1){
    gene.ps <- sapply(1:length(G),function(x) lrtest.GE(Y,G[[x]],E,X))
  }else{gene.ps <- unlist(mclapply(1:length(G),function(x) lrtest.GE(Y,G[[x]],E,X),mc.cores=numcores))}
  STTs <- c(0.01,0.05,0.1,0.15,0.2,1/exp(1))  ## I search multiple STTs including Fisher's
  path.p <- sapply(STTs,function(x) gamma.p(gene.ps,STT=x)) 
  path.p <- c(path.p,min(path.p)) ## also take minimum p-value over search
  names(path.p) <- paste("Gamma",c(STTs[1:(length(STTs)-1)],"Fisher","min"),sep="_")
  
  ###############################
  ## perform parametric bootstrap
  ###############################
  Gall <- do.call(cbind,Glist) #make large genotype matrix for whole pathway
  gene.GRM <-list(as.matrix(Gall%*%t(Gall))) #specify covariance of random effect
  DAT <- data.frame(Y,X,E) #needed for GMMAT 
  
  if (fam=="binomial"){
    null.model<-try(glmmkin(Y~1+X+E,data=DAT,kins=gene.GRM,family=binomial(link="logit")))
    pis <- null.model$fitted.values
    n <- length(Y)
    if(numcores==1){
      null.ps <- t(sapply(1:B,function(y){
        null.gene.ps <- sapply(1:length(G), function(x){Yb <- rbinom(n,size=1,p=pis); lrtest.GE(Yb,G[[x]],E,X)})
        null.path.p <- sapply(STTs,function(x) gamma.p(null.gene.ps,STT=x)) 
        null.path.p <- c(null.path.p,min(null.path.p)) ## also take minimum p-value over search
      }))
      dimnames(null.ps)[[2]] <- paste("Gamma",c(STTs[1:(length(STTs)-1)],"Fisher","min"),sep="_")
    }else{
      null.ps <- do.call(rbind,mclapply(1:B,function(y){
        null.gene.ps <- sapply(1:length(G), function(x){Yb <- rbinom(n,size=1,p=pis); lrtest.GE(Yb,G[[x]],E,X)})
        null.path.p <- sapply(STTs,function(x) gamma.p(null.gene.ps,STT=x)) 
        null.path.p <- c(null.path.p,min(null.path.p)) ## also take minimum p-value over search
      },mc.cores=numcores))
      dimnames(null.ps)[[2]] <- paste("Gamma",c(STTs[1:(length(STTs)-1)],"Fisher","min"),sep="_")
    }
  }
  if (fam=="gaussian"){
    null.model<-try(glmmkin(Y~1+X+E,data=DAT,kins=gene.GRM,family=gaussian(link="identity")))
    mu <- null.model$fitted.values
    n <- length(Y)
    sig <- sqrt(sum((Y-mu)^2)/n)
    if(numcores==1){
      null.ps <- t(sapply(1:B,function(y){
        null.gene.ps <- sapply(1:length(G), function(x){Yb <- rnorm(n,mu,sig); lrtest.GE(Yb,G[[x]],E,X)})
        null.path.p <- sapply(STTs,function(x) gamma.p(null.gene.ps,STT=x)) 
        null.path.p <- c(null.path.p,min(null.path.p)) ## also take minimum p-value over search
      }))
      dimnames(null.ps)[[2]] <- paste("Gamma",c(STTs[1:(length(STTs)-1)],"Fisher","min"),sep="_")
    }else{
      null.ps <- do.call(rbind,mclapply(1:B,function(y){
        null.gene.ps <- sapply(1:length(G), function(x){Yb <- rnorm(n,mu,sig); lrtest.GE(Yb,G[[x]],E,X)})
        null.path.p <- sapply(STTs,function(x) gamma.p(null.gene.ps,STT=x)) 
        null.path.p <- c(null.path.p,min(null.path.p)) ## also take minimum p-value over search
      },mc.cores=numcores))
      dimnames(null.ps)[[2]] <- paste("Gamma",c(STTs[1:(length(STTs)-1)],"Fisher","min"),sep="_")
    }
  }
  corrected.path.p <- sapply(1:length(path.p),function(x) mean(path.p[x]>null.ps[,x]))  ## compare observed to null 
  return(list(asym=path.p,parboot=corrected.path.p))
}

