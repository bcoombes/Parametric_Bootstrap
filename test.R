####### Example implementation of the code
### Simulate simple Y, G, E, and X to demonstrate the code above
set.seed(1930)
snp.map <- data.frame(snp=paste("rs",1:1000,sep=""),gene=paste("Gene",rep(1:20,each=50),sep=""),stringsAsFactors = F)
G <- sapply(1:1000,function(x) rbinom(500,2,0.3)) #create 1000 independent SNPs for 500 subjects
dimnames(G)[[2]] <- snp.map[,1] #names SNPs according to above fake rsIDs
binary <- rbinom(500,1,0.5) #binary variable used for either Y or E
X <- cbind(rbinom(500,1,0.5),runif(500,20,80)) #covariates representing sex and age
mean.model <- G[,1:5]%*%rep(0.5,5)+0.05*binary+X%*%c(0.2,0.01) + (G[,1:5]*binary)%*%rep(0.5,5) #sample mean model with interaction with onl$
normvar <- rnorm(500,mean=mean.model,sd=1) #normal variable dependent on the other

summary(glm(normvar~G[,1:5]*binary+X,family="gaussian")) ##check the correct disease model
summary(glm(binary~G[,1:5]*normvar+X,family="binomial")) ##check the reversed disease model

## map SNPs to genes and set up for analysis
gene.mats <- Geno.matrix.to.list(G,snp.map,0.8)
Glist <- gene.mats$Glist
PClist <- gene.mats$PClist

#binary version with only 100 bootstraps
parboot.GE(Y=binary,Glist = Glist,PClist = PClist,E = normvar,X=X,B = 100,numcores = 10) #parallelized
parboot.GE(Y=binary,Glist = Glist,PClist = PClist,E = normvar,X=X,B = 10,numcores = 1) #no parallelization

#quantitative version with only 100 bootstraps
parboot.GE(Y=normvar,Glist = Glist,PClist = PClist,E = binary,X=X,B = 100,numcores = 10) #parallelized
parboot.GE(Y=normvar,Glist = Glist,PClist = PClist,E = binary,X=X,B = 10,numcores = 1) #no parallelization

