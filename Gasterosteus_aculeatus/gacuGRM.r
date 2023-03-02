#######################################
### RECOMBINATION RATE HERITABILITY ###
#######################################

library(snpReady)
library(MCMCglmm)

### I - Data preparation


geno = read.table('../../relatedness/3spHeadered.vcf.012',na.strings = "-1",stringsAsFactor=F)#load data
geno.ids=read.table('../../relatedness/3spHeadered.vcf.012.indv',stringsAsFactor=F)
geno$V1=geno.ids$V1
marker.names=read.table("../../relatedness/3spHeadered.vcf.012.pos",stringsAsFactors = F,header=F,sep="\t")
##
n.markers.used=ncol(geno)-1

pheno = read.table('../../data/gacuCrossovers.csv',header = TRUE,sep="\t",stringsAsFactor=F)
data.pheno=pheno

### II - Prepare genotype data 
pardata=geno
##
X1 <- pardata[,2:ncol(pardata)]
colnames(X1) <- paste0(marker.names$V1,":",marker.names$V2)#markers 
X3 = as.matrix(X1)
X4 = matrix(data=X3, nrow=nrow(X1), ncol=ncol(X1))
X4 = apply(X3, 2, as.numeric)
colnames(X4)=paste0(marker.names$V1,":",marker.names$V2)#markers
rownames(X4)=pardata$V1
## 2a - Make GRM

GacuHel <- raw.data(X4, frame="wide", base=FALSE, sweep.sample= 0.95,call.rate=0.50, maf=0.01, imput=TRUE, imput.type = "mean")
GacuHelG = G.matrix(GacuHel$M.clean, method="VanRaden", format="wide", plot = F)
G_GACUHEL1 = GacuHelG$Ga
#Keep parents only
G_GACUHEL_parents=G_GACUHEL1[which(rownames(G_GACUHEL1) %in% c(data.pheno$MALE,data.pheno$FEMALE)),which(colnames(G_GACUHEL1) %in% c(data.pheno$MALE,data.pheno$FEMALE))]

## 2b - Format GRM to sparse for MCMCglmm
N=nrow(G_GACUHEL_parents) #number of individuals in the dataset
i <- rep(1:N,rep(N,N))
j <- rep(1:N,N)


s <-spMatrix(N,N,i,j,as.vector(G_GACUHEL_parents))
snpreadyGRM<-solve(s)
class(snpreadyGRM) <- "dgCMatrix"
rownames(snpreadyGRM) <- snpreadyGRM@Dimnames[[1]] <- rownames(G_GACUHEL_parents)#with(pardata,V1)
rownames(snpreadyGRM) <- snpreadyGRM@Dimnames[[2]] <- rownames(G_GACUHEL_parents) #with(pardata,V1)

#save(snpreadyGRM,file="gacuGRM.rds") 


