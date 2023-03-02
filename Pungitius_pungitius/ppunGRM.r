##################################
### Genetic relatedness matrix ###
##################################

library(snpReady)
library(SNPRelate)
library(tidyverse)
library(stringr)
library(MCMCglmm)
library(QGglmm)
library(purrr)

### I - Data preparation

# 1a - Import all raw data ###

geno = read.table("helsinkiHeaderedSorted.012",na.strings = "-1",stringsAsFactor=F)#load data
geno.ids=read.table("helsinkiHeaderedSorted.012.indv",stringsAsFactor=F)
geno$V1=geno.ids$V1
geno$V1=gsub(x=geno$V1,pattern="42-f-NA",replacement = "42-f")
marker.names=read.table("helsinkiHeaderedSorted.012.pos",stringsAsFactors = F,header=F,sep="\t")
##
n.markers.used=ncol(geno)-1

##
#PHENO DATA
data.pheno=read.table("Ppun_Helsinki_Crossovers_toAnimalModel.csv", header=T,stringsAsFactors = F, sep="\t")#From collectCrossovers.r
##
### II - Prepare genotype data 
pardata=geno

##
X1 <- pardata[,2:ncol(pardata)]
colnames(X1) <- paste0(marker.names$V1,":",marker.names$V2)#markers 
X3 = as.matrix(X1)
X4 = matrix(data=X3, nrow=nrow(X1), ncol=ncol(X1))
X4 = apply(X3, 2, as.numeric)
colnames(X4)= paste0(marker.names$V1,":",marker.names$V2)#markers
rownames(X4)=pardata$V1

## 2a - Make GRM

PpunHel <- raw.data(X4, frame="wide", base=FALSE, sweep.sample= 0.95,call.rate=0.50, maf=0.01, imput=T, imput.type = "mean")
PpunHelG = G.matrix(PpunHel$M.clean, method="VanRaden", format="wide", plot = F)
G_PpunHEL1 = PpunHelG$Ga

#Keep parents only
G_PpunHEL1_parents=G_PpunHEL1[which(rownames(G_PpunHEL1) %in% data.pheno$animal),which(colnames(G_PpunHEL1) %in% data.pheno$animal)]

## 2b - Format GRM to sparse for MCMCglmm

N=nrow(G_PpunHEL1_parents) #number of individuals in the dataset
i <- rep(1:N,rep(N,N))
j <- rep(1:N,N)

s <-spMatrix(N,N,i,j,as.vector(G_PpunHEL1_parents))
snpreadyGRM<-solve(s)
class(snpreadyGRM) <- "dgCMatrix"
rownames(snpreadyGRM) <- snpreadyGRM@Dimnames[[1]] <- rownames(G_PpunHEL1_parents)
rownames(snpreadyGRM) <- snpreadyGRM@Dimnames[[2]] <- rownames(G_PpunHEL1_parents)

#save(snpreadyGRM,file="PpunGRM.rds") 
