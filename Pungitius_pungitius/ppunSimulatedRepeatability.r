###################################################
### RECOMBINATION RATE REPEATABILITY SIMULATION ###
###################################################
library(MCMCglmm)
###
pascal.likelihoods=function(n){
  likelihoods=c()
  for(i in 0:n){
    likelihoods=c(likelihoods,choose(n,i)*0.5^n)
  }
  return(likelihoods)
}

generate.gamete.yf=function(probs){
  n.cos=sample(x=0:9,size=1,prob=probs)
  obs.co=sample(x=0:n.cos,size=1,prob=pascal.likelihoods(n.cos))
  return(obs.co)
}

### I - Data preparation
pheno.long = read.table('Ppun_Helsinki_Crossovers_toAnimalModel.csv',header = TRUE,sep="\t",stringsAsFactor=F)

fem=subset(pheno.long,parentsex=="Female")
fem=droplevels(fem)

male=subset(pheno.long,parentsex=="Male")
male=droplevels(male)

#Maximum-likelihood estimates
ml.paternal=read.table("",header=T,stringsAsFactors = F,sep="\t")#Kivikoski et al. Supplementary table 2
ml.maternal=read.table("",header=T,stringsAsFactors = F,sep="\t") #Kivikoski et al. Supplementary table 1


### III - Run MCMCglmmm
BURN = 10000
THIN = 50
NITT = 110000
prior.C<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))

#Generate data #Takes about 10 sec per 10 
female.pop.grand.mean=c()
female.population.variance=c()
female.variance.of.parental.means=c()
Nsimulations=1000
for(i in 1:Nsimulations){
  crossovers=matrix(nrow=21,ncol=934)
  for (r in 1:21){
    crossovers[r,]=replicate(expr=generate.gamete.yf(as.numeric(ml.maternal[r,paste0("MLEpRestricted",0:9)])),n=934)
  }
  fem[,paste0("simulatedOffspring",i)]=colSums(crossovers[-12,])
  fem$tmp=colSums(crossovers[-12,])
  var.parent.means=var(aggregate(tmp ~ animal,mean,data=fem)$tmp)
  female.pop.grand.mean=c(female.pop.grand.mean,mean(mean(fem$tmp)))
  female.variance.of.parental.means=c(female.variance.of.parental.means,var.parent.means)
  female.population.variance=c(female.population.variance,var(fem$tmp))
}

prior=prior.C
#FEMALE REPEATABILITY

Va.out=c()
Va.posterior.mean.out=c()
Va.posterior.median.out=c()
Va.posterior.mode.out=c()

Vpe.out=c()
Vpe.posterior.mean.out=c()
Vpe.posterior.median.out=c()
Vpe.posterior.mode.out=c()

Ve.out=c()
Ve.posterior.mean.out=c()
Ve.posterior.median.out=c()
Ve.posterior.mode.out=c()

R.out=c()
R.posterior.mean.out=c()
R.posterior.median.out=c()
R.posterior.mode.out=c()

Vp.out=c()
Vp.posterior.mean.out=c()
Vp.posterior.median.out=c()
Vp.posterior.mode.out=c()

for(i in 1:Nsimulations){
  fem$COtmp=fem[,paste0("simulatedOffspring",i)]
  ppun.heritability.female.ranomized <- MCMCglmm(COtmp ~ 1,
                                                 random = ~ ID,
                                                 rcov = ~ units,
                                                 data = fem,
                                                 prior = prior,
                                                 family = c("gaussian"),pl=T,pr=T,
                                                 nitt = NITT, thin = THIN, burnin = BURN,verbose=F)

  Vpe=ppun.heritability.female.ranomized[["VCV"]][,"ID"]
  Vpe.posterior.mean=mean(Vpe)
  Vpe.posterior.median=median(Vpe)
  Vpe.posterior.mode=posterior.mode(Vpe)

  Ve=ppun.heritability.female.ranomized[["VCV"]][,"units"]
  Ve.posterior.mean=mean(Ve)
  Ve.posterior.median=median(Ve)
  Ve.posterior.mode=posterior.mode(Ve)

  R=ppun.heritability.female.ranomized[["VCV"]][,"ID"]/rowSums(ppun.heritability.female.ranomized[["VCV"]])
  R.posterior.mean=mean(R)
  R.posterior.median=median(R)
  R.posterior.mode=posterior.mode(R)

  Vp=rowSums(ppun.heritability.female.ranomized[["VCV"]])
  Vp.posterior.mean=mean(Vp)
  Vp.posterior.median=median(Vp)
  Vp.posterior.mode=posterior.mode(Vp)

  Vpe.out=c(Vpe.out,Vpe)
  Vpe.posterior.mean.out=c(Vpe.posterior.mean.out,Vpe.posterior.mean)
  Vpe.posterior.median.out=c(Vpe.posterior.median.out,Vpe.posterior.median)
  Vpe.posterior.mode.out=c(Vpe.posterior.mode.out,Vpe.posterior.mode)

  Ve.out=c(Ve.out,Ve)
  Ve.posterior.mean.out=c(Ve.posterior.mean.out,Ve.posterior.mean)
  Ve.posterior.median.out=c(Ve.posterior.median.out,Ve.posterior.median)
  Ve.posterior.mode.out=c(Ve.posterior.mode.out,Ve.posterior.mode)

  R.out=c(R.out,R)
  R.posterior.mean.out=c(R.posterior.mean.out,R.posterior.mean)
  R.posterior.median.out=c(R.posterior.median.out,R.posterior.median)
  R.posterior.mode.out=c(R.posterior.mode.out,R.posterior.mode)

  Vp.out=c(Vp.out,Vp)
  Vp.posterior.mean.out=c(Vp.posterior.mean.out,Vp.posterior.mean)
  Vp.posterior.median.out=c(Vp.posterior.median.out,Vp.posterior.median)
  Vp.posterior.mode.out=c(Vp.posterior.mode.out,Vp.posterior.mode)
}
female.out=data.frame(Species="P.pungitius",Sex="Female",Iteration=1:Nsimulations,
                    Vpe.mean=Vpe.posterior.mean.out, Vpe.median=Vpe.posterior.median.out, Vpe.mode=Vpe.posterior.mode.out,
                    Ve.mean=Ve.posterior.mean.out, Ve.median=Ve.posterior.median.out, Ve.mode=Ve.posterior.mode.out,
                    R.mean=R.posterior.mean.out, R.median=R.posterior.median.out, R.mode=R.posterior.mode.out,
                    Vp.mean=Vp.posterior.mean.out,Vp.median=Vp.posterior.median.out,Vp.mode=Vp.posterior.mode.out,
                    GrandMean=female.pop.grand.mean, VarianceOfParentalMeans=female.variance.of.parental.means)

#write.table(file="PpunSimulatedRepeatabilityFemale.csv",female.out,col.names = T, row.names = F,sep="\t", quote = F)




##SIMULATE MALE DATA
male.pop.grand.mean=c()
male.variance.of.parental.means=c()
male.population.variance=c()
Nsimulations=1000
for(i in 1:Nsimulations){
  crossovers=matrix(nrow=21,ncol=934)
  for (r in 1:21){
    crossovers[r,]=replicate(expr=generate.gamete.yf(as.numeric(ml.paternal[r,paste0("MLEpRestricted",0:9)])),n=934)
  }
  male[,paste0("simulatedOffspring",i)]=colSums(crossovers[-12,])
  male$tmp=colSums(crossovers[-12,])
  var.parent.means=var(aggregate(tmp ~ animal,mean,data=male)$tmp)
  male.pop.grand.mean=c(male.pop.grand.mean,mean(mean(male$tmp)))
  male.variance.of.parental.means=c(male.variance.of.parental.means,var.parent.means)
  male.population.variance=c(male.population.variance,var(male$tmp))
  
}

###
#MALE REPEATABILITY

Va.out=c()
Va.posterior.mean.out=c()
Va.posterior.median.out=c()
Va.posterior.mode.out=c()

Vpe.out=c()
Vpe.posterior.mean.out=c()
Vpe.posterior.median.out=c()
Vpe.posterior.mode.out=c()

Ve.out=c()
Ve.posterior.mean.out=c()
Ve.posterior.median.out=c()
Ve.posterior.mode.out=c()

R.out=c()
R.posterior.mean.out=c()
R.posterior.median.out=c()
R.posterior.mode.out=c()

Vp.out=c()
Vp.posterior.mean.out=c()
Vp.posterior.median.out=c()
Vp.posterior.mode.out=c()

for(i in 1:Nsimulations){
  male$COtmp=male[,paste0("simulatedOffspring",i)]
  ppun.heritability.male.ranomized <- MCMCglmm(COtmp ~ 1,
                                                 random = ~ ID,
                                                 rcov = ~ units,
                                                 data = male,
                                                 prior = prior,
                                                 family = c("gaussian"),pl=T,pr=T,
                                                 nitt = NITT, thin = THIN, burnin = BURN,verbose=F)

  Vpe=ppun.heritability.male.ranomized[["VCV"]][,"ID"]
  Vpe.posterior.mean=mean(Vpe)
  Vpe.posterior.median=median(Vpe)
  Vpe.posterior.mode=posterior.mode(Vpe)
  
  Ve=ppun.heritability.male.ranomized[["VCV"]][,"units"]
  Ve.posterior.mean=mean(Ve)
  Ve.posterior.median=median(Ve)
  Ve.posterior.mode=posterior.mode(Ve)
  
  R=(ppun.heritability.male.ranomized[["VCV"]][,"ID"])/rowSums(ppun.heritability.male.ranomized[["VCV"]])
  R.posterior.mean=mean(R)
  R.posterior.median=median(R)
  R.posterior.mode=posterior.mode(R)
  
  Vp=rowSums(ppun.heritability.male.ranomized[["VCV"]])
  Vp.posterior.mean=mean(Vp)
  Vp.posterior.median=median(Vp)
  Vp.posterior.mode=posterior.mode(Vp)
  
  Vpe.out=c(Vpe.out,Vpe)
  Vpe.posterior.mean.out=c(Vpe.posterior.mean.out,Vpe.posterior.mean)
  Vpe.posterior.median.out=c(Vpe.posterior.median.out,Vpe.posterior.median)
  Vpe.posterior.mode.out=c(Vpe.posterior.mode.out,Vpe.posterior.mode)

  Ve.out=c(Ve.out,Ve)
  Ve.posterior.mean.out=c(Ve.posterior.mean.out,Ve.posterior.mean)
  Ve.posterior.median.out=c(Ve.posterior.median.out,Ve.posterior.median)
  Ve.posterior.mode.out=c(Ve.posterior.mode.out,Ve.posterior.mode)

  R.out=c(R.out,R)
  R.posterior.mean.out=c(R.posterior.mean.out,R.posterior.mean)
  R.posterior.median.out=c(R.posterior.median.out,R.posterior.median)
  R.posterior.mode.out=c(R.posterior.mode.out,R.posterior.mode)
 
  Vp.out=c(Vp.out,Vp)
  Vp.posterior.mean.out=c(Vp.posterior.mean.out,Vp.posterior.mean)
  Vp.posterior.median.out=c(Vp.posterior.median.out,Vp.posterior.median)
  Vp.posterior.mode.out=c(Vp.posterior.mode.out,Vp.posterior.mode)
}

male.out=data.frame(Species="P.pungitius", Sex="Male", Iteration=1:Nsimulations, 
                    Vpe.mean=Vpe.posterior.mean.out, Vpe.median=Vpe.posterior.median.out, Vpe.mode=Vpe.posterior.mode.out,
                    Ve.mean=Ve.posterior.mean.out, Ve.median=Ve.posterior.median.out, Ve.mode=Ve.posterior.mode.out,
                    R.mean=R.posterior.mean.out, R.median=R.posterior.median.out, R.mode=R.posterior.mode.out,
		    Vp.mean=Vp.posterior.mean.out,Vp.median=Vp.posterior.median.out,Vp.mode=Vp.posterior.mode.out,
		    GrandMean=male.pop.grand.mean, VarianceOfParentalMeans=male.variance.of.parental.means)
#write.table(file="PpunSimulatedRepeatabilityMale.csv",male.out,col.names = T, row.names = F,sep="\t", quote = F)


##
#THE SUPPLEMENT TABLE2
##

s1=data.frame(Statistic=c("Population grand mean", "Population variance", "Variance of parent means", "Repeatability"), 
           Female=c(paste(round(as.numeric(quantile(female.pop.grand.mean,probs=c(0.025,0.975))),4),collapse="-"),
                    paste(round(as.numeric(quantile(female.out$Vp.mode,probs=c(0.025,0.975))),4),collapse="-"), 
                    paste(round(as.numeric(quantile(female.variance.of.parental.means,probs=c(0.025,0.975))),4),collapse="-"),
                    paste(round(as.numeric(quantile(female.out$R.mode,probs=c(0.025,0.975))),4),collapse="-")),
           
           Male=c(paste(round(as.numeric(quantile(male.pop.grand.mean,probs=c(0.025,0.975))),4),collapse="-"),
                  paste(round(as.numeric(quantile(male.out$Vp.mode,probs=c(0.025,0.975))),4),collapse="-"), 
                  paste(round(as.numeric(quantile(male.variance.of.parental.means,probs=c(0.025,0.975))),4),collapse="-"),
                  paste(round(as.numeric(quantile(male.out$R.mode,probs=c(0.025,0.975))),4),collapse="-")))

#write.table(file="TableS1_PpunSimulatedRepeatability.csv",s1,col.names = T, row.names = F,sep="\t", quote = F)




