###################################################
### RECOMBINATION RATE REPEATABILITY SIMULATION ###
###################################################
library(MCMCglmm)
sessionInfo()

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
pheno = read.table('../../data/gacuCrossovers.csv',header = TRUE,sep="\t",stringsAsFactor=F)

pheno.maternal=aggregate(MATERNALCOUNT ~ OFFSPRING + MALE + FEMALE + SEX,sum, data=subset(pheno,CHR!="chrXIX"))
pheno.paternal=aggregate(PATERNALCOUNT ~ OFFSPRING + MALE + FEMALE + SEX, sum, data=subset(pheno,CHR!="chrXIX"))
pheno.long=data.frame(animal=c(pheno.maternal$FEMALE,pheno.paternal$MALE), offspring=c(pheno.maternal$OFFSPRING,pheno.paternal$OFFSPRING),
                      CO=c(pheno.maternal$MATERNALCOUNT, pheno.paternal$PATERNALCOUNT),
                      offspringsex=c(pheno.maternal$SEX,pheno.paternal$SEX),
                      parentsex=rep(c("Female","Male"),each=nrow(pheno.paternal)))
pheno.long$ID=pheno.long$animal
ggplot(pheno.long,aes(x=reorder(animal, CO, FUN=mean), y=CO, col=parentsex)) + geom_boxplot()

fem=subset(pheno.long,parentsex=="Female")
fem=droplevels(fem)

male=subset(pheno.long,parentsex=="Male")
male=droplevels(male)

#Maximum-likelihood estimates
ml.paternal=read.table("",header=T,stringsAsFactors = F,sep="\t")#Kivikoski et al. 2022 supplementary table 4
ml.maternal=read.table("",header=T,stringsAsFactors = F,sep="\t")#Kivikoski et al. 2022 supplementary table 3


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
  crossovers=matrix(nrow=21,ncol=517)
  for (r in 1:21){
    crossovers[r,]=replicate(expr=generate.gamete.yf(as.numeric(ml.maternal[r,paste0("MLEpRestricted",0:9)])),n=517)
  }
  fem[,paste0("simulatedOffspring",i)]=colSums(crossovers[-19,])
  fem$tmp=colSums(crossovers[-19,])
  var.parent.means=var(aggregate(tmp ~ animal,mean,data=fem)$tmp)
  female.pop.grand.mean=c(female.pop.grand.mean,mean(mean(fem$tmp)))
  female.variance.of.parental.means=c(female.variance.of.parental.means,var.parent.means)
  female.population.variance=c(female.population.variance,var(fem$tmp))
}

#hist(female.variance.of.parental.means)
#hist(female.population.variance)


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
  gacu.heritability.female.ranomized <- MCMCglmm(COtmp ~ 1,
                                                 random = ~ ID,
                                                 rcov = ~ units,
                                                 data = fem,
                                                 prior = prior,
                                                 family = c("gaussian"),pl=T,pr=T,
                                                 nitt = NITT, thin = THIN, burnin = BURN,verbose=F)

  
  Vpe=gacu.heritability.female.ranomized[["VCV"]][,"ID"]
  Vpe.posterior.mean=mean(Vpe)
  Vpe.posterior.median=median(Vpe)
  Vpe.posterior.mode=posterior.mode(Vpe)

  Ve=gacu.heritability.female.ranomized[["VCV"]][,"units"]
  Ve.posterior.mean=mean(Ve)
  Ve.posterior.median=median(Ve)
  Ve.posterior.mode=posterior.mode(Ve)

  R=gacu.heritability.female.ranomized[["VCV"]][,"ID"]/rowSums(gacu.heritability.female.ranomized[["VCV"]])
  R.posterior.mean=mean(R)
  R.posterior.median=median(R)
  R.posterior.mode=posterior.mode(R)

  Vp=rowSums(gacu.heritability.female.ranomized[["VCV"]])
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
female.out=data.frame(Species="G.aculeatus",Sex="Female",Iteration=1:Nsimulations,
                    #Va.mean=Va.posterior.mean.out, Va.median=Va.posterior.median.out, Va.mode=Va.posterior.mode.out,
                    Vpe.mean=Vpe.posterior.mean.out, Vpe.median=Vpe.posterior.median.out, Vpe.mode=Vpe.posterior.mode.out,
                    Ve.mean=Ve.posterior.mean.out, Ve.median=Ve.posterior.median.out, Ve.mode=Ve.posterior.mode.out,
                    R.mean=R.posterior.mean.out, R.median=R.posterior.median.out, R.mode=R.posterior.mode.out,
                    Vp.mean=Vp.posterior.mean.out,Vp.median=Vp.posterior.median.out,Vp.mode=Vp.posterior.mode.out,
                    GrandMean=female.pop.grand.mean, VarianceOfParentalMeans=female.variance.of.parental.means)

#write.table(file="GacuSimulatedRepeatabilityFemale.csv",female.out,col.names = T, row.names = F,sep="\t", quote = F)




##SIMULATE MALE DATA
male.pop.grand.mean=c()
male.variance.of.parental.means=c()
male.population.variance=c()
Nsimulations=1000
for(i in 1:Nsimulations){
  crossovers=matrix(nrow=21,ncol=517)
  for (r in 1:21){
    crossovers[r,]=replicate(expr=generate.gamete.yf(as.numeric(ml.paternal[r,paste0("MLEpRestricted",0:9)])),n=517)
  }
  male[,paste0("simulatedOffspring",i)]=colSums(crossovers[-19,])
  male$tmp=colSums(crossovers[-19,])
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
  gacu.heritability.male.ranomized <- MCMCglmm(COtmp ~ 1,
                                                 random = ~ ID,
                                                 rcov = ~ units,
                                                 data = male,
                                                 prior = prior,
                                                 family = c("gaussian"),pl=T,pr=T,
                                                 nitt = NITT, thin = THIN, burnin = BURN,verbose=F)

  Vpe=gacu.heritability.male.ranomized[["VCV"]][,"ID"]
  Vpe.posterior.mean=mean(Vpe)
  Vpe.posterior.median=median(Vpe)
  Vpe.posterior.mode=posterior.mode(Vpe)
  
  Ve=gacu.heritability.male.ranomized[["VCV"]][,"units"]
  Ve.posterior.mean=mean(Ve)
  Ve.posterior.median=median(Ve)
  Ve.posterior.mode=posterior.mode(Ve)

  R=(gacu.heritability.male.ranomized[["VCV"]][,"ID"])/rowSums(gacu.heritability.male.ranomized[["VCV"]])
  R.posterior.mean=mean(R)
  R.posterior.median=median(R)
  R.posterior.mode=posterior.mode(R)
  
  Vp=rowSums(gacu.heritability.male.ranomized[["VCV"]])
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

male.out=data.frame(Species="G.aculeatus", Sex="Male", Iteration=1:Nsimulations, 
                    Vpe.mean=Vpe.posterior.mean.out, Vpe.median=Vpe.posterior.median.out, Vpe.mode=Vpe.posterior.mode.out,
                    Ve.mean=Ve.posterior.mean.out, Ve.median=Ve.posterior.median.out, Ve.mode=Ve.posterior.mode.out,
                    R.mean=R.posterior.mean.out, R.median=R.posterior.median.out, R.mode=R.posterior.mode.out,
		    Vp.mean=Vp.posterior.mean.out,Vp.median=Vp.posterior.median.out,Vp.mode=Vp.posterior.mode.out,
		    GrandMean=male.pop.grand.mean, VarianceOfParentalMeans=male.variance.of.parental.means)
#write.table(file="GacuSimulatedRepeatabilityMale.csv",male.out,col.names = T, row.names = F,sep="\t", quote = F)


##
#THE SUPPLEMENT TABLE2
##

s2=data.frame(Statistic=c("Population grand mean", "Population variance", "Variance of parent means", "Repeatability"), 
           Female=c(paste(round(as.numeric(quantile(female.pop.grand.mean,probs=c(0.025,0.975))),4),collapse="-"),
                    paste(round(as.numeric(quantile(female.out$Vp.mode,probs=c(0.025,0.975))),4),collapse="-"), 
                    paste(round(as.numeric(quantile(female.variance.of.parental.means,probs=c(0.025,0.975))),4),collapse="-"),
                    paste(round(as.numeric(quantile(female.out$R.mode,probs=c(0.025,0.975))),4),collapse="-")),
           
           Male=c(paste(round(as.numeric(quantile(male.pop.grand.mean,probs=c(0.025,0.975))),4),collapse="-"),
                  paste(round(as.numeric(quantile(male.out$Vp.mode,probs=c(0.025,0.975))),4),collapse="-"), 
                  paste(round(as.numeric(quantile(male.variance.of.parental.means,probs=c(0.025,0.975))),4),collapse="-"),
                  paste(round(as.numeric(quantile(male.out$R.mode,probs=c(0.025,0.975))),4),collapse="-")))

#write.table(file="TableS2_GacuSimulatedRepeatability.csv",s2,col.names = T, row.names = F,sep="\t", quote = F)




