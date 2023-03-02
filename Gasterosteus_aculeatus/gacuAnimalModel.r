library(snpReady)
library(MCMCglmm)

#######
load("gacuGRM.rds")#snpreadyGRM
pheno = read.table('../../data/gacuCrossovers.csv',header = TRUE,sep="\t",stringsAsFactor=F)
pheno.maternal=aggregate(MATERNALCOUNT ~ OFFSPRING + MALE + FEMALE + SEX,sum, data=subset(pheno,CHR!="chrXIX"))#Exclude the sex chromsome 19
pheno.paternal=aggregate(PATERNALCOUNT ~ OFFSPRING + MALE + FEMALE + SEX, sum, data=subset(pheno,CHR!="chrXIX"))#Exclude the sex chromsome 19
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


BURN <- 100000; THIN <- 500
NITT = 1100000

prior.C <- list(R = list(V = 1, nu = 0.002),
                G = list(G1 = list(V = 1, nu = 0.002),
                         G2 = list(V = 1, nu = 0.002)))


gacu.heritability.both.C <- MCMCglmm(CO ~ 1 + parentsex,
                                   random = ~ animal+ID,
                                   rcov = ~ units,
                                   ginverse = list(animal = snpreadyGRM),
                                   data = data.pheno,
                                   prior = prior.C,
                                   family = c("gaussian"),
                                   nitt = NITT, thin = THIN, burnin = BURN,verbose=F)                 

save(gacu.heritability.both.C,file=paste("results/PriorC-Sex-NITT",NITT,".rds",sep="-"))

gacu.heritability.female.C <- MCMCglmm(CO ~ 1 ,
                                       random = ~ animal+ID,
                                       rcov = ~ units,
                                       ginverse = list(animal = snpreadyGRM),
                                       data = fem,
                                       prior = prior.C,
                                       family = c("gaussian"),
                                       nitt = NITT, thin = THIN, burnin = BURN,verbose=F)                 

save(gacu.heritability.female.C,file=paste("results/Female-PriorC-NITT",NITT,".rds",sep="-"))

gacu.heritability.male.C <- MCMCglmm(CO ~ 1 ,
                                   random = ~ animal+ID,
                                   rcov = ~ units,
                                   ginverse = list(animal = snpreadyGRM),
                                   data = male,
                                   prior = prior.C,
                                   family = c("gaussian"),
                                   nitt = NITT, thin = THIN, burnin = BURN,verbose=F)                 

save(gacu.heritability.male.C,file=paste("results/Male-PriorC-NITT",NITT,".rds",sep="-"))

