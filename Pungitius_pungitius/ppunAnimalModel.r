library(MCMCglmm)

#######
data.pheno=read.table("Ppun_Helsinki_Crossovers_toAnimalModel.csv", header=T, stringsAsFactors = F,sep="\t")
load("PpunGRM.rds")

data.pheno.female=subset(data.pheno,parentsex=="Female")
data.pheno.male=subset(data.pheno,parentsex=="Male")

BURN <- 100000; THIN <- 500
NITT = 1100000



############################


prior.C <- list(R = list(V = 1, nu = 0.002),
                G = list(G1 = list(V = 1, nu = 0.002),
                         G2 = list(V = 1, nu = 0.002)))

#OVERALL HERITABILITY


ppun.heritability.both.C <- MCMCglmm(CO ~ 1 + parentsex,
                                   random = ~ animal+ID,
                                   rcov = ~ units,
                                   ginverse = list(animal = snpreadyGRM),
                                   data = data.pheno,
                                   prior = prior.C,
                                   family = c("gaussian"),
                                   nitt = NITT, thin = THIN, burnin = BURN,verbose=F)                 

save(ppun.heritability.both.C,file=paste("results/Both-PriorC-NITT",NITT,".rds",sep="-"))



ppun.heritability.female.C <- MCMCglmm(CO ~ 1 ,
                                       random = ~ animal+ID,
                                       rcov = ~ units,
                                       ginverse = list(animal = snpreadyGRM),
                                       data = data.pheno.female,
                                       prior = prior.C,
                                       family = c("gaussian"),
                                       nitt = NITT, thin = THIN, burnin = BURN,verbose=F)                 

save(ppun.heritability.female.C,file=paste("results/Female-PriorC-NITT",NITT,".rds",sep="-"))


ppun.heritability.male.C <- MCMCglmm(CO ~ 1 ,
                                   random = ~ animal+ID,
                                   rcov = ~ units,
                                   ginverse = list(animal = snpreadyGRM),
                                   data = data.pheno.male,
                                   prior = prior.C,
                                   family = c("gaussian"),
                                   nitt = NITT, thin = THIN, burnin = BURN,verbose=F)                 

save(ppun.heritability.male.C,file=paste("results/Male-PriorC-NITT",NITT,".rds",sep="-"))

