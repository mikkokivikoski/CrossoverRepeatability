setwd("/scratch/project_2000878/211027crossOverHeritability/ppun/230125analyses/heritability/results")
library(MCMCglmm)
library(ggplot2)

NITT=1100000


load(paste0("Both-PriorC-NITT-",NITT,"-.rds"))
load(paste0("Female-PriorC-NITT-",NITT,"-.rds"))
load(paste0("Male-PriorC-NITT-",NITT,"-.rds"))

summary(ppun.heritability.both.C)
summary(ppun.heritability.female.C)
summary(ppun.heritability.male.C)


###########

R_noVf.repeated=c()
R_wVf.repeated=c()
h2_noVf.repeated=c()
h2_wVf.repeated=c()
Vp_wVf.repeated=c()
Vp_noVf.repeated=c()
Vf.repeated=c()
va.repeated=c()
vpe.repeated=c()
vr.repeated=c()

repeated.model.random.effects=data.frame()
repeated.model.fixed.effects=data.frame()

for(i in 1:3){
  print(i)
  if(i==1){
    tmp=ppun.heritability.both.C
    s="Both"
  } else if (i==2){
    tmp=ppun.heritability.female.C
    s="Female"
  } else if (i==3){
    tmp=ppun.heritability.male.C
    s="Male"
  }
  Vp_noVf = rowSums(tmp[["VCV"]])
  vf=rep(0,2000)
  vf_mode=0
  vf_int=data.frame(lower=0,upper=0)
  vf_eff=0
  ###
  if(s=="Both"){
    compute_varpred <- function(beta, design_matrix) {
      var(as.vector(design_matrix %*% beta))
    }
    X <- tmp[["X"]]
    vf <- apply(tmp[["Sol"]], 1, compute_varpred, design_matrix = X)
    vf_mode=posterior.mode(as.mcmc(vf))
    vf_int=data.frame(HPDinterval(as.mcmc(vf)),row.names=NULL)
    vf_eff=effectiveSize(as.mcmc(vf))
  }
  Vp_wVf <- Vp_noVf + vf
  h2_noVf=tmp[["VCV"]][ , "animal"] /Vp_noVf
  R_noVf=(tmp[["VCV"]][ , "animal"] + tmp[["VCV"]][ , "ID"]) /Vp_noVf
  
  h2_wVf=tmp[["VCV"]][ , "animal"] /Vp_wVf
  R_wVf=(tmp[["VCV"]][ , "animal"] + tmp[["VCV"]][ , "ID"]) /Vp_wVf
  ###
  tmp.random=data.frame(Analysis="Prior C",Sex=s,Type="Random",Estimate=c("Vf","Va","Vpe","Ve","Vp_noVf","Vp_wVf","h2_noVf","h2_wVf","R_noVF","R_wVF"))
  tmp.random$Post.mode=as.numeric(c(vf_mode,posterior.mode(tmp[["VCV"]]),posterior.mode(as.mcmc(Vp_noVf)),posterior.mode(as.mcmc(Vp_wVf)),
                                    posterior.mode(h2_noVf),posterior.mode(h2_wVf), posterior.mode(R_noVf), posterior.mode(R_wVf)))
  tmp.random$Post.median=as.numeric(c(median(vf),apply(tmp[["VCV"]],2,median),median(Vp_noVf),median(Vp_wVf),median(h2_noVf),median(h2_wVf),
                                      median(R_noVf),median(R_wVf)))
  tmp.random$Post.mean=as.numeric(c(mean(vf),apply(tmp[["VCV"]],2,mean),mean(Vp_noVf),mean(Vp_wVf),mean(h2_noVf),mean(h2_wVf),mean(R_noVf),mean(R_wVf)))
  tmp.random=cbind(tmp.random,rbind(vf_int,data.frame(HPDinterval(tmp[["VCV"]]),row.names=NULL),
                                    data.frame(HPDinterval(as.mcmc(Vp_noVf)),row.names=NULL),data.frame(HPDinterval(as.mcmc(Vp_wVf)),row.names=NULL),
                                    data.frame(HPDinterval(h2_noVf),row.names=NULL),data.frame(HPDinterval(h2_wVf),row.names=NULL),
                                    data.frame(HPDinterval(R_noVf),row.names=NULL), data.frame(HPDinterval(R_wVf),row.names=NULL)))
  tmp.random$Eff.samp=c(vf_eff,effectiveSize(tmp[["VCV"]]),effectiveSize(Vp_noVf),effectiveSize(Vp_wVf),effectiveSize(h2_noVf),effectiveSize(h2_wVf),
                        effectiveSize(R_noVf), effectiveSize(R_wVf))
  colnames(tmp.random)[c(8,9)]=c("l.95.CI","u.95.CI")
  tmp.random$pMCMC=NA
  tmp.random$DIC=tmp$DIC
  
  r=colnames(tmp$Sol)
  variables=r[!(grepl(pattern="ID.[0-9]",x=r)) & !(grepl(pattern="animal.",x=r))]
  
  
  if(s=="Both"){
    tmp.fixed=data.frame(Analysis="Prior C",Sex=s,Type="Fixed",Estimate=c(variables,"MALE"))
    male=rowSums(tmp[["Sol"]][,c(1,2)])#Get the values for males 
    tmp.fixed$Post.mode=c(as.numeric(posterior.mode(tmp[["Sol"]][,variables])),posterior.mode(as.mcmc(male)))
    tmp.fixed$Post.median=c(as.numeric(apply(tmp[["Sol"]][,variables],2,median)),median(male))
    tmp.fixed=cbind(tmp.fixed,rbind(as.data.frame(summary(tmp)$solutions,row.names=NULL),c(mean(male),HPDinterval(as.mcmc(male)),effectiveSize(rowSums(tmp[["Sol"]][,c(1,2)])),NA)))
    
  } else{
    tmp.fixed=data.frame(Analysis="Prior C",Sex=s,Type="Fixed",Estimate=variables)
    tmp.fixed$Post.mode=as.numeric(posterior.mode(tmp[["Sol"]][,variables]))
    tmp.fixed$Post.median=as.numeric(median(tmp[["Sol"]][,variables]))
    tmp.fixed=cbind(tmp.fixed,as.data.frame(summary(tmp)$solutions,row.names=NULL))
  }
  tmp.fixed$DIC=tmp$DIC
  colnames(tmp.fixed)=colnames(tmp.random)
  
  repeated.model.random.effects=rbind(repeated.model.random.effects,tmp.random)
  repeated.model.fixed.effects=rbind(repeated.model.fixed.effects,tmp.fixed)
  
  R_noVf.repeated=c(R_noVf.repeated,R_noVf)
  R_wVf.repeated=c(R_wVf.repeated,R_wVf)
  h2_noVf.repeated=c(h2_noVf.repeated,h2_noVf)
  h2_wVf.repeated=c(h2_wVf.repeated,h2_wVf)
  Vp_wVf.repeated=c(Vp_wVf.repeated,Vp_wVf)
  Vp_noVf.repeated=c(Vp_noVf.repeated,Vp_noVf)
  
  Vf.repeated=c(Vf.repeated,vf)
  va.repeated=c(va.repeated,tmp[["VCV"]][ , "animal"])
  vpe.repeated=c(vpe.repeated,tmp[["VCV"]][ , "ID"])
  vr.repeated=c(vr.repeated,tmp[["VCV"]][ , "units"])
}


plot.table=data.frame(Sex=rep(c("Both","Female","Male"),each=2000),Value=c(R_noVf.repeated,R_wVf.repeated, h2_noVf.repeated,h2_wVf.repeated, va.repeated, vpe.repeated, vr.repeated, Vp_noVf.repeated), Estimate=rep(c("Repeatability_noVf","Repeatability_wVf", "Heritability_noVf", "Heritability_wVf","Va","Vpe","Vr","Vp"),each=6000))
plot.table$Estimate=factor(plot.table$Estimate,levels=c("Repeatability_noVf","Repeatability_wVf","Heritability_noVf","Heritability_wVf","Va","Vpe","Vr","Vp"),labels=c(expression(italic(R)),expression(italic(R)[italic("Vf")]), expression(italic("h")^italic("2")),expression(italic("h")[italic("Vf")]^italic("2")),expression(italic("V")[italic("a")]), expression(italic("V")[italic("pe")]), expression(italic("V")[italic("r")]), expression(italic("V")[italic("p")])))
plot.table.S1=subset(plot.table, Sex=="Both" & Estimate !="italic(\"V\")[italic(\"p\")]")
plot.table.S1=rbind(plot.table.S1, subset(plot.table, Sex=="Female" & Estimate %in% c("italic(R)","italic(\"h\")^italic(\"2\")","italic(\"V\")[italic(\"a\")]", "italic(\"V\")[italic(\"pe\")]","italic(\"V\")[italic(\"r\")]")))
plot.table.S1=rbind(plot.table.S1, subset(plot.table, Sex=="Male" & Estimate %in% c("italic(R)","italic(\"h\")^italic(\"2\")","italic(\"V\")[italic(\"a\")]", "italic(\"V\")[italic(\"pe\")]","italic(\"V\")[italic(\"r\")]")))



ggplot(plot.table.S1, aes(Value,after_stat(ndensity)))+ geom_freqpoly() + 
  facet_grid(Sex ~ Estimate,scales="free",labeller=label_parsed) + theme_bw() + theme(panel.grid.minor = element_blank(),axis.text.x =element_text(size=5)) + 
  ylab("Normalized density")+ xlab("Estimate")
ggsave(paste0("FigS1_PpunPosterior_",NITT,"-PriorC.pdf"),width = 120,height=80,units="mm")
#Va-Vpe correlation
plot.table2=data.frame(Sex=rep(c("Both","Female","Male"),each=2000),Va=va.repeated, Vpe=vpe.repeated, Species="P. pungigitus")
ggplot(plot.table2, aes(x=Vpe,y=Va,))+ geom_point(alpha=0.1) + facet_wrap(~ Sex,scales="free") + theme_bw()+ theme(panel.grid.minor = element_blank()) + 
  ylab(expression(italic("V")[italic("a")]))+ xlab(expression(italic("V")[italic("pe")]))
write.table(plot.table2,file=paste0("FigS3_PpunVaVpe_",NITT,"-PriorC.csv"),sep="\t",row.names=F,col.names=T,quote=F)
ggsave(paste0("FigS3_PpunVaVpe_",NITT,"-PriorC-.pdf"),width = 120,height=80,units="mm")

#Write result tables

write.table(rbind(repeated.model.fixed.effects,repeated.model.random.effects),file=paste0(NITT,"-ppunHeritabilityEstimates-PriorC.csv"),sep="\t",row.names=F,col.names=T,quote=F)

rounded.out=rbind(repeated.model.fixed.effects,repeated.model.random.effects)
rounded.out[,5:9]=round(rounded.out[,5:9],4)
write.table(rounded.out,file=paste0(NITT,"-ppunHeritabilityEstimatesRounded-PriorC.csv"),sep="\t",row.names=F,col.names=T,quote=F)

############
