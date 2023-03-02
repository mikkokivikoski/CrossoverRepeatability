IN = read.table("Ppungitius_Helsinki_linkagemap.txt",sep="\t",stringsAsFactors=F,skip=3,header=F)#From Kivikoski et al. 2022
header = read.table("Ppungitius_Helsinki_linkagemap.txt",sep="\t",stringsAsFactors=F,nrows=3,header=F) #From Kivikoski et al. 2022
table.out=data.frame()
raw.offspring = data.frame()

lgs=paste0("LG",1:21)
for(LG in lgs){
  d=subset(IN,V1==LG)
  for(i in 6:ncol(d)){
    MALE = unlist(strsplit(header[1,i],split="f"))[2]
    FEMALE = paste0(unlist(strsplit(header[1,i],"-"))[1],"-",unlist(strsplit(header[1,i],"-"))[2])
    FAMILY= header[1,i]
    OFFSPRING = header[2,i]
    SEX = header[3,i]
    
    tmp = data.frame(d$V1,d$V2,d$V4,d$V5,unlist(sapply(strsplit(d[,i]," "),function(x) as.numeric(x[1]))),unlist(sapply(strsplit(d[,i]," "),function(x) as.numeric(x[2]))))
    colnames(tmp) = c("CHR","SITE","PATERNALMAP","MATERNALMAP","PATERNAL","MATERNAL")
    
    sites.paternal = c()
    precise.sites.paternal = c()
    precise.sites.paternal.CM = c()
    sites.maternal = c()
    precise.sites.maternal = c()
    precise.sites.maternal.CM = c()
    intervals.paternal = c()
    intervals.maternal = c()
    
    first.paternal.non.zero=min(which(tmp$PATERNAL!=0))
    first.maternal.non.zero=min(which(tmp$MATERNAL!=0))
    
    previous.site.MATERNAL = tmp[first.maternal.non.zero,"SITE"]
    previous.site.MATERNAL.CM = tmp[first.maternal.non.zero,"MATERNALMAP"]
    previous.site.PATERNAL = tmp[first.paternal.non.zero,"SITE"]
    previous.site.PATERNAL.CM = tmp[first.paternal.non.zero,"PATERNALMAP"]
    previous.gt.PATERNAL = tmp[first.paternal.non.zero,"PATERNAL"]
    previous.gt.MATERNAL = tmp[first.maternal.non.zero,"MATERNAL"]
    haplotypes.paternal=c(previous.gt.PATERNAL)
    haplotypes.maternal=c(previous.gt.MATERNAL)
    
    for (ROW in 1:nrow(tmp)){
      paternal.gt = tmp[ROW,"PATERNAL"]
      maternal.gt = tmp[ROW,"MATERNAL"]
      if(ROW>=first.paternal.non.zero){
        #PATERNAL                         
        if (paternal.gt == previous.gt.PATERNAL){
          previous.site.PATERNAL = tmp[ROW,"SITE"]
          previous.site.PATERNAL.CM = tmp[ROW,"PATERNALMAP"]
          
        } else if (paternal.gt != previous.gt.PATERNAL & paternal.gt != 0){
          precise.sites.paternal = c(precise.sites.paternal, floor(mean(c(previous.site.PATERNAL,tmp[ROW,"SITE"]))))
          precise.sites.paternal.CM = c(precise.sites.paternal.CM, mean(c(previous.site.PATERNAL.CM,tmp[ROW,"PATERNALMAP"])))
          previous.gt.PATERNAL = paternal.gt
          haplotypes.paternal=c(haplotypes.paternal,previous.gt.PATERNAL)
          previous.site.PATERNAL = tmp[ROW,"SITE"]
          previous.site.PATERNAL.CM = tmp[ROW,"PATERNALMAP"]
        }
      }
      if(ROW>=first.maternal.non.zero){
        #MATERNAL   
        if (maternal.gt == previous.gt.MATERNAL){
          previous.site.MATERNAL = tmp[ROW,"SITE"]
          previous.site.MATERNAL.CM = tmp[ROW,"MATERNALMAP"]
          
        } else if (maternal.gt != previous.gt.MATERNAL & maternal.gt != 0){
          precise.sites.maternal = c(precise.sites.maternal, floor(mean(c(previous.site.MATERNAL,tmp[ROW,"SITE"]))))
          precise.sites.maternal.CM = c(precise.sites.maternal.CM, mean(c(previous.site.MATERNAL.CM,tmp[ROW,"MATERNALMAP"])))
          previous.gt.MATERNAL = maternal.gt
          haplotypes.maternal=c(haplotypes.maternal,previous.gt.MATERNAL)
          previous.site.MATERNAL = tmp[ROW,"SITE"]
          previous.site.MATERNAL.CM = tmp[ROW,"MATERNALMAP"]
        }
      }
    } #FOR INDIVIDUAL ENDS
    
    raw.offspring.tmp = data.frame(OFFSPRING=OFFSPRING, MALE=MALE, FEMALE=FEMALE,CHR=LG,
                                   MATERNALCOUNT=length(precise.sites.maternal),PATERNALCOUNT=length(precise.sites.paternal),MATERNALSITES=paste(precise.sites.maternal,collapse=","),
                                   SEX=SEX,PATERNALSITES=paste(precise.sites.paternal,collapse=","), PATERNALSITESCM=paste(precise.sites.paternal.CM,collapse=","), MATERNALSITESCM=paste(precise.sites.maternal.CM,collapse=","),
                                   PATERNALHAPLOTYPES=paste(haplotypes.paternal,collapse=","),MATERNALHAPLOTYPES=paste(haplotypes.maternal,collapse=","))
    raw.offspring=rbind(raw.offspring,raw.offspring.tmp)
    
  }#FOR DATA ENDS
} #FOR LG ENDS

max.co.count=max(c(raw.offspring$MATERNALCOUNT,raw.offspring$PATERNALCOUNT))
for(i in 1:max.co.count){
  raw.offspring[[paste0("MATERNALSITE",i)]]=sapply(strsplit(as.character(raw.offspring$MATERNALSITES),split=","),function(x) x[i])
  raw.offspring[[paste0("MATERNALSITE",i,"CM")]]=sapply(strsplit(as.character(raw.offspring$MATERNALSITESCM),split=","),function(x) x[i])
  raw.offspring[[paste0("PATERNALSITE",i)]]=sapply(strsplit(as.character(raw.offspring$PATERNALSITES),split=","),function(x) x[i])
  raw.offspring[[paste0("PATERNALSITE",i,"CM")]]=sapply(strsplit(as.character(raw.offspring$PATERNALSITESCM),split=","),function(x) x[i])
}
table.out=raw.offspring

#write.table(table.out,file="Ppun_Helsinki_Crossovers.csv",col.names=T,sep="\t", row.names=F, quote=F)

autosomalCO=subset(table.out,CHR!="LG12" & !(MALE %in% c("69-m-1", "72-m-2", "88-m-2", "66-m-2")))#Exclude the sex chromosome and four male parents with one offspring 

pheno.maternal=aggregate(MATERNALCOUNT ~ OFFSPRING + MALE + FEMALE + SEX,sum, data=autosomalCO)
pheno.paternal=aggregate(PATERNALCOUNT ~ OFFSPRING + MALE + FEMALE + SEX, sum, data=autosomalCO)
pheno.long=data.frame(animal=c(pheno.maternal$FEMALE,pheno.paternal$MALE), offspring=c(pheno.maternal$OFFSPRING,pheno.paternal$OFFSPRING),
                      CO=c(pheno.maternal$MATERNALCOUNT, pheno.paternal$PATERNALCOUNT),
                      offspringsex=c(pheno.maternal$SEX,pheno.paternal$SEX),
                      parentsex=rep(c("Female","Male"),each=nrow(pheno.paternal)))
pheno.long$animal=gsub(x=pheno.long$animal, pattern="-f[0-9]",replacement = "-f")
pheno.long$animal=gsub(x=pheno.long$animal, pattern="-f[0-9]",replacement = "-f")
pheno.long[pheno.long$animal=="-NA42-m-1", "animal"]="42-m-1"
pheno.long[pheno.long$animal=="-NA42-m-2", "animal"]="42-m-2"

pheno.long$ID=pheno.long$animal

#write.table(pheno.long,file="Ppun_Helsinki_Crossovers_toAnimalModel.csv",col.names=T,sep="\t", row.names=F, quote=F)
