if(!require(MCMCglmm)){
  install.packages("MCMCglmm", repos = "https://mirror.rcg.sfu.ca/mirror/CRAN/")
  library(MCMCglmm)
}
####Subset populations####
#Change file path to wherever i.multifull.csv is
i.multi <- read.csv("/scratch/henrygeo/RNAproject/Scripts/i.multifull.csv", header = T)
i.multi$Dam<-as.factor(i.multi$Dam)
i.multi$GH.area<-as.factor(i.multi$GH.area)
i.multi$Pop<-factor(i.multi$Pop, levels = c("PENN", "MARY", "HOFF", "ELLR"))
i.multi.all <- as.data.frame(subset(i.multi, select= c("Dam", "Pop", "GH.area", "SM.s", "GR.s", "FT.s", "CW.s", "ASD.s")))

# G matrix for PENNSYLVANIA population 
i.multi.penn <- subset(i.multi.all, Pop=="PENN")
i.multi.penn <- droplevels(i.multi.penn)

# G matrix for MARYLAND population 
i.multi.mary <- subset(i.multi.all, Pop=="MARY")
i.multi.mary <- droplevels(i.multi.mary)


# G matrix for HOFFMAN, NC population
i.multi.hoff <- subset(i.multi.all, Pop=="HOFF")
i.multi.hoff <- droplevels(i.multi.hoff)


# G matrix for ELLERBE, NC population 
i.multi.ellr <- subset(i.multi.all, Pop=="ELLR")
i.multi.ellr <- droplevels(i.multi.ellr)

#### P matrix ####
phen.cov <- list()
phen.cov.tab <- list()
pop.list <- levels(i.multi$Pop)
for (i in 1:4) {
  i.multi.sub.pcov <- as.matrix(subset(i.multi[i.multi$Pop==pop.list[i],], select= c("SM.s", "GR.s", "FT.s", "CW.s", "ASD.s")))
  phen.cov[[i]] <- cov(i.multi.sub.pcov, use="complete.obs")
}
names(phen.cov)<-pop.list
i.multi.sub.pcov <- as.matrix(subset(i.multi, select= c("SM.s", "GR.s", "FT.s", "CW.s", "ASD.s")))
phen.cov[[5]] <- cov(i.multi.sub.pcov, use="complete.obs")
names(phen.cov)[[5]]<-"all"

#### Prior and parameters ####
P<-phen.cov
#Same prior as for the observed G matrices
priorse<- list(R = list(V1 = list(V =diag(5), nu=4.002),
                   V2 = list(V =diag(5), nu=4.002),
                   V3 = list(V =diag(5), nu=4.002),
                   V4 = list(V =diag(5), nu=4.002)), 
          G = list(G1 = list(V=P[[1]]*100/8,nu=4.002),
                   G2 = list(V=P[[2]]*100/8,nu=4.002),
                   G3 = list(V=P[[3]]*100/8,nu=4.002),
                   G4 = list(V=P[[4]]*100/8,nu=4.002),
                   G5 = list(V=P[[1]]*100/8,nu=4.002),
                   G6 = list(V=P[[2]]*100/8,nu=4.002),
                   G7 = list(V=P[[3]]*100/8,nu=4.002),
                   G8 = list(V=P[[4]]*100/8,nu=4.002)))
nitt<-20000
burnin<-5000
thin<-100

#### Loop to run randomized within pop G models ####
for (i in 1:1000){
  
#shuffle individuals across dams
  i.multi.penn$randDam <- sample(i.multi.penn$Dam, length(i.multi.penn$Dam), replace = F)
  i.multi.mary$randDam <- sample(i.multi.mary$Dam, length(i.multi.mary$Dam), replace = F)
  i.multi.hoff$randDam <- sample(i.multi.hoff$Dam, length(i.multi.hoff$Dam), replace = F)
  i.multi.ellr$randDam <- sample(i.multi.ellr$Dam, length(i.multi.ellr$Dam), replace = F)
  i.multi.rand <- data.frame(rbind(i.multi.penn, i.multi.mary, i.multi.hoff, i.multi.ellr))
  mcmcdiag<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                     random=~us(at.level(Pop, "PENN"):trait):randDam + us(at.level(Pop, "MARY"):trait):randDam +
                       us(at.level(Pop, "HOFF"):trait):randDam +us(at.level(Pop, "ELLR"):trait):randDam +
                       us(at.level(Pop, "PENN"):trait):GH.area + us(at.level(Pop, "MARY"):trait):GH.area +
                       us(at.level(Pop, "HOFF"):trait):GH.area +us(at.level(Pop, "ELLR"):trait):GH.area,family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), data=i.multi.rand,
                     rcov=~us(at.level(Pop, "PENN"):trait):units + us(at.level(Pop, "MARY"):trait):units +
                       us(at.level(Pop, "HOFF"):trait):units +us(at.level(Pop, "ELLR"):trait):units, nitt=nitt, prior=priorse, verbose=T, thin=thin, burnin=burnin) 
  saveRDS(mcmcdiag, file = paste("MCMCoutput/Diag",i,"poprand.RDS",sep="_"))
}  
