######################################
#The gmatrix() function is courtesy of
#Brechann McGoey, from McGoey and Stinchcombe. Evol Appl, 2021
#Additional code from GA Henry
######################################
#The gmatrix() function is courtesy of
#Brechann McGoey, from McGoey and Stinchcombe. Evol Appl, 2021
#Additional code from GA Henry

## Multivariate models for Ipomoea
library(MCMCglmm)
library(MasterBayes)
library(coda)
library(bestNormalize)
#### Functions ####
gmat.list <- function(x){list(matrix(x,nrow=5,ncol=5,byrow=F))}
Gmat_HPD <-  function (Gmat) {
  corG2 <- lapply(Gmat, as.vector)
  corG2 <- do.call(rbind,corG2)
  corG2 <- as.mcmc(corG2)
  fhpd <- function (x) { x <- as.mcmc(x); HPDinterval(x, prob=0.95) }
  hpd <- round(apply(corG2,2, fhpd ),4)
  int.res <- paste(hpd[1,],hpd[2,],sep=" , ")
  mat.int <- matrix(int.res,nrow=5,ncol=5, byrow=F)
  return(mat.int)
}
phen.sd<-function (dat){
  
  traits<-c("SM.s","GR.s","FT.s","CW.s","ASD.s")
  subtraits <- dat[traits]
  
  sds<-apply(subtraits,2,sd,na.rm=T)
  
  return(sds)
  
}
# Import data ####
#Change file path to correct location
#i.multi<-read.csv("C:/Users/yourpath/i.multifull.csv", header = T)
i.multi$Pop <- factor(i.multi$Pop, levels = c("PENN", "MARY", "HOFF", "ELLR")) # Set order of levels
#i.multi<-na.omit(i.multi) #Should be done already
i.multi$Dam<-as.factor(i.multi$Dam)
i.multi$GH_area<-as.factor(i.multi$GH_area)

#These have already been done and are saved in i.multifull.csv
#i.multi$SM.s <- scale(i.multi$SM_g)
#i.multi$GR.s <- scale(i.multi$GR)
#i.multi$FT.s <- orderNorm(i.multi$DTF)$x.t
#i.multi$CW.s <- scale(i.multi$CW)
#i.multi$ASD.s <- scale(i.multi$ASD)

#P matrix####
phen.cov <- list()
phen.cov.tab <- list()
pop.list <- levels(i.multi$Pop)
m=4
for (i in 1:m) {
  i.multi.sub.pcov <- as.matrix(subset(i.multi[i.multi$Pop==pop.list[i],], select= c("SM.s", "GR.s", "FT.s", "CW.s", "ASD.s")))
  phen.cov[[i]] <- cov(i.multi.sub.pcov, use="complete.obs")
}
names(phen.cov)<-pop.list
#For overall P matrix
i.multi.sub.pcov <- as.matrix(subset(i.multi, select= c("SM.s", "GR.s", "FT.s", "CW.s", "ASD.s")))
phen.cov[[5]] <- cov(i.multi.sub.pcov, use="complete.obs")

P<-phen.cov
# Estimate G matrices ####
#### Priors list for testing ####
priors<-list (a = (list(R = list(V1 = list(V =diag(5), nu=4.002),
                                 V2 = list(V =diag(5), nu=4.002),
                                 V3 = list(V =diag(5), nu=4.002),
                                 V4 = list(V =diag(5), nu=4.002)), 
                      G = list(G1 = list(V=P[[1]]*100/4,nu=4.002),
                               G2 = list(V=P[[2]]*100/4,nu=4.002),
                               G3 = list(V=P[[3]]*100/4,nu=4.002),
                               G4 = list(V=P[[4]]*100/4,nu=4.002),
                               G5 = list(V=P[[1]]*100/4,nu=4.002),
                               G6 = list(V=P[[2]]*100/4,nu=4.002),
                               G7 = list(V=P[[3]]*100/4,nu=4.002),
                               G8 = list(V=P[[4]]*100/4,nu=4.002)))),
             b = (list(R = list(V1 = list(V =diag(5), nu=4.002),
                                V2 = list(V =diag(5), nu=4.002),
                                V3 = list(V =diag(5), nu=4.002),
                                V4 = list(V =diag(5), nu=4.002)), 
                       G = list(G1 = list(V=P[[1]]*100/6,nu=4.002),
                                G2 = list(V=P[[2]]*100/6,nu=4.002),
                                G3 = list(V=P[[3]]*100/6,nu=4.002),
                                G4 = list(V=P[[4]]*100/6,nu=4.002),
                                G5 = list(V=P[[1]]*100/6,nu=4.002),
                                G6 = list(V=P[[2]]*100/6,nu=4.002),
                                G7 = list(V=P[[3]]*100/6,nu=4.002),
                                G8 = list(V=P[[4]]*100/6,nu=4.002)))),
             c = (list(R = list(V1 = list(V =diag(5), nu=0.002),
                                V2 = list(V =diag(5), nu=0.002),
                                V3 = list(V =diag(5), nu=0.002),
                                V4 = list(V =diag(5), nu=0.002)), 
                       G = list(G1 = list(V=P[[1]]*100/4,nu=4.002),
                                G2 = list(V=P[[2]]*100/4,nu=4.002),
                                G3 = list(V=P[[3]]*100/4,nu=4.002),
                                G4 = list(V=P[[4]]*100/4,nu=4.002),
                                G5 = list(V=P[[1]]*100/4,nu=4.002),
                                G6 = list(V=P[[2]]*100/4,nu=4.002),
                                G7 = list(V=P[[3]]*100/4,nu=4.002),
                                G8 = list(V=P[[4]]*100/4,nu=4.002)))),
             d = (list(R = list(V1 = list(V =diag(5), nu=0.002),
                                V2 = list(V =diag(5), nu=0.002),
                                V3 = list(V =diag(5), nu=0.002),
                                V4 = list(V =diag(5), nu=0.002)), 
                       G = list(G1 = list(V=P[[1]]*100/6,nu=4.002),
                                G2 = list(V=P[[2]]*100/6,nu=4.002),
                                G3 = list(V=P[[3]]*100/6,nu=4.002),
                                G4 = list(V=P[[4]]*100/6,nu=4.002),
                                G5 = list(V=P[[1]]*100/6,nu=4.002),
                                G6 = list(V=P[[2]]*100/6,nu=4.002),
                                G7 = list(V=P[[3]]*100/6,nu=4.002),
                                G8 = list(V=P[[4]]*100/6,nu=4.002)))),
             e = (list(R = list(V1 = list(V =diag(5), nu=4.002),
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
                                G8 = list(V=P[[4]]*100/8,nu=4.002)))),
             f = (list(R = list(V1 = list(V =diag(5), nu=0.002),
                                V2 = list(V =diag(5), nu=0.002),
                                V3 = list(V =diag(5), nu=0.002),
                                V4 = list(V =diag(5), nu=0.002)), 
                       G = list(G1 = list(V=P[[1]]*100/4,nu=0.002),
                                G2 = list(V=P[[2]]*100/4,nu=0.002),
                                G3 = list(V=P[[3]]*100/4,nu=0.002),
                                G4 = list(V=P[[4]]*100/4,nu=0.002),
                                G5 = list(V=P[[1]]*100/4,nu=0.002),
                                G6 = list(V=P[[2]]*100/4,nu=0.002),
                                G7 = list(V=P[[3]]*100/4,nu=0.002),
                                G8 = list(V=P[[4]]*100/4,nu=0.002)))),
             g = (list(R = list(V1 = list(V =diag(5), nu=0.002),
                                V2 = list(V =diag(5), nu=0.002),
                                V3 = list(V =diag(5), nu=0.002),
                                V4 = list(V =diag(5), nu=0.002)),
                       G = list(G1 = list(V=P[[1]]*100/6,nu=0.002),
                                G2 = list(V=P[[2]]*100/6,nu=0.002),
                                G3 = list(V=P[[3]]*100/6,nu=0.002),
                                G4 = list(V=P[[4]]*100/6,nu=0.002),
                                G5 = list(V=P[[1]]*100/6,nu=0.002),
                                G6 = list(V=P[[2]]*100/6,nu=0.002),
                                G7 = list(V=P[[3]]*100/6,nu=0.002),
                                G8 = list(V=P[[4]]*100/6,nu=0.002)))),
             h = (list(R = list(V1 = list(V =diag(5), nu=0.002),
                                V2 = list(V =diag(5), nu=0.002),
                                V3 = list(V =diag(5), nu=0.002),
                                V4 = list(V =diag(5), nu=0.002)),
                      G = list(G1 = list(V=P[[1]]*100/8,nu=0.002),
                               G2 = list(V=P[[2]]*100/8,nu=0.002),
                               G3 = list(V=P[[3]]*100/8,nu=0.002),
                               G4 = list(V=P[[4]]*100/8,nu=0.002),
                               G5 = list(V=P[[1]]*100/8,nu=0.002),
                               G6 = list(V=P[[2]]*100/8,nu=0.002),
                               G7 = list(V=P[[3]]*100/8,nu=0.002),
                               G8 = list(V=P[[4]]*100/8,nu=0.002)))),
            i = (list(R = list(V1 = list(V =diag(5), nu=0.002),
                               V2 = list(V =diag(5), nu=0.002),
                               V3 = list(V =diag(5), nu=0.002),
                               V4 = list(V =diag(5), nu=0.002)),
                      G = list(G1 = list(V =diag(5),nu=0.002),
                               G2 = list(V =diag(5),nu=0.002),
                               G3 = list(V =diag(5),nu=0.002),
                               G4 = list(V =diag(5),nu=0.002),
                               G5 = list(V =diag(5),nu=0.002),
                               G6 = list(V =diag(5),nu=0.002),
                               G7 = list(V =diag(5),nu=0.002),
                               G8 = list(V =diag(5),nu=0.002)))))

#### Prior e ended up being best ####
#See DIC summary table
priore<- list(R = list(V1 = list(V =diag(5), nu=4.002),
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

##### set variables for gmatrix() ####
#Use fewer iterations for testing or it will take a very very long time
nitt<-13000
burnin<-3000
thin<-100

#### function to iterate through the priors list to run models for all ####
gmatrix<-function(dat,priors, nitt, burnin, thin, verbose){
  
  #for loop through priors
  for(i in names(priors)){  
    
    
    mod<-MCMCglmm(formula(paste('c(SM.s,GR.s,FT.s,CW.s,ASD.s)*10 ~','trait-1',sep="")),  random=~us(at.level(Pop, "PENN"):trait):Dam + us(at.level(Pop, "MARY"):trait):Dam +
                    us(at.level(Pop, "HOFF"):trait):Dam +us(at.level(Pop, "ELLR"):trait):Dam +
                    us(at.level(Pop, "PENN"):trait):GH.area + us(at.level(Pop, "MARY"):trait):GH.area +
                    us(at.level(Pop, "HOFF"):trait):GH.area +us(at.level(Pop, "ELLR"):trait):GH.area,family=c("gaussian","gaussian","gaussian","gaussian","gaussian"),
                  data=dat,rcov=~us(at.level(Pop, "PENN"):trait):units + us(at.level(Pop, "MARY"):trait):units +
                    us(at.level(Pop, "HOFF"):trait):units +us(at.level(Pop, "ELLR"):trait):units, nitt=nitt, prior=priors[[i]],verbose=verbose,  thin=thin,  
                  burnin=burnin)
    
  
    #name will be the population, model used, prior used separated by _
    
    name=paste(i,sep="_")
    print(name)
    print(mod$DIC)
    #save all output as the name assigned above
    saveRDS(mod,paste(name,'all.RDS',sep=""))
    
    #open a pdf file to save autocorr plots to check for autocorrelation
    
    pdf(file=paste("Autocorr_",i,".pdf"),onefile=TRUE)
    
    autocorr.plot(mod$Sol, auto.layout = TRUE,lag.max=50, sub=paste("S",name,sep="_"))
    autocorr.plot(mod$VCV, auto.layout = TRUE, lag.max=50, sub=paste("V",name,sep="_"))
    graphics.off()
    
    GmcmcPenn <- lapply(apply(mod$VCV[,1:25]/100, 1, gmat.list), "[[", 1)
    GmcmcMary <- lapply(apply(mod$VCV[,26:50]/100, 1, gmat.list), "[[", 1)
    GmcmcHoff <- lapply(apply(mod$VCV[,51:75]/100, 1, gmat.list), "[[", 1)
    GmcmcEllr <- lapply(apply(mod$VCV[,76:100]/100, 1, gmat.list), "[[", 1)
    
    #save just Gmcmc
    saveRDS(GmcmcPenn,paste("Penn", name,'G.RDS',sep=""))
    saveRDS(GmcmcMary,paste("Mary", name,'G.RDS',sep=""))
    saveRDS(GmcmcHoff,paste("Hoff", name,'G.RDS',sep=""))
    saveRDS(GmcmcEllr,paste("Ellr", name,'G.RDS',sep=""))
    
    print(length(mod$VCV))
    
    #set up VCV check plots
    num.plots <- 200 #200 plots the G, and E VCVs
    my.plots <- vector(num.plots, mode='list')
    my.plots2<- vector(num.plots,mode='list')
    
    for(j in 1:num.plots){
      plot(mod$VCV[,j], main=paste("VCV_",name,j,sep=""))
      my.plots[[j]] <- recordPlot()
      
      
    }
    graphics.off()
    
    #save all the VCV plots in one pdf file
    pdf(file= paste("VCV_",name,".pdf"), onefile=TRUE)
    for (k in (1:num.plots)){
      replayPlot(my.plots[[k]])
    }  
    
    graphics.off()
    
  }
  
}


#Each population using above to generate RDS file 
P<-phen.cov
mcmctest <- list(gmatrix(dat=i.multi, priors = priors, nitt = nitt, burnin = burnin, thin = thin, verbose = TRUE))

#### Full model for after best prior is found ####
nitt = 505000
burnin = 5000
thin = 500
library(MCMCglmm)

P<-phen.cov
mcmcDiagFull<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                   random=~us(at.level(Pop, "PENN"):trait):Dam + us(at.level(Pop, "MARY"):trait):Dam +
                     us(at.level(Pop, "HOFF"):trait):Dam + us(at.level(Pop, "ELLR"):trait):Dam + 
                     us(at.level(Pop, "PENN"):trait):GH.area+ us(at.level(Pop, "MARY"):trait):GH.area +
                     us(at.level(Pop, "HOFF"):trait):GH.area +us(at.level(Pop, "ELLR"):trait):GH.area, 
                   family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), data = i.multi,
                   rcov=~us(at.level(Pop, "PENN"):trait):units + us(at.level(Pop, "MARY"):trait):units + 
                    us(at.level(Pop, "HOFF"):trait):units +us(at.level(Pop, "ELLR"):trait):units, 
                  nitt=nitt, prior=priore, verbose=F, thin=thin, burnin=burnin) 
saveRDS(mcmcDiagFull, file="mcmc.diag.full.RDS")


#Extract G matrices ####
#divided by 100 to adjust for the *10 of the response variables 
# (*10 is used to improve estimates of small variances)

Gmcmc.penn<- lapply(apply(mcmcDiagFull$VCV[,1:25]/100, 1, gmat.list), "[[", 1) 
Gmcmc.mary<- lapply(apply(mcmcDiagFull$VCV[,26:50]/100, 1, gmat.list), "[[", 1) 
Gmcmc.hoff<- lapply(apply(mcmcDiagFull$VCV[,51:75]/100, 1, gmat.list), "[[", 1) 
Gmcmc.ellr<- lapply(apply(mcmcDiagFull$VCV[,76:100]/100, 1, gmat.list), "[[", 1) 


saveRDS(Gmcmc.penn, file="Gmcmc.penn.RDS")
saveRDS(Gmcmc.mary, file="Gmcmc.mary.RDS")
saveRDS(Gmcmc.hoff, file="Gmcmc.hoff.RDS")
saveRDS(Gmcmc.ellr, file="Gmcmc.ellr.RDS")

# Calculate mean G matrix for each model ####

Gmcmc.mean.penn<- matrix(colMeans(mcmcDiagFull$VCV[,1:25]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.mary<- matrix(colMeans(mcmcDiagFull$VCV[,26:50]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.hoff<- matrix(colMeans(mcmcDiagFull$VCV[,51:75]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.ellr<- matrix(colMeans(mcmcDiagFull$VCV[,76:100]),nrow=5,ncol=5,byrow=F)/100

# Get posterior 95% HPD intervals for standardized G matrices
GHPD.penn<-Gmat_HPD(Gmcmc.penn)
GHPD.mary<-Gmat_HPD(Gmcmc.mary)
GHPD.hoff<-Gmat_HPD(Gmcmc.hoff)
GHPD.ellr<-Gmat_HPD(Gmcmc.ellr)

GHPD.all<-list(GHPD.penn,GHPD.mary,GHPD.hoff,GHPD.ellr)
write.csv(GHPD.all, file = "HPDAlldiag.csv")

### Assessing Effective Sample size and autocorr ####
sum(effectiveSize(mcmcDiagFull$VCV[,1:100])>=850)

sum(autocorr.diag(mcmcDiagFull$VCV[,1:25])>0.1 & autocorr.diag(mcmcDiagFull$VCV[,1:25])<0.9999)
sum(autocorr.diag(mcmcDiagFull$VCV[,26:50])>0.1 & autocorr.diag(mcmcDiagFull$VCV[,26:50])<0.9999)
sum(autocorr.diag(mcmcDiagFull$VCV[,51:75])>0.1 & autocorr.diag(mcmcDiagFull$VCV[,51:75])<0.9999)
sum(autocorr.diag(mcmcDiagFull$VCV[,76:100])>0.1 & autocorr.diag(mcmcDiagFull$VCV[,76:100])<0.9999)

