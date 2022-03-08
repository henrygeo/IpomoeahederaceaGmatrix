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
#P matrix####
phen.sd<-function (dat){
  
  traits<-c("SM.s","GR.s","FT.s","CW.s","ASD.s")
  subtraits <- dat[traits]
  
  sds<-apply(subtraits,2,sd,na.rm=T)
  
  return(sds)
  
}
phen.cov <- list()
phen.cov.tab <- list()
pop.list <- levels(i.multi$Pop)
m=4
for (i in 1:m) {
  i.multi.sub.pcov <- as.matrix(subset(i.multi[i.multi$Pop==pop.list[i],], select= c("SM.s", "GR.s", "FT.s", "CW.s", "ASD.s")))
  phen.cov[[i]] <- cov(i.multi.sub.pcov, use="complete.obs")
}
names(phen.cov)<-pop.list

# Import data ####
#Change file path to correct location
#i.multi<-read.csv("C:/Users/yourpath/i.multifull.csv", header = T)
i.multi$Pop <- factor(i.multi$Pop, levels = c("PENN", "MARY", "HOFF", "ELLR")) # Set order of levels
#i.multi<-na.omit(i.multi) #Should be done already
i.multi$Dam<-as.factor(i.multi$Dam)

#These have already been done in i.multifull.csv
#i.multi$SM.s <- scale(i.multi$SM_g)
#i.multi$GR.s <- scale(i.multi$GR)
#i.multi$FT.s <- orderNorm(i.multi$DTF)$x.t
#i.multi$CW.s <- scale(i.multi$CW)
#i.multi$ASD.s <- scale(i.multi$ASD)


#### Subset populations ####
# G matrix for PENNSYLVANIA population 
i.multi.penn <- subset(i.multi, Pop=="PENN")
i.multi.penn <- droplevels(i.multi.penn)

# G matrix for MARYLAND population 
i.multi.mary <- subset(i.multi, Pop=="MARY")
i.multi.mary <- droplevels(i.multi.mary)


# G matrix for HOFFMAN, NC population
i.multi.hoff <- subset(i.multi, Pop=="HOFF")
i.multi.hoff <- droplevels(i.multi.hoff)


# G matrix for ELLERBY, NC population 
i.multi.ellr <- subset(i.multi, Pop=="ELLR")
i.multi.ellr <- droplevels(i.multi.ellr)

# Estimate G matrices ####
#### Priors list for testing ####
priors<-list (a=(list(R=list(V=diag(5), nu=4.001), 
                      G = list(G1 = list(V=diag(5),nu=4.001),
                               G2 = list(V=diag(5),nu=4.001)
                      ))),b=(list(R=list(V=diag(5), nu=0.002), 
                                  G = list(G1 = list(V=diag(5),nu=0.002),
                                           G2 = list(V=diag(5),nu=0.002)
                                  ))),c= (list(R=list(V=diag(5), nu=4.001), 
                                               G = list(G1 = list(V=P/2,nu=4.001),
                                                        G2 = list(V=P/2,nu=4.001)
                                               ))),d= (list(R=list(V=diag(5), nu=0.002), 
                                                            G = list(G1 = list(V=P*10/2,nu=0.002),
                                                                     G2 = list(V=P*10/2,nu=0.002)
                                                            ))), e= (list(R=list(V=diag(5), nu=4.001), 
                                                                          G = list(G1 = list(V=P*10/4,nu=4.001),
                                                                                   G2 = list(V=P*10/4,nu=4.001)
                                                                          ))),f= (list(R=list(V=diag(5), nu=0.002), 
                                                                                       G = list(G1 = list(V=P*10/4,nu=0.002),
                                                                                                G2 = list(V=P*10/4,nu=0.002)
                                                                                       ))),g= (list(R=list(V=diag(5), nu=4.001), 
                                                                                                    G = list(G1 = list(V=P*10/3,nu=4.001),
                                                                                                             G2 = list(V=P*10/3,nu=4.001)
                                                                                                   ))),h = (list(R=list(V=diag(5), nu=4.001),
                                                                                                                G=list(G1=list(V=P*10/4, nu=0.002),
                                                                                                                       G2=list(V=diag(5), nu=0.002))))
                                                                                       )


#### Prior g ended up being best ####
#See DIC summary table
priorsg<- list(R=list(V=diag(5), nu=4.001), 
                        G = list(G1 = list(V=P*10/3,nu=4.001),
                                 G2 = list(V=P*10/3,nu=4.001)))


##### set variables for gmatrix() ####
#Use fewer iterations for testing or it will take a very very long time
nitt<-1003000
burnin<-3000
thin<-100


pop<-list(i.multi.penn, i.multi.mary, i.multi.hoff, i.multi.ellr)
names(pop)<-c("Penn", "Mary", "Hoff", "Ellr")
model1='trait-1'

#function to estimate G matrix
#arguments are dat --> the data
#             model --> can dictate whether to include fixed effects or not
#             priors --> needs to be a list
#             nitt --> number of itterations
#             burnin --> number from beginning to be thrown away
#             thin -->intervals at which the Markov chain is stored
#             verbose --> print progress statements

#### function to iterate through the priors list to run models for all ####
gmatrix<-function(dat,model,priors, nitt, burnin, thin, verbose){
  
  list1<-list()
   
  print(model)
  popname<-levels(dat[[2]])


  #for loop through priors
  for(i in names(priors)){  
    
    
    mod<-MCMCglmm(formula(paste('c(SM.s,GR.s,FT.s,CW.s,ASD.s)*10 ~',model,sep="")),  random=~us(trait):Dam + us(trait):GH.area,  
                  family=c("gaussian","gaussian","gaussian","gaussian","gaussian"),
                  data=dat,rcov=~us(trait):units, nitt=nitt, prior=priors[[i]],verbose=verbose,  thin=thin,  
                  burnin=burnin)
    
  
    #name will be the population, model used, prior used separated by _
    
    name=paste(popname,model,i,sep="_")
    print(name)
    print(mod$DIC)
    #save all output as the name assigned above
    saveRDS(mod,paste(name,'all.RDS',sep=""))
    
    #open a pdf file to save autocorr plots to check for autocorrelation
    pdf(file=paste("Autocorr_",popname,i,".pdf"),onefile=TRUE)
    
    autocorr.plot(mod$Sol, auto.layout = TRUE,lag.max=50, sub=paste("S",name,sep="_"))
    autocorr.plot(mod$VCV, auto.layout = TRUE, lag.max=50, sub=paste("V",name,sep="_"))
    graphics.off()
    
    
    Gmcmc <- lapply(apply(mod$VCV[,1:25]/100, 1, f), "[[", 1)
    #save jut Gmcmc
    saveRDS(Gmcmc,paste(name,'G.RDS',sep=""))
    Gmcmc.mean <- matrix(colMeans(mod$VCV[,1:25]),nrow=5,ncol=5,byrow=F)/100
    GmatHPD1<-Gmat_HPD(Gmcmc)
    print(length(mod$VCV))
    
    
    #set up VCV check plots
    num.plots <- 50
    my.plots <- vector(num.plots, mode='list')
    my.plots2<- vector(num.plots,mode='list')
    
    for(j in 1:50){
      plot(mod$VCV[,j], main=paste("VCV_",name,j,sep=""))
      my.plots[[j]] <- recordPlot()
      
      
    }
    graphics.off()
    
    #save all the VCV plots in one pdf file
    pdf(file= paste("VCV_",name,".pdf"), onefile=TRUE)
    for (k in (1:50)){
      replayPlot(my.plots[[k]])
    }  
    
    graphics.off()
    
  }
  
}


#Each population using above to generate RDS file 
P<-phen.cov[["PENN"]]
mcmcPenntest<-list(gmatrix(dat=i.multi.penn,model=model1,priors=priors,nitt=nitt,burnin=burnin,thin=thin,verbose=F))

P<-phen.cov[["MARY"]]
mcmcMarytest<-list(gmatrix(dat=i.multi.mary,model=model1,priors=priors,nitt=nitt,burnin=burnin,thin=thin,verbose=F))                        

P<-phen.cov[["HOFF"]]
mcmcHofftest<-list(gmatrix(dat=i.multi.hoff,model=model1,priors=priors,nitt=nitt,burnin=burnin,thin=thin,verbose=F))

P<-phen.cov[["ELLR"]]
mcmcEllrtest<-list(gmatrix(dat=i.multi.ellr,model=model1,priors=priors,nitt=nitt,burnin=burnin,thin=thin,verbose=F))

#### Full model for after best prior is found ####
library(MCMCglmm)
P<-phen.cov[["PENN"]]
mcmcPenn<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
         random=~us(trait):Dam + us(trait):GH.area,family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), data=i.multi.penn,
         rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
saveRDS(mcmcPenn, file="mcmc.penn.full.RDS")

P<-phen.cov[["MARY"]]
mcmcMary<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                   random=~us(trait):Dam + us(trait):GH.area,family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), data=i.multi.mary,
                   rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
saveRDS(mcmcMary, file="mcmc.mary.full.RDS")

P<-phen.cov[["HOFF"]]
mcmcHoff<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                   random=~us(trait):Dam + us(trait):GH.area,family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), data=i.multi.hoff,
                   rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
saveRDS(mcmcHoff, file="mcmc.hoff.full.RDS")

P<-phen.cov[["ELLR"]]
mcmcEllr<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                   random=~us(trait):Dam + us(trait):GH.area,family=c("gaussian","gaussian","gaussian","gaussian","gaussian"), data=i.multi.ellr,
                   rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
saveRDS(mcmcEllr, file="mcmc.ellr.full.RDS")


## Extract VCV matrices ####
#for 5 traits 1:25 is G, 26:50 is env (GH.area), 51:75 is error
VCV.penn<- mcmcPenn$VCV[,1:75]/100
VCV.mary<- mcmcMary$VCV[,1:75]/100
VCV.hoff<- mcmcHoff$VCV[,1:75]/100
VCV.ellr<- mcmcEllr$VCV[,1:75]/100
saveRDS(VCV.penn, file="VCV.penn.RDS")
saveRDS(VCV.mary, file="VCV.mary.RDS")
saveRDS(VCV.hoff, file="VCV.hoff.RDS")
saveRDS(VCV.ellr, file="VCV.ellr.RDS")

#Extract G matrices ####
#divided by 100 to adjust for the *10 of the response variables 
# (*10 is used to improve estimates of small variances)

Gmcmc.penn<- lapply(apply(mcmcPenn$VCV[,1:25]/100, 1, gmat.list), "[[", 1) 
Gmcmc.mary<- lapply(apply(mcmcMary$VCV[,1:25]/100, 1, gmat.list), "[[", 1) 
Gmcmc.hoff<- lapply(apply(mcmcHoff$VCV[,1:25]/100, 1, gmat.list), "[[", 1) 
Gmcmc.ellr<- lapply(apply(mcmcEllr$VCV[,1:25]/100, 1, gmat.list), "[[", 1) 

saveRDS(Gmcmc.penn, file="Gmcmc.penn.RDS")
saveRDS(Gmcmc.mary, file="Gmcmc.mary.RDS")
saveRDS(Gmcmc.hoff, file="Gmcmc.hoff.RDS")
saveRDS(Gmcmc.ellr, file="Gmcmc.ellr.RDS")

# Calculate mean G matrix for each model ####
Gmcmc.mean.penn<- matrix(colMeans(mcmcPenn$VCV[,1:25]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.mary<- matrix(colMeans(mcmcMary$VCV[,1:25]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.hoff<- matrix(colMeans(mcmcHoff$VCV[,1:25]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.ellr<- matrix(colMeans(mcmcEllr$VCV[,1:25]),nrow=5,ncol=5,byrow=F)/100

# or median
Gmcmc.median.penn<- matrix(apply(mcmcPenn$VCV[,1:25],2,median),nrow=5,ncol=5,byrow=F)/100
Gmcmc.median.mary<- matrix(apply(mcmcMary$VCV[,1:25],2,median),nrow=5,ncol=5,byrow=F)/100
Gmcmc.median.hoff<- matrix(apply(mcmcHoff$VCV[,1:25],2,median),nrow=5,ncol=5,byrow=F)/100
Gmcmc.median.ellr<- matrix(apply(mcmcEllr$VCV[,1:25],2,median),nrow=5,ncol=5,byrow=F)/100

# Get posterior 95% HPD intervals for standardized G matrices
GHPD.penn<-Gmat_HPD(Gmcmc.penn)
GHPD.mary<-Gmat_HPD(Gmcmc.mary)
GHPD.hoff<-Gmat_HPD(Gmcmc.hoff)
GHPD.ellr<-Gmat_HPD(Gmcmc.ellr)

GHPD.all<-list(GHPD.penn,GHPD.mary,GHPD.hoff,GHPD.ellr)
write.csv(GHPD.all, file = "HPDAll.csv")

### Assessing Effective Sample size and autocorr ####
sum(effectiveSize(mcmc.penn.full$VCV)>8500)
sum(effectiveSize(mcmc.mary.full$VCV)>8500)
sum(effectiveSize(mcmc.hoff.full$VCV)>8500)
sum(effectiveSize(mcmc.ellr.full$VCV)>8500)

sum(autocorr.diag(mcmc.penn.full$VCV[,1:25])>0.1 & autocorr.diag(mcmc.penn.full$VCV[,1:25])<0.9999)
sum(autocorr.diag(mcmc.mary.full$VCV[,1:25])>0.1 & autocorr.diag(mcmc.mary.full$VCV[,1:25])<0.9999)
sum(autocorr.diag(mcmc.hoff.full$VCV[,1:25])>0.1 & autocorr.diag(mcmc.hoff.full$VCV[,1:25])<0.9999)
sum(autocorr.diag(mcmc.ellr.full$VCV[,1:25])>0.1 & autocorr.diag(mcmc.ellr.full$VCV[,1:25])<0.9999)
