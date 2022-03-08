library(MCMCglmm)
library(mcmcplots)
library(forcats)
library(dplyr)
library(plyr)
####Subset populations####
#Change file path of choose file path for i.multifull.csv
#i.multi <- read.csv(file.choose(), header = T)
i.multi.all <- as.data.frame(subset(i.multi, select= c("Dam", "Pop", "GH.area", "SM.s", "GR.s", "FT.s", "CW.s", "ASD.s")))
names(i.multi.all)[2]<-"randDam"
Dams<-unique(subset(i.multi.all, select = Dam))

Penn_levels = as.character(c(101:150))
Mary_levels = as.character(c(201:250))
Hoff_levels = as.character(c(301:350))
Ellr_levels = as.character(c(401:434, 436:450)) #dam 435 died prior to seed set

poplev<-data.frame("Pop" = c(rep("PENN", length(Penn_levels)),rep("MARY", length(Mary_levels)),rep("HOFF", length(Hoff_levels)),rep("ELLR", length(Ellr_levels)) ))
poplev$Dam<-c(Penn_levels,Mary_levels,Hoff_levels,Ellr_levels)
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
P<-phen.cov$all
#Same prior as for the observed G matrices
priorsg<- list(R=list(V=diag(5), nu=4.001), 
               G = list(G1 = list(V=P*10/3,nu=4.001),
                        G2 = list(V=P*10/3,nu=4.001)))
nitt<-20000
burnin<-5000
thin<-100



#### Loop to run randomized within pop G models ####
#You can run for shorter stints by changing range of i
for (i in 1:1000){
  
  #shuffle individuals across dams
  i.multi.penn$randDam <- sample(i.multi.penn$Dam, length(i.multi.penn$Dam), replace = F)
  i.multi.mary$randDam <- sample(i.multi.mary$Dam, length(i.multi.mary$Dam), replace = F)
  i.multi.hoff$randDam <- sample(i.multi.hoff$Dam, length(i.multi.hoff$Dam), replace = F)
  i.multi.ellr$randDam <- sample(i.multi.ellr$Dam, length(i.multi.ellr$Dam), replace = F)
  
  ##Penn
  mcmcPenn<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                     random=~us(trait):randDam + us(trait):GH.area,family=rep("gaussian", 5), data=i.multi.penn,
                     rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
  
  saveRDS(mcmcPenn, paste("Penn",i,"poprand.RDS",sep=""))
  
  ##Mary 
  
  mcmcMary<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                     random=~us(trait):randDam + us(trait):GH.area,family=rep("gaussian",5), data=i.multi.mary,
                     rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
  
  saveRDS(mcmcMary, paste("Mary",i,"poprand.RDS",sep=""))
  
  ##Hoff
  
  mcmcHoff<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                     random=~us(trait):randDam + us(trait):GH.area,family=rep("gaussian",5), data=i.multi.hoff,
                     rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
  
  saveRDS(mcmcHoff, paste("Hoff",i,"poprand.RDS",sep=""))
  
  ##Ellr
  
  mcmcEllr<-MCMCglmm(c(SM.s,GR.s,FT.s,CW.s,ASD.s) *10 ~ trait-1,
                     random=~us(trait):randDam + us(trait):GH.area,family=rep("gaussian",5), data=i.multi.ellr,
                     rcov=~us(trait):units, nitt=nitt, prior=priorsg, verbose=T, thin=thin, burnin=burnin) 
  
  saveRDS(mcmcEllr, paste("Ellr",i,"poprand.RDS",sep=""))
  
}  



#### Loading them in and evaluating convergence ####
#Example with 5 samples from Penn. change for other chains and populations
Penn1rand <- readRDS("C:/Users/georg/Penn1poprand.RDS")
Penn2rand <- readRDS("C:/Users/georg/Penn2poprand.RDS")
Penn3rand <- readRDS("C:/Users/georg/Penn3poprand.RDS")
Penn4rand <- readRDS("C:/Users/georg/Penn4poprand.RDS")
Penn5rand <- readRDS("C:/Users/georg/Penn5poprand.RDS")
P.mcmc.list<-mcmc.list(Penn1rand$Sol, Penn2rand$Sol, Penn3rand$Sol,Penn4rand$Sol,Penn5rand$Sol)
par(mfrow = c(4,4))
rmeanplot(P.mcmc.list, plot.title = "Pennsylvania")
gelman.diag(P.mcmc.list) #should be close to 1 (1.02 in this case)


#### Loading them in and sampling from posterior ####

## import_multiple_csv_files_to_R
#make sure directory is same as where the null output from above is stored.
#setwd("C:/Users/georg")
dir()

# list all RDS files from the current directory
# use regex to import all null model files
list.files(pattern="poprand.RDS$") 

# create a list from these files
list.filenames<-list.files(pattern="poprand.RDS$")
list.filenames

# create an empty list that will serve as a container to receive the incoming files
list.data<-list()

list.data[[1]]<-readRDS(list.filenames[1])$VCV

# create a loop to read in the data
for (i in 1:length(list.filenames)){
  list.data[[i]]<-readRDS(list.filenames[i])$VCV
}

# add the names of your data to the list
names(list.data)<-list.filenames
m=4
ll<-1000
rand.VCV<-list()
rand.VCV[[1]]<-list.data[1:ll]
rand.VCV[[2]]<-list.data[(ll+1):(ll*2)]
rand.VCV[[3]]<-list.data[(ll*2+1):(ll*3)]
rand.VCV[[4]]<-list.data[(ll*3+1):(ll*4)]
names(rand.VCV)<-c("Ellr", "Hoff", "Mary", "Penn")  
#saveRDS(rand.VCV, file = "randompopVCV.RDS")
#This file is large 
n=5
m=4
rand.Garray<-array(,c(n,n,m,ll)) #this was alphabetical, usual order is the opposite so this is fixed here
for(i in 1:ll){
  j<-150
  rand.Garray[,,4,i]<-matrix(rand.VCV[[1]][[i]][j,1:25],nrow=5,ncol=5)/100
  rand.Garray[,,3,i]<-matrix(rand.VCV[[2]][[i]][j,1:25],nrow=5,ncol=5)/100
  rand.Garray[,,2,i]<-matrix(rand.VCV[[3]][[i]][j,1:25],nrow=5,ncol=5)/100
  rand.Garray[,,1,i]<-matrix(rand.VCV[[4]][[i]][j,1:25],nrow=5,ncol=5)/100
}

saveRDS(rand.Garray, file = "randompopGarray.RDS")
#this is the file for analyses

# set wd back to whatever you want if you've changed it

#setwd("C:/Users/georg")

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

#for the randomized G matrices

rGHPD.penn<-Gmat_HPD(rand.Garray[,,1,])
rGHPD.mary<-Gmat_HPD(rand.Garray[,,2,])
rGHPD.hoff<-Gmat_HPD(rand.Garray[,,3,])
rGHPD.ellr<-Gmat_HPD(rand.Garray[,,4,])

rGHPD.all<-list(rGHPD.penn,rGHPD.mary,rGHPD.hoff,rGHPD.ellr)
write.csv(rGHPD.all, file = "randHPDAll.csv")
