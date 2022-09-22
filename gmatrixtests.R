#################################
#G matrix stability of clinally diverging populations of an annual weed 
#Georgia A. Henry and John R. Stinchcombe
#
#Much of the code below is modified from: 
#Aguirre et al., Heredity, 2014 and Puentes et al., Evolution, 2016
#
#Modifications and additions made by G.A.H.
# georgia.henry@mail.utoronto.ca
################################
#### Packages to load in ####
library(matrixcalc)
library(MCMCglmm)
library(psych)
library(MasterBayes)
library(bestNormalize)
library(gdata)
#### Functions to load in ####
gmat.list <- function(x){list(matrix(x,nrow=5,ncol=5,byrow=F))} #Puentes
R.proj <- function(Gs,p,vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]]
  rand.vec <-matrix(,vec,n)
  for (i in 1:vec){
    b <- runif(n,-1,1)
    rand.vec[i,] <- b/(sqrt(sum(b^2)))
  }
  #generate unit length random vectors  
  proj<- function(G,b) t(b) %*% G %*% (b)
  #internal function to do projection
  G.proj <- array(,c(MCMCsamp, m, vec))
  colnames(G.proj) <- dimnames(Gs)[[3]]
  for (i in 1:vec){
    G.proj[,,i]<- t(apply(Gs, 3:4, proj, b = rand.vec[i,]))
  }
  #project each random vector through each MCMC sample of each G
  prs <- cbind(rep(1:m, each = m), 1:m) 
  prs.comp <- prs[prs[,1] < prs[,2], , drop = FALSE] 
  #setting up an index for HPD comparisons
  proj.score <-matrix(,vec,((m^2 - m)/2))
  for (k in 1:vec){
    HPD.int <- HPDinterval(as.mcmc(G.proj[,,k]), prob = p)
    proj.score[k,] <- ifelse(HPD.int[prs.comp[,1],1] > HPD.int[prs.comp[,2],2] | HPD.int[prs.comp[,2],1] > HPD.int[prs.comp[,1],2],1,0) 
  }
  #for a given random vector, examine if the HPD intervals of any pair of G matrices overlap
  vec.score <-cbind(rand.vec, proj.score)
  colnames(vec.score) <- c(1:n, paste(dimnames(Gs)[[3]][prs.comp[, 1]], ".vs.", dimnames(Gs)[[3]][prs.comp[, 2]], sep = ""))
  #collate the random vectors and the outcome of their projection on the G matrices
  sig.vec <- subset(vec.score, rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0) 
  #collate just the random vectors that resulted in significant differences in variance
  if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger vec or lower p"); eig.R <- "Na"}
  else{
    eig.R <- eigen(cov(sig.vec[,1:n]))
    rownames(eig.R$vectors) <- dimnames(Gs)[[1]]
    colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
  }  
  #eigen analysis of the R matrix
  list(G.proj = G.proj, vec.score = vec.score, eig.R = eig.R)
} #Aguirre
kr.subspace <- function(Gs, vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]] 
  if(length(vec) != m){stop("vec must have length = m")}
  h <- function (g, v){
    AA <- array(, c(n, n, m))  
    for (k in 1:m){
      g.vec <- eigen(g[,,k])$vectors[,1:(v[k])] 
      AA[,,k] <- g.vec %*% t(g.vec)
    }
    H <- apply(AA, 1:2, sum)
    list(H = H, AA = AA)
  }
  #internal function to calculate AA and H
  MCMC.H <- array(, c(n, n, MCMCsamp))
  dimnames(MCMC.H) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[4]])      
  MCMC.AA <- array(, c(n, n, m, MCMCsamp))
  dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]], dimnames(Gs)[[4]])
  for (i in 1:MCMCsamp){
    kr <- h(Gs[,,,i], v = vec)
    MCMC.H[,,i] <- kr$H
    MCMC.AA[,,,i] <- kr$AA
  }	
  #calculate AA and H for the ith MCMC sample of the G array		
  avH <- apply(MCMC.H, 1:2, mean)
  rownames(avH) <- dimnames(Gs)[[1]]
  colnames(avH) <- dimnames(Gs)[[1]]
  #calculate the posterior mean H
  avAA <- apply(MCMC.AA, 1:3, mean)
  dimnames(avAA) <- list(dimnames(Gs)[[1]], dimnames(Gs)[[1]], dimnames(Gs)[[3]])
  #calculate the posterior mean AA
  avH.vec <- eigen(avH)$vectors
  #eigenanalysis of posterior mean H	
  proj<- function(a, b) t(b) %*% a %*% b
  #internal function to do projection
  avH.theta <- matrix(, n, m)
  for (i in 1:n){
    for (i in 1:n){
      avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #angles between the eigenvectors posterior mean H and the posterior mean subspaces of each population
  MCMC.H.val <- matrix(, MCMCsamp, n)
  colnames(MCMC.H.val) <- paste("h", 1:n, sep="")
  for (i in 1:n){
    MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean H 
  MCMC.H.theta <- array(, c(n, m, MCMCsamp))
  rownames(MCMC.H.theta) <- paste("h", 1:n, sep="")
  colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
  for(i in 1:n){
    for(j in 1:MCMCsamp){
      MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  #posterior distribution of the angles between the eigenvectors of posterior mean H and the MCMC samples of the subspaces of each population
  list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
} #Aguirre
covtensor <- function(Gs){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  neigten <- n*(n+1)/2 
  #Number of eigentensors
  MCMC.S <- array(,c(neigten, neigten, MCMCsamp))
  dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
  for (k in 1:MCMCsamp){
    MCMCG <- Gs[,,,k] 
    MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
    #find the variances of the kth G and store them 
    MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
    #find the covariances of the kth G and store them
    MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
    #fill the upper left quadrant of the kth S
    MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
    #fill the lower right quadrant of the kth S
    MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
    #fill the upper right quadrant of the kth S
    MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
    #fill the lower left quadrant of the kthS
  }  
  av.S <- apply(MCMC.S, 1:2, mean)
  #posterior mean S
  av.S.val <- eigen(av.S)$values
  #eigenvalues of posterior mean S 
  av.S.vec <- eigen(av.S)$vectors
  #eigenvalues of posterior mean S
  eTmat <- array(, c(n, n, neigten))
  dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
  for (i in 1:neigten){
    emat <- matrix(0, n, n) 
    lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
    emat <- emat + t(emat)
    diag(emat) <- av.S.vec[1:n,i]
    eTmat[,,i] <- emat 
  }
  #construct the second-order eigentensors of posterior mean S
  eT.eigen <- array(, c(n+1, n, neigten))
  for (i in 1:neigten){
    eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
    #Eigenvalues of the ith eigentensor
    eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
    #Eigenvectors of the ith eigentensor
    eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = T), i]
  }
  MCMC.S.val <- matrix(, MCMCsamp, neigten)
  colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
  for (i in 1:MCMCsamp){
    for(j in 1:neigten){
      MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
    }
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
  av.G.coord <- array(, c(m, neigten, 1))
  dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the jth avG for the eigentensors of posterior mean S
  MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
  dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
  tensor.summary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
  colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
  rownames(tensor.summary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
} #Aguirre
proj<- function(G, b){ t(b) %*% G %*% (b)}
vec.corr<-function(z1=z1,z2=z2){
  (sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}
cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X,use="complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R
}
eff.dim<-function(G){
  sum(eigen(G)$values)/eigen(G)$values[1]
}

#### Observed G matrices ####
#Change file path to correct one
mcmcDiagFull <- readRDS("mcmc.diag.full.RDS")

### Posterior samples of G
Gmcmc.penn<- lapply(apply(mcmcDiagFull$VCV[,1:25]/100, 1, gmat.list), "[[", 1) 
Gmcmc.mary<- lapply(apply(mcmcDiagFull$VCV[,26:50]/100, 1, gmat.list), "[[", 1) 
Gmcmc.hoff<- lapply(apply(mcmcDiagFull$VCV[,51:75]/100, 1, gmat.list), "[[", 1) 
Gmcmc.ellr<- lapply(apply(mcmcDiagFull$VCV[,76:100]/100, 1, gmat.list), "[[", 1) 

### Mean and median of posterior samples
Gmcmc.mean.penn<- matrix(colMeans(mcmcDiagFull$VCV[,1:25]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.mary<- matrix(colMeans(mcmcDiagFull$VCV[,26:50]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.hoff<- matrix(colMeans(mcmcDiagFull$VCV[,51:75]),nrow=5,ncol=5,byrow=F)/100
Gmcmc.mean.ellr<- matrix(colMeans(mcmcDiagFull$VCV[,76:100]),nrow=5,ncol=5,byrow=F)/100

mean.Gs<-list(Gmcmc.mean.penn, Gmcmc.mean.mary, Gmcmc.mean.hoff, Gmcmc.mean.ellr)
write.csv(mean.Gs, file = "MeanGtable.csv")  

#Gmcmc.median.penn <- matrix(apply(mcmcDiagFull$VCV[,1:25]),2,median),nrow=5,ncol=5,byrow=F)/100
#Gmcmc.median.mary <- matrix(apply(mcmcDiagFull$VCV[,26:50]),2,median),nrow=5,ncol=5,byrow=F)/100
#Gmcmc.median.hoff <- matrix(apply(mcmcDiagFull$VCV[,51:75]),2,median),nrow=5,ncol=5,byrow=F)/100
#Gmcmc.median.ellr <- matrix(apply(mcmcDiagFull$VCV[,76:100]),2,median),nrow=5,ncol=5,byrow=F)/100


#### Random G matrices ####
#loading in rand.Garray produced in randomnull.R#
#Change file path to correct one
rand.Garray<-readRDS("randomdiagGarray.RDS")
#The *10 has already been factored out
### Mean and median of Random G Array posterior samples
rand.mean.penn<- apply(rand.Garray[,,1,], 1:2, mean)
rand.mean.mary<- apply(rand.Garray[,,2,], 1:2, mean)
rand.mean.hoff<- apply(rand.Garray[,,3,], 1:2, mean)
rand.mean.ellr<- apply(rand.Garray[,,4,], 1:2, mean)
rmean.Gs<-list(rand.mean.penn, rand.mean.mary, rand.mean.hoff, rand.mean.ellr)
write.csv(rmean.Gs, file = "randMeanGtable.csv")  

rand.median.penn<-matrix(apply(rand.Garray[,,1,], 1:2, median),nrow=5,ncol=5,byrow=F)
rand.median.mary<-matrix(apply(rand.Garray[,,2,], 1:2, median),nrow=5,ncol=5,byrow=F)
rand.median.hoff<-matrix(apply(rand.Garray[,,3,], 1:2, median),nrow=5,ncol=5,byrow=F)
rand.median.ellr<-matrix(apply(rand.Garray[,,4,], 1:2, median),nrow=5,ncol=5,byrow=F)

rand.Gmcmc.penn<- lapply(apply(rand.Garray[,,1,], 3, list), "[[", 1) 
rand.Gmcmc.mary<- lapply(apply(rand.Garray[,,2,], 3, list), "[[", 1) 
rand.Gmcmc.hoff<- lapply(apply(rand.Garray[,,3,], 3, list), "[[", 1) 
rand.Gmcmc.ellr<- lapply(apply(rand.Garray[,,4,], 3, list), "[[", 1) 

rpop.Gs <- list(rand.Gmcmc.penn, rand.Gmcmc.mary, rand.Gmcmc.hoff, rand.Gmcmc.ellr)
rGarray.pop <- array( , c(n,n,m,ll))
dimnames(rGarray.pop) <- list(traitnames,traitnames,Gnames, c(1:1000))
for (i in 1:m) {
  for (j in 1:ll) {
    rGarray.pop[,,i,j] <- rpop.Gs[[i]][[j]]
  }
}
randGmat_HPD<-function (randGmat) {
  fhpd <- function (x) { x <- as.mcmc(x); HPDinterval(x, prob=0.95)}
  hpd <- round(apply(randGmat,1:2, fhpd ),3)
  int.res <- paste(hpd[1,,],hpd[2,,],sep=" , ")
  mat.int <- matrix(int.res,nrow=5,ncol=5, byrow=F)
  return(mat.int)
}
rGHPD.penn<-randGmat_HPD(rGarray.pop[,,1,])
rGHPD.mary<-randGmat_HPD(rGarray.pop[,,2,])
rGHPD.hoff<-randGmat_HPD(rGarray.pop[,,3,])
rGHPD.ellr<-randGmat_HPD(rGarray.pop[,,4,])

rGHPD.all<-list(rGHPD.penn,rGHPD.mary,rGHPD.hoff,rGHPD.ellr)
write.csv(rGHPD.all, file = "randHPDAll.csv")


#### Assign variables and constants ####

#clinemax<-as.vector(c(0.265198245, 0.506647609, -0.345393991, 0.627613055, -0.399728915))#SD clinemax
clinemaxsd<-as.vector(c(0.4298,-0.5351,-0.3528,0.2625,-0.7384))
MCMCsamp<-length(Gmcmc.penn)
m=4 #number of populations
n=5 #number of traits
r=3 #number of random effects in model: dam, greenhouse area, and residual
traitnames <- c("Seed mass", "Growth", "Flow time", "Corolla", "AS dist.") #trait names
Gnames <- c("Penn","Mary", "Hoff","Ellr")#population names
pop.Gs <- list(Gmcmc.penn, Gmcmc.mary, Gmcmc.hoff, Gmcmc.ellr)

#Construct an array with posterior Gs
Garray.pop <- array( , c(n,n,m,MCMCsamp))
dimnames(Garray.pop) <- list(traitnames,traitnames,Gnames)
for (i in 1:m) {
  for (j in 1:MCMCsamp) {
    Garray.pop[,,i,j] <- pop.Gs[[i]][[j]]
  }
}

#### G matrix dimensionality ####
### mean gmaxes
mean.g.list<-array(c(Gmcmc.mean.penn, Gmcmc.mean.mary, Gmcmc.mean.hoff, Gmcmc.mean.ellr), c(5,5,4))
dimnames(mean.g.list)<-list(traitnames, traitnames, pop.list)
eigenGs<-apply(mean.g.list, 3, FUN=eigen)
mean.gmaxes<-array(,c(5,4))
dimnames(mean.gmaxes)<-list(traitnames, pop.list)
for (i in 1:m){
  mean.gmaxes[,i]<-eigenGs[[i]][[2]][,1]
}
write.csv(mean.gmaxes, file = "Mean.gmaxes.csv")
### Trace ####
trace.array<-array(,c(4,MCMCsamp))
hpd.trace<-array(,c(m,2))
for (k in 1:m) {
  for (i in 1:MCMCsamp){
  trace.array[k,i]<-tr(Garray.pop[,,k,i])
  }
  hpd.trace[k,]<-HPDinterval(as.mcmc(trace.array[k,]))
  }

trace.matrix<-array(,c(4,3))
trace.matrix[,1]<-apply(trace.array, 1, mean)
for (i in 1:m){
 trace.matrix[i,2]<-hpd.trace[i,1]
 trace.matrix[i,3]<-hpd.trace[i,2]
}
colnames(trace.matrix)<-c("Trace","low","high")
rownames(trace.matrix)<-c("Penn", "Mary","Hoff","Ellr")
View(trace.matrix)
write.csv(trace.matrix, file="Gmatrixtrace95.csv")

rtrace.array<-array(,c(4,1000))
rhpd.trace<-array(,c(m,2))
for (k in 1:m) {
  for (i in 1:1000){
    rtrace.array[k,i]<-tr(rand.Garray[,,k,i])
  }
  rhpd.trace[k,]<-HPDinterval(as.mcmc(rtrace.array[k,]))
}

rtrace.matrix<-array(,c(4,3))
rtrace.matrix[,1]<-apply(rtrace.array, 1, mean)
for (i in 1:m){
  rtrace.matrix[i,2]<-rhpd.trace[i,1]
  rtrace.matrix[i,3]<-rhpd.trace[i,2]
}
colnames(rtrace.matrix)<-c("Trace","low","high")
rownames(rtrace.matrix)<-c("Penn", "Mary","Hoff","Ellr")
View(rtrace.matrix)
write.csv(rtrace.matrix, file="Gmatrixtracerand95.csv")

combs<-list(c(1,3), c(1,4), c(2,3), c(2,4))
tracediff<-array(,c(4,MCMCsamp))
for (i in 1:length(combs)){
  cm<-combs[[i]]
  tracediff[i,]<-trace.array[cm[1],]/trace.array[cm[2],]
}
dimnames(tracediff)<-list(c("Penn v. Hoff", "Penn v. Ellr", "Mary v. Hoff", "Mary v. Ellr"),NULL)

#ROPE of +/- 0.10 --> within ROPE the traces are probably the same
hpd.trdif<-array(,c(4,2))
for (i in 1:length(combs)){
  hpd.trdif[i,]<- HPDinterval(as.mcmc(tracediff[i,]))
}
rownames(hpd.trdif)<-list("Penn v. Hoff", "Penn v. Ellr", "Mary v. Hoff", "Mary v. Ellr")
for (i in 1:length(combs)){
  hist(tracediff[i,], col = "#2B6364", xlim = c(0,1.5), main=dimnames(tracediff)[[1]][i])
  abline(v=0.9, lwd = 2, lty="dashed", col = "#C66234")
  abline(v=1.1, lwd = 2, lty="dashed", col = "#C66234")
  abline(v=c(hpd.trdif[i,1],hpd.trdif[i,2]), col = "#519ea0", lwd = 2)
}  
  
### Genetic variation in the direction of maximum clinal divergence  ####
proj.array<-array(,c(m,MCMCsamp))
hpd.proj<-array(,c(m,2))
for (k in 1:m) {
  for (i in 1:MCMCsamp){
    proj.array[k,i]<-proj(Garray.pop[,,k,i], clinemaxsd)
  }
  hpd.proj[k,]<-HPDinterval(as.mcmc(proj.array[k,]))
}

proj.matrix<-array(,c(m,3))
proj.matrix[,1]<-apply(proj.array, 1, mean)
for (i in 1:m){
  proj.matrix[i,2]<-hpd.proj[i,1]
  proj.matrix[i,3]<-hpd.proj[i,2]
}
colnames(proj.matrix)<-c("Proj","low","high")
rownames(proj.matrix)<-c("Penn", "Mary","Hoff","Ellr")
View(proj.matrix)
write.csv(proj.matrix, file="Gmatrixproj95diag.csv")
#dividing by the amount of variation in Gmax standardizes the value
projstd.array<-array(,c(m,MCMCsamp))
hpdstd.proj<-array(,c(m,2))
for (k in 1:m) {
  for (i in 1:MCMCsamp){
    projstd.array[k,i]<-proj(Garray.pop[,,k,i], clinemaxsd)/eigen(Garray.pop[,,k,i])$values[1]#change
  }
  hpdstd.proj[k,]<-HPDinterval(as.mcmc(projstd.array[k,]))
}

projstd.matrix<-array(,c(m,3))
projstd.matrix[,1]<-apply(projstd.array, 1, mean)
for (i in 1:m){
  projstd.matrix[i,2]<-hpdstd.proj[i,1]
  projstd.matrix[i,3]<-hpdstd.proj[i,2]
}
colnames(projstd.matrix)<-c("Proj(std)","low","high")
rownames(projstd.matrix)<-c("Penn", "Mary","Hoff","Ellr")
View(projstd.matrix)
write.csv(projstd.matrix, file="Gmatrixproj95stddiag.csv")

p=1
theta.array<-array(,c(m,MCMCsamp))
hpd.theta<-array(,c(m,2))
for (k in 1:m) {
  for (i in 1:MCMCsamp){
    theta.array[k,i]<-theta(eigen(Garray.pop[,,k,i])$vectors[,p], clinemaxsd)
  }
  hpd.theta[k,]<-HPDinterval(as.mcmc(theta.array[k,]))
}

#### Correlation between gmax and clinemax ####
# random vectors
vec = MCMCsamp
rand.vec <-matrix(,n,vec)
for (i in 1:vec){
  b <- runif(n,-1,1)
  rand.vec[,i] <- b/(sqrt(sum(b^2)))
}
#For gmax p=1, can change to get the other PC comparisons
rMCMC.cline.corr <- array(, c(n, m, MCMCsamp))
rcline.hpd<-array(,c(n,2,m))
for (p in 1:n){
  for (k in 1:m){
    for(j in 1:MCMCsamp){
      rMCMC.cline.corr[p,k,j] <- vec.corr(eigen(Garray.pop[,,k,j])$vectors[,p], rand.vec[,j])
    }
    }
}

for (i in 1:n){
  for (j in 1:m){
    rcline.hpd[i,,j]<-HPDinterval(as.mcmc(rMCMC.cline.corr[i,j,]))
  }
}
rmean.cline.corr<-t(apply(rMCMC.cline.corr, 1:2, mean))

MCMC.cline.corr <- array(, c(n, m, MCMCsamp))
cline.hpd<-array(,c(n,2,m))
for (p in 1:n){
  for (k in 1:m){
    for(j in 1:MCMCsamp){
      MCMC.cline.corr[p,k,j] <- vec.corr(eigen(Garray.pop[,,k,j])$vectors[,p], clinemaxsd)
    }
  }
}

for (i in 1:n){
  for (j in 1:m){
    cline.hpd[i,,j]<-HPDinterval(as.mcmc(MCMC.cline.corr[i,j,]))
  }
}
mean.cline.corr<-t(apply(MCMC.cline.corr, 1:2, mean))

#Figure 1
#win.metafile("Figure_1.wmf")
pdf("Figure_1.pdf")
#par(mfrow=c(2,2),mar=c(4,4,4,4), oma = c(0,0.5,0,0))
n=1
for(p in 1:n){
  plot(1:m-0.1, rmean.cline.corr[,p],type="p",xlab="",ylab="",pch=20,cex=1.5, lwd = 2, col="darkgrey",xaxt="n",frame.plot=F, ylim = c(-1,1.3), xlim = c(0.8, 4.2), font = 2)
  mtext(side=2, line = 3, font = 2, "Correlation")
  points(1:m+0.1, mean.cline.corr[,p],type="p",xlab="",pch=20,cex=1.5, col="#2B6364")
  axis(1,at=1:4,labels=Gnames, font = 2)
  #axis(2, at=c(-1,-0.5,0,0.5,1),font = 2)
  arrows(1:m+0.1,  mean.cline.corr[,p],1:m+0.1,  cline.hpd[p,1,],length=0.1,angle=90,lwd = 2, col="#2B6364")
  arrows(1:m-0.1, rmean.cline.corr[,p],1:m-0.1, rcline.hpd[p,1,],length=0.1,angle=90,lwd = 2, col="darkgrey")
  arrows(1:m+0.1,  mean.cline.corr[,p],1:m+0.1,  cline.hpd[p,2,],length=0.1,angle=90,lwd = 2, col="#2B6364")
  arrows(1:m-0.1, rmean.cline.corr[,p],1:m-0.1, rcline.hpd[p,2,],length=0.1,angle=90,lwd = 2, col="darkgrey")
  mtext(paste("PC",p, sep = ""), side = 3, at = 1)
  #if(p==3){legend(3.3,1.5, legend = c("Clinemax", "Random"), lty = 1, col = c("#2B6364", "darkgrey"), pch = 20, cex = 1, lwd = 2, bty = "n")}
}
legend(3,1.5, legend = c("Clinemax", "Random"), lty = 1, col = c("#2B6364", "darkgrey"), pch = 20, cex = 1.2, lwd = 2, bty = "n")
dev.off()
n=5
###Random Skewers ####
set.seed(555)
MCMC.R.proj <- R.proj(Garray.pop, p = 0.95, vec = 1000)
MCMC.R.proj$vec.score[1:n,(n+1):(n+((m^2 - m)/2))]
table(rowSums(MCMC.R.proj$vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0 )
Rmat<-lapply(MCMC.R.proj$eig.R, round, digits = 3)
write.csv(Rmat, file = "Rmatrix.csv")
R.vec.proj <- array(, c(MCMCsamp, m, n))
for (i in 1:n){
  R.vec.proj[,,i] <- t(apply(Garray.pop, 3:4, proj, b = MCMC.R.proj$eig.R$vectors[,i]))
}
#Genetic variance in each population in the direction of the eigenvectors of R

HPD.R.vec.proj <- array(, c(m, 2, n))
for (i in 1:n){
  HPD.R.vec.proj[,,i] <- HPDinterval(as.mcmc(R.vec.proj[,,i]), prob = 0.95)    
}

rse1<-MCMC.R.proj[["eig.R"]]

# Figure 2
win.metafile("Figure_2.wmf")

par(mfrow=c(2,3), cex=1.2, lwd = 2, mar = c(2,4,2,1))
n=5
for (i in 1:n){
  plot(1:m,colMeans(R.vec.proj[,,i]),ylab="Vg",xlab="",pch=20, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),xaxt="n",frame.plot=F,xlim=c(0,4),ylim=c(0,1.5))
  axis(1,at=1:m,labels=Gnames)
  arrows(1:m,colMeans(R.vec.proj[,,i]),1:m,HPD.R.vec.proj[,1,i],length=0.1,angle=90, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"), lwd = 2)
  arrows(1:m,colMeans(R.vec.proj[,,i]),1:m,HPD.R.vec.proj[,2,i], length=0.1,angle=90, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"), lwd = 2)
  mtext(paste("RS",colnames(Rmat$vectors)[i],sep=""),side=3,at=0.3,font=1)
}
plot.new()
legend("left", legend = Gnames, cex = 1.2, pch = 16, lty = 1, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),bty = "n" )
dev.off()

###Krzanowski's subspace analysis####
### with % var >90 ####

val <- matrix(, n, m)
for (i in 1:m){
  avG <- apply(Garray.pop, 1:3, mean)
  val[,i] <- round(cumsum(t(eigen(avG[,,i])$values))/sum(eigen(avG[,,i])$values)*100)
}

n.vec <- rep(min(apply(ifelse(round(val,1) < 90, 1, 0), 2, sum)+1), m)
MCMCsamp = 1000
MCMCG.kr.rand90 <- kr.subspace(rand.Garray, vec = n.vec)

MCMCG.kr90 <- kr.subspace(Garray.pop, vec = n.vec)

round(eigen(MCMCG.kr90$avH)$vectors, 2)
H.eigen<-eigen(MCMCG.kr90$avH)
rand.H.eigen<-eigen(MCMCG.kr.rand90$avH)
(round(apply(MCMCG.kr90$MCMC.H.theta, 1:2, mean), 1))

### angle between each population and S
penn.H.theta90<-MCMCG.kr90$MCMC.H.theta[,1,]
mary.H.theta90<-MCMCG.kr90$MCMC.H.theta[,2,]
hoff.H.theta90<-MCMCG.kr90$MCMC.H.theta[,3,]
ellr.H.theta90<-MCMCG.kr90$MCMC.H.theta[,4,]

penn.H90<-array(,c(MCMCsamp,n))
for (i in 1:5) {
  int.penn<-penn.H.theta90[i,]
  penn.H90[,i]<-int.penn
}
mary.H90<-array(,c(MCMCsamp,n))
for (i in 1:5) {
  int.mary<-mary.H.theta90[i,]
  mary.H90[,i]<-int.mary
}
hoff.H90<-array(,c(MCMCsamp,n))
for (i in 1:5) {
  int.hoff<-hoff.H.theta90[i,]
  hoff.H90[,i]<-int.hoff
}
ellr.H90<-array(,c(MCMCsamp,n))
for (i in 1:5) {
  int.ellr<-ellr.H.theta90[i,]
  ellr.H90[,i]<-int.ellr
}
n=4
HPD.H.penn90<-HPDinterval(as.mcmc(penn.H90), 0.95)[1:n,]
HPD.H.mary90<-HPDinterval(as.mcmc(mary.H90), 0.95)[1:n,]
HPD.H.hoff90<-HPDinterval(as.mcmc(hoff.H90), 0.95)[1:n,]
HPD.H.ellr90<-HPDinterval(as.mcmc(ellr.H90), 0.95)[1:n,]

###Figures for krza90####

HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr90$MCMC.H.val))[1:n,], HPDinterval(as.mcmc(MCMCG.kr.rand90$MCMC.H.val))[1:n,])
HPD.H.theta<- cbind(HPD.H.penn90,HPD.H.mary90,HPD.H.hoff90,HPD.H.ellr90)

H.table<-cbind(colMeans(MCMCG.kr90$MCMC.H.val)[1:n],HPDinterval(as.mcmc(MCMCG.kr90$MCMC.H.val))[1:n,] )
write.csv(H.table, file = "H_table.csv")

#Figure 3

#win.metafile("Figure_3.wmf")
pdf("Figure_3.pdf")
par(mfrow=c(1,1))
plot((1:n)-0.1,colMeans(MCMCG.kr90$MCMC.H.val)[1:n], col = c("#2B6364","#2B6364", "#2B6364"), type="p",xlab="",ylab="",pch=20,cex=1.5,xaxt="n", yaxt = "n",frame.plot=F,ylim=c(0,m+0.7),xlim=c(0.2,n+0.5))
mtext(side= 2, line = 0,"Lambda", font = 2)
axis(1, at=1:n,labels=c(paste("h",rep(1:n),sep="")), lwd = 2, font = 2)
axis(2, line = -4, at=0:4, labels=0:4, lwd = 2, font = 2)
points((1:n)+0.1,colMeans(MCMCG.kr.rand90$MCMC.H.val)[1:n],type="p",pch=20,cex=1.5, col = c("black","black","black"))
arrows((1:n)-0.1,colMeans(MCMCG.kr90$MCMC.H.val)[1:n], col = c("#2B6364","#2B6364", "#2B6364"),(1:n)-0.1,HPD.H.val[,1],length=0.1,angle=90, lwd =2)
arrows((1:n)-0.1,colMeans(MCMCG.kr90$MCMC.H.val)[1:n],col = c("#2B6364","#2B6364", "#2B6364"),(1:n)-0.1,HPD.H.val[,2],length=0.1,angle=90, lwd =2)
arrows((1:n)+0.1,colMeans(MCMCG.kr.rand90$MCMC.H.val)[1:n],col = c("black","black","black"),(1:n)+0.1,HPD.H.val[,3],length=0.1,angle=90,lwd =2)
arrows((1:n)+0.1,colMeans(MCMCG.kr.rand90$MCMC.H.val)[1:n],col = c("black","black","black"),(1:n)+0.1,HPD.H.val[,4],length=0.1,angle=90,lwd =2)
legend(3.2, 1,legend=c("Observed","Randomised"),lty=c(1,1), cex = 1.3, pch = 20, bty="n", bg="white", col = c("#2B6364", "black"))
dev.off()


#Supplementary Figure 3
mean.thetas<-apply(MCMCG.kr90$MCMC.H.theta, 1:2, mean)
mean.thetas<-(mean.thetas)[1:n,]
mean.thetas<-t(mean.thetas)
par(cex = 2, mar = c(2,4,1,1))
plot  ((1:n)-0.15,  mean.thetas[1,],type="p",xlab="",ylab="Theta",pch=16,cex=1, col="#311B36",xaxt="n", yaxp = c(0, 90, 3), frame.plot=F, ylim = c(0,100), xlim = c(0.8, 4.3), lwd = 2)
points((1:n)-0.05, mean.thetas[2,],type="p",xlab="",,pch=16,cex=1, col=c("#7A6FAE"))
points((1:n)+0.05, mean.thetas[3,],type="p",xlab="",pch=16,cex=1, col="#C66234" )
points((1:n)+0.15,  mean.thetas[4,],type="p",xlab="",pch=16,cex=1, col="#8F4024"  )
axis(1,at=1:n,labels=c(paste("h",rep(1:n),sep="")), lwd=2)

arrows((1:n)-0.15,mean.thetas[1,],(1:n)-0.15, HPD.H.penn90[,1],length=0.05,angle=90, col="#311B36", lwd = 2)
arrows((1:n)-0.05,mean.thetas[2,],(1:n)-0.05, HPD.H.mary90[,1],length=0.05,angle=90, col="#7A6FAE", lwd = 2)
arrows((1:n)+0.05,mean.thetas[3,],(1:n)+0.05, HPD.H.hoff90[,1],length=0.05,angle=90, col="#C66234", lwd = 2)
arrows((1:n)+0.15,mean.thetas[4,],(1:n)+0.15, HPD.H.ellr90[,1],length=0.05,angle=90, col="#8F4024", lwd = 2)
arrows((1:n)-0.15,mean.thetas[1,],(1:n)-0.15, HPD.H.penn90[,2],length=0.05,angle=90, col="#311B36", lwd = 2)
arrows((1:n)-0.05,mean.thetas[2,],(1:n)-0.05, HPD.H.mary90[,2],length=0.05,angle=90, col="#7A6FAE", lwd = 2)
arrows((1:n)+0.05,mean.thetas[3,],(1:n)+0.05, HPD.H.hoff90[,2],length=0.05,angle=90, col="#C66234", lwd = 2)
arrows((1:n)+0.15,mean.thetas[4,],(1:n)+0.15, HPD.H.ellr90[,2],length=0.05,angle=90, col="#8F4024", lwd = 2)
legend(0.5, 90, xjust = -1, legend = c("Penn", "Mary", "Hoff", "Ellr"), col=c("#311B36","#7A6FAE","#C66234","#8F4024"), bty = "n", pch = 16, lty = 1)


### Covariance Tensor ####
n=5
MCMCsamp=1000
MCMC.covtensor <- covtensor(Garray.popd)
nnonzero <- min(n*(n+1)/2,m-1)

#randomized
MCMCsamp=1000
MCMC.covtensor.rand <- covtensor(rand.Garray)

##HPD eval
HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero]), prob=0.95))
Eigentensor<-round(cbind(unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]), HPD.eT.val[,1:2], unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),HPD.eT.val[,3:4]), 5)
write.csv(Eigentensor, file = "Eigentensors.csv")

tens.summary<-(round(MCMC.covtensor$tensor.summary[1:(n*3),2:dim(MCMC.covtensor$tensor.summary)[2]], 4))
write.csv(tens.summary, file = "covtensor_summary.csv")
#e11 vector
e11 <- c(as.numeric(MCMC.covtensor$tensor.summary[1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
#e21 vector
e21 <- c(as.numeric(MCMC.covtensor$tensor.summary[(n+1),3:dim(MCMC.covtensor$tensor.summary)[2]]))
#e12 vector
e12 <- c(as.numeric(MCMC.covtensor$tensor.summary[2,3:dim(MCMC.covtensor$tensor.summary)[2]]))
#e22 vector
e22 <- c(as.numeric(MCMC.covtensor$tensor.summary[(n+2),3:dim(MCMC.covtensor$tensor.summary)[2]]))

#e11 vector rand
r.e11 <- c(as.numeric(MCMC.covtensor.rand$tensor.summary[1,3:dim(MCMC.covtensor.rand$tensor.summary)[2]]))
#e21 vector rand
r.e21 <- c(as.numeric(MCMC.covtensor.rand$tensor.summary[(n+1),3:dim(MCMC.covtensor.rand$tensor.summary)[2]]))
#e31 vector rand
r.e31 <- c(as.numeric(MCMC.covtensor.rand$tensor.summary[((n*2)+1),3:dim(MCMC.covtensor.rand$tensor.summary)[2]]))

#genetic variance along e1, 1 for each MCMC sample of each replicate line
e11.proj <- apply(Garray.pop, 3:4, proj, b = e11)
HPD.e11 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.95)
e12.proj <- apply(Garray.pop, 3:4, proj, b = e12)
HPD.e12 <- HPDinterval(t(as.mcmc(e12.proj)),prob = 0.95)
#genetic variance along e1, 1 for each MCMC sample of each replicate line
e21.proj <- apply(Garray.pop, 3:4, proj, b = e21)
HPD.e21 <- HPDinterval(t(as.mcmc(e21.proj)),prob = 0.95)
e22.proj <- apply(Garray.pop, 3:4, proj, b = e22)
HPD.e22 <- HPDinterval(t(as.mcmc(e22.proj)),prob = 0.95)

#Mean variances along e1 and e2
av.e11.proj<- apply(e11.proj, 1, mean)
av.e21.proj<- apply(e21.proj, 1, mean)

av.e12.proj<- apply(e12.proj, 1, mean)
av.e22.proj<- apply(e22.proj, 1, mean)

#Correlation between e11 and the first eigenvector from Random Skewers
round(vec.corr(e11, rse1$vectors[,1]),3)

###figures for tensor####
HPD.tensor.coord <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord) <- list(Gnames,c("lower","upper"), paste("E",1:3,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.95)[1:2]
  }
}


n=4
#Figure 4
tiff("Figure_4.tiff")
layout.matrix <- matrix(c(0,1,1,1,0, 2,0, 3), nrow = 4, ncol = 2)
layout(mat = layout.matrix,
       heights = c(0.2, 1, 0.2, 1), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns
layout.show(3)
par(mar = c(3,3,0,0), mai = c(0.5,1,0.5,0))
plot((1:nnonzero)-0.1,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),  xlab="",ylab="Alpha",pch=16,xaxt="n",  frame.plot=F, xlim=c(0.8,3.5),ylim=c(0,0.25),col = ("#2B6364"))
axis(1,at=1:nnonzero,labels=c(paste("E",rep(1:nnonzero),sep="")))
points((1:nnonzero)+0.1, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),pch=16,cex=1)
arrows((1:nnonzero)-0.1, col = ("#2B6364"),  unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)-0.1,HPD.eT.val[,1],length=0.08,angle=90)
arrows((1:nnonzero)-0.1, col = ("#2B6364"), unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)-0.1,HPD.eT.val[,2],length=0.08,angle=90)
arrows((1:nnonzero)+0.1, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)+0.1,HPD.eT.val[,3],length=0.08,angle=90,lty=5)
arrows((1:nnonzero)+0.1, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),(1:nnonzero)+0.1,HPD.eT.val[,4],length=0.08,angle=90,lty=5)
legend(2, 0.23, cex = 1.2, legend=c("Observed","Randomised"),lty=c(1,5),pch=c(16,16),bty="n", col = c("#2B6364", "black"))
mtext(LETTERS[1], side =3, line = 0, at = 1, font = 2)
par(mar = c(3,3,0,0), mai = c(0.5,0.5,0.2,0.2))
for (k in 1:2){  
  plot(1:m,MCMC.covtensor$av.G.coord[,k,],col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024") , ylab="Vg",xlab="",pch=16,xaxt="n",frame.plot=F,xlim=c(0.9,4.5),ylim=c(floor(min(HPD.tensor.coord[,,k])),ceiling(max(HPD.tensor.coord[,,k]))),main = "")
  axis(1,at=1:m,labels=Gnames)
  arrows(1:m,MCMC.covtensor$av.G.coord[,k,],1:m,HPD.tensor.coord[,1,k],length=0.1,angle=90, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024") )
  arrows(1:m,MCMC.covtensor$av.G.coord[,k,],1:m,HPD.tensor.coord[,2,k],length=0.1,angle=90, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024") )
  mtext(dimnames(MCMC.covtensor$av.G.coord)[[2]][k],side=3, line = 0,at=3.4,font=1)
  mtext(LETTERS[k+1], side = 3, line = 0, at = 1, font = 2)
}
dev.off()
##Figure 5
#win.metafile("Figure_5.wmf")
pdf("Figure_5.pdf")
par(mfrow=c(2,2))
plot(1:m,rowMeans(e11.proj),ylab="Vg",xlab="",pch=16, cex = 1, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),xaxt="n",frame.plot=F,xlim=c(0.75,4),ylim=c(0,ceiling(max(HPD.e11))))
axis(1,at=1:m,labels=Gnames)
arrows(1:m,rowMeans(e11.proj),1:m,HPD.e11[,1],length=0.05,angle=90,col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"))
arrows(1:m,rowMeans(e11.proj),1:m,HPD.e11[,2],length=0.05,angle=90,col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"))
mtext("e11",side=3,at=1.1,font=2)
plot(1:m,rowMeans(e12.proj),ylab="Vg",xlab="",pch=16, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),xaxt="n",frame.plot=F,xlim=c(0.75,4),ylim=c(0,ceiling(max(HPD.e12))))
axis(1,at=1:m,labels=Gnames)
arrows(1:m,rowMeans(e12.proj),1:m,HPD.e12[,1],length=0.05,angle=90, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"))
arrows(1:m,rowMeans(e12.proj),1:m,HPD.e12[,2],length=0.05,angle=90, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"))
mtext("e12",side=3,at=1.1,font=2)
legend("topright", legend = Gnames, pch = 16, lty = 1, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),bty = "n" )

plot(1:m,rowMeans(e21.proj),ylab="Vg",xlab="",pch=16,cex = 1,col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),xaxt="n",frame.plot=F,xlim=c(0.75,4),ylim=c(0,ceiling(max(HPD.e21))))
axis(1,at=1:m,labels=Gnames)
arrows(1:m,rowMeans(e21.proj),1:m,HPD.e21[,1],length=0.05, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),angle=90)
arrows(1:m,rowMeans(e21.proj),1:m,HPD.e21[,2],length=0.05, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),angle=90)
mtext("e21",side=3,at=1.1,font=2)
plot(1:m,rowMeans(e22.proj),ylab="Vg",xlab="",pch=16,cex = 1,col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),xaxt="n",frame.plot=F,xlim=c(0.75,4),ylim=c(0,ceiling(max(HPD.e22))))
axis(1,at=1:m,labels=Gnames)
arrows(1:m,rowMeans(e22.proj),1:m,HPD.e22[,1],length=0.05, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),angle=90)
arrows(1:m,rowMeans(e22.proj),1:m,HPD.e22[,2],length=0.05, col=c("#311B36", "#7A6FAE", "#C66234", "#8F4024"),angle=90)
mtext("e22",side=3,at=1.1,font=2)
dev.off()

