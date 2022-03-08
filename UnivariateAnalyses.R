##################################################################################################################
#G matrix stability of clinally diverging populations of an annual weed
#Georgia A. Henry and John R. Stinchcombe
#Univariate analyses modified from Puentes et al., Evolution, 2016
##################################################################################################################

# Packages and functions to load first ####
library(MCMCglmm)
library(nlme)
library(ggplot2)
library(bestNormalize)
meanstd<- function (col){
  ave<-mean(col)
  stdval<-vector(, length = length(col))
  for (i in 1:length(col)){
    stdval[i]<-(col[i]-ave)/ave
  }
  stdval
}

# Index function using the mean std Cline max coeffs
Index.func<-function (SM, GR, FT, CW, ASD){
  (SM*0.252 + GR*-0.5066 + FT*-0.3454 + CW*0.6271 + ASD*-0.3997)
}
# Change GH.area to Mat.area for SM_g
# Wrapper function for mcmc
uni_table_mcmc <- function (dat, resp, resp_means, nitt=152000,burnin=2000) {
  var.env <- numeric();	var.dam <- numeric(); var.resid <- numeric()
  sign.env <- numeric();	sign.dam <- numeric();heritability <- numeric()
  mean.var <- numeric(); se.var <- numeric()
  formula = "Pop"
  ss<-eval( parse( text = as.expression( paste( resp, "~", formula, sep="" ))))
  # parameter expansion priors
  prior <- list(R=list(V=diag(4), nu=4.002),
                G=list(G1=list(V=diag(4), nu=0.002),
                       G2=list(V=diag(4), nu=0.002)))
  
  # sire = sires nested in Population with unique ids
  # dam = dams nested in sire and Population with unique ids
  mod.full <- MCMCglmm(ss, random=~idh(Pop):Dam + idh(Pop):GH.area,
                       rcov=~idh(Pop):units, data=dat, nitt=nitt, burnin=burnin,thin=100,
                       prior=prior)
  
  # Run models for each Population to examine delta DIC for variance components  
 Pop.name <- levels( dat[[ "Pop" ]])
 for (i in 1:4) {
  
  dat1<-dat[dat$Pop==Pop.name[i], ]
  mean.var[i]<-mean(dat1[[resp_means]])
  se.var[i]<-sd(dat1[[resp_means]])/sqrt(length(dat1[[resp_means]]))
  
  ss <- eval( parse( text = as.expression( paste( resp, "~", "1", sep = "" ))))
  
  
  priorP <- list(R=list(V=diag(1), nu=0.002),
                 G=list(G1=list(V=diag(1), nu=0.002),
                        G2=list(V=diag(1), nu=0.002)))    
  priorPred <- list(R=list(V=diag(1), nu=0.002),
                    G=list(G1=list(V=diag(1), nu=0.002)))    
  
  mod.start <- MCMCglmm(ss, random=~Dam + GH.area,
                        data=dat1, nitt=nitt, burnin=burnin,thin=100,
                        prior=priorP)
  mod.red <- MCMCglmm(ss, random=~ Dam,
                      data=dat1, nitt=nitt, burnin=burnin,thin=100,
                      prior=priorPred)
  
  sign.env[i] <- mod.start$DIC - mod.red$DIC
  
  mod.red <- MCMCglmm(ss, random=~ GH.area,
                      data=dat1, nitt=nitt, burnin=burnin,thin=100,
                      prior=priorPred)
  sign.dam[i] <- mod.start$DIC - mod.red$DIC
  
}
  
  Ndam <- sapply(with(dat, tapply(Dam, Pop, unique)), length)
  Nplant <- c( xtabs( ~ Pop, dat) )
  var.dam <- apply(mod.full$VCV[,1:4], 2, mean)
  var.env <- apply(mod.full$VCV[,5:8], 2, mean)
  var.resid <- apply(mod.full$VCV[,9:12], 2, mean)
  heritability <- (var.dam) / (var.dam + var.env +  var.resid)
  
  # Get uncertainty around h2
 heritabilityPENN <- (mod.full$VCV[,1]) / 
   (mod.full$VCV[,1] +  mod.full$VCV[,5] + mod.full$VCV[,9] )
 heritabilityMARY <- (mod.full$VCV[,2]) / 
   (mod.full$VCV[,2] +  mod.full$VCV[,6] + mod.full$VCV[,10] )
 heritabilityHOFF <- (mod.full$VCV[,3]) / 
   (mod.full$VCV[,3] +  mod.full$VCV[,7] + mod.full$VCV[,11] )
 heritabilityELLR <- (mod.full$VCV[,4]) / 
   (mod.full$VCV[,4] +  mod.full$VCV[,8] + mod.full$VCV[,12] )
  herUnPENN <- HPDinterval(heritabilityPENN)
  herUnMARY <- HPDinterval(heritabilityMARY)
  herUnHOFF <- HPDinterval(heritabilityHOFF)
  herUnELLR <- HPDinterval(heritabilityELLR)
  herUn <- matrix(c(herUnPENN, herUnMARY, herUnHOFF, herUnELLR),ncol=2,byrow=TRUE)
  # Multcomp of trait means with 95% CI
  i1 <- HPDinterval(mod.full$Sol[,1])
  i2 <- HPDinterval(mod.full$Sol[,2])
  i3 <- HPDinterval(mod.full$Sol[,3])
  i4 <- HPDinterval(mod.full$Sol[,4])
  
  m1 <- mean(mod.full$Sol[,1])
  m2 <- mean(mod.full$Sol[,2])
  m3 <- mean(mod.full$Sol[,3])
  m4 <- mean(mod.full$Sol[,4])

  inter <- matrix(c(i1,i2,i3,i4), ncol=2,byrow=TRUE)
  means <- matrix(c(m1,m2,m3,m4), ncol=1,byrow=TRUE)  

  return( list("response" = resp, "Pop" = Pop.name, "mean" = round(mean.var, 2), "se" = round(se.var, 2), 
               "meanHPD" = inter, "meanest"=means,  "deltaDIC_env" = round(sign.env, 3), "Dam_Sigma" = round(var.dam, 3), 
               "deltaDIC_dam" = round(sign.dam, 3),  
               "H2" = round(heritability, 3), "HPD_H2" = herUn, "No. Dams" = Ndam, "No. Plants" = Nplant ))
}

# Wrapper function for lme
uni_table_lme <- function (dat, resp, resp_means) {
  var.env <- numeric();  var.dam <- numeric(); var.resid <- numeric()
  sign.env <- numeric();	sign.dam <- numeric();heritability <- numeric()
  mean.var <- numeric(); se.var <- numeric()
  formula = "Pop"
  ss<-eval( parse( text = as.expression( paste( resp, "~", formula, sep="" ))))
  mod.full <- lme(ss, random = ~ 0 + Pop|Dam/GH.area, 
                  data = dat, weights=varIdent(form=~1|Pop), control=lmeControl(msMaxIter = 100) )
  
  # run models for each population to test variance components with LRT
  pop.name <- levels( dat[[ "Pop" ]])
  for (i in 1:4) {
    
    dat1<-dat[dat$Pop==pop.name[i], ]
    mean.var[i]<-mean(dat1[[resp_means]])
    se.var[i]<-sd(dat1[[resp_means]])/sqrt(length(dat1[[resp_means]]))
    
    ss <- eval( parse( text = as.expression( paste( resp, "~", "1", sep = "" ))))
    mod.start <- lme(ss, random = ~ 1 |Dam/GH.area, 
                     data = dat1, control=lmeControl(msMaxIter = 50) )
    
    mod.red <- lme(ss, random = ~ 1 |GH.area, 
                   data = dat1, control=lmeControl(msMaxIter = 50) )
    
    sign.dam[i] <- (1 - pchisq(2 * (logLik(mod.start,REML=TRUE) - logLik(mod.red,REML=TRUE)), 1 )) / 2
    
    mod.red <- lme(ss, random = ~ 1 |Dam, 
                   data = dat1, control=lmeControl(msMaxIter = 50) )
    sign.env[i] <- (1 - pchisq(2 * (mod.start$logLik - mod.red$logLik), 1 )) / 2
    
  }
  

  Ndam <- sapply(with(dat, tapply(Dam, Pop, unique)), length)
  Nplant <- c( xtabs( ~ Pop, dat) )
  var.dam <- as.numeric(VarCorr(mod.full)[ 2:5 ])
  var.env <- as.numeric(VarCorr(mod.full)[ 7:10 ])
  var.resid <- c(as.numeric(VarCorr(mod.full)[ 11]) * (coef(mod.full$modelStruct$varStruct, unconstrained=FALSE, allCoef=TRUE))^2)
  #var.resid <- var.resid[c(4,3,2,1)] # the weights argument in lme() dont keep the pop level order. This fixes that.
  heritability <- (var.dam) / (var.dam + var.env + var.resid)
  
  library(multcomp) # To test differences between population trait means
  # The ghlt() function requires that the model terms are in the global environment
  # Hence the repetition of "ss" and rm()
  ss <<- eval( parse( text = as.expression( paste( resp, "~", formula, sep="" ))))
  mult.sum <- summary(glht(mod.full, linfct = mcp(Pop = "Tukey")) )
  mult.letter <- cld(mult.sum)$mcletters$Letters
  rm(list = ls(envir=globalenv())[grep("ss", ls(envir=globalenv()))], envir = globalenv())
  
  
  return( list("response" = resp, "Pop" = pop.name, "mean" = round(mean.var, 2), "se" = round(se.var, 2),
               "F_stat" = round(anova(mod.full)$'F-value'[2], 2), "P-value" = round(anova(mod.full)$'p-value'[2], 3), "Sign. diff" = mult.letter, 
               "Dam_Sigma" = round(var.dam, 3), "P_dam" = round(sign.dam, 3), "Error_sigma" = round(var.resid, 3), 
               "H2" = round(heritability, 3),  "No. Dams" = Ndam, "No. Plants" = Nplant ))
}

# Import data ####
i.multi <- read.csv(file.choose(), header = T)
i.multi$Pop <- factor(i.multi$Pop, levels = c("PENN", "MARY","HOFF", "ELLR")) # Set order of levels
i.multi$Dam<-as.factor(i.multi$Dam)
i.multi$GH.area<-as.factor(i.multi$GH.area)
i.multi$Mat.area<-as.factor(i.multi$Mat.area)

#Below has already been done for i.multifull.csv
#i.multi<-na.omit(i.multi)
#i.multi$FT.s<-orderNorm(i.multi$FT)$x.t
#i.multi$SM.s <- scale(i.multi$SM_g)
#i.multi$GR.s <- scale(i.multi$GR)
#i.multi$CW.s <- scale(i.multi$CW)
#i.multi$ASD.s <- scale(i.multi$ASD)
#mean.std<-subset(i.multi, select = c(Pop, Dam, GH.area, Mat.area))
#mean.std<-cbind(mean.std, meanstd(i.multi$SM_g), meanstd(i.multi$GR), meanstd(i.multi$DTF), meanstd(i.multi$CW), meanstd(i.multi$ASD))
#i.multi$Index <- Index.func(mean.std$SM, mean.std$GR, mean.std$DTF, mean.std$CW, mean.std$ASD)

#Best prior determined via DIC score, summary available with data
prior <- list(R=list(V=diag(4), nu=4.002),
              G=list(G1=list(V=diag(4), nu=0.002),
                     G2=list(V=diag(4), nu=0.002)))

traitnames<-c("Seed mass", "Growth rate", "Flowering date", "Corolla width", "A-S distance")

### Run MCMC models ####
# resp = response used for estimating genetic variances
# reps_means = response used to calculate the Population trait means

set.nitt = 152000
set.burnin = 2000
tab.mod1 <- uni_table_mcmc( i.multi, resp = "SM_g", resp_means = "SM_g", 
                           nitt=set.nitt,burnin=set.burnin)
tab.mod2 <- uni_table_mcmc( i.multi, resp = "GR", resp_means = "GR", 
                            nitt=set.nitt,burnin=set.burnin)
tab.mod3 <- uni_table_mcmc( i.multi, resp = "DTF", resp_means = "DTF", 
                            nitt=set.nitt,burnin=set.burnin)
tab.mod4 <- uni_table_mcmc( i.multi, resp = "CW", resp_means = "CW", 
                            nitt=set.nitt,burnin=set.burnin)
tab.mod5 <- uni_table_mcmc( i.multi, resp = "ASD", resp_means = "ASD", 
                            nitt=set.nitt,burnin=set.burnin)
tab.mod6 <- uni_table_mcmc( i.multi, resp = "Index", resp_means = "Index", 
                            nitt=set.nitt,burnin=set.burnin)

# Put together the table
i.multi.table <- rbind(
  do.call(cbind.data.frame, tab.mod1),  	
  do.call(cbind.data.frame, tab.mod2),		
  do.call(cbind.data.frame, tab.mod3),		
  do.call(cbind.data.frame, tab.mod4),		
  do.call(cbind.data.frame, tab.mod5),
  do.call(cbind.data.frame, tab.mod6))		
View(i.multi.table)# Print table

#save table
write.csv(i.multi.table, file = "i.multi_table_mcmc.csv")

### Run REML models ####
# Reproduce REML analyses presented as supplemental material
# resp = response used for estimating genetic variances
# reps_means = response used to calculate the Population trait means
tab.mod1LME <- uni_table_lme( i.multi, resp = "SM.s", resp_means = "SM.s")
tab.mod2LME <- uni_table_lme( i.multi, resp = "GR.s", resp_means = "GR.s")
tab.mod3LME <- uni_table_lme( i.multi, resp = "FT.s", resp_means = "FT.s")
tab.mod4LME <- uni_table_lme( i.multi, resp = "CW.s", resp_means = "CW.s")
tab.mod5LME <- uni_table_lme( i.multi, resp = "ASD.s", resp_means = "ASD.s")
tab.mod6LME <- uni_table_lme( i.multi, resp = "Index", resp_means = "Index")

# Put table together
i.multi.tableLME <- rbind(
  do.call(cbind.data.frame, tab.mod1LME),		
  do.call(cbind.data.frame, tab.mod2LME),		
  do.call(cbind.data.frame, tab.mod3LME),		
  do.call(cbind.data.frame, tab.mod4LME),		
  do.call(cbind.data.frame, tab.mod5LME),
  do.call(cbind.data.frame, tab.mod6LME))		
View(i.multi.tableLME)

#save table
write.csv(i.multi.tableLME, file = "i.multi_table_reml.csv")

###H2 for figure S1 ####
hert.list<-list()
hert.list[[1]]<-list(H2=matrix(tab.mod1$H2), HPD=tab.mod1[["HPD_H2"]])
hert.list[[2]]<-list(H2=matrix(tab.mod2$H2), HPD=tab.mod2[["HPD_H2"]])
hert.list[[3]]<-list(H2=matrix(tab.mod3$H2), HPD=tab.mod3[["HPD_H2"]])
hert.list[[4]]<-list(H2=matrix(tab.mod4$H2), HPD=tab.mod4[["HPD_H2"]])
hert.list[[5]]<-list(H2=matrix(tab.mod5$H2), HPD=tab.mod5[["HPD_H2"]])
hert.list[[6]]<-list(H2=matrix(tab.mod6$H2), HPD=tab.mod6[["HPD_H2"]])
names(hert.list)<-c("Seed mass", "Growth rate", "Flowering time", "Corolla width", "Anther-Stigma dist." , "Index")



# m is number of populations
# t is number of traits (5 plus index)
m=4
t=6
# HPD overlap Comparisons for each population hertiability estimate
prs <- cbind(rep(1:m, each = m), 1:m) 
prs.comp <- prs[prs[,1] < prs[,2], , drop = FALSE]
#setting up an index for hert
hertoverlap.score <-matrix(,t,((m^2 - m)/2))
for (j in 1:t){
  hertoverlap.score[j,] <- ifelse(hert.list[[j]][["HPD"]][prs.comp[,1],1] > hert.list[[j]][["HPD"]][prs.comp[,2],2] | hert.list[[j]][["HPD"]][prs.comp[,2],1] > hert.list[[j]][["HPD"]][prs.comp[,1],2],0,1) 
}
colnames(meanoverlap.score)<-c("P-M","P-H","P-E","M-H","M-E","H-H")
t=6

# Figure S1 ####

par(mfrow=c(2,3))

for (i in 1:t) {
  plot(hert.list[[i]][["H2"]], ylim = c(0,1), xaxt="n", pch=16, xlim = c(0.8, 4.2), xlab ="", ylab=names(hert.list)[[i]])
  axis(1,at=1:4,labels=Gnames)
  arrows(1:4,hert.list[[i]][["H2"]],1:4,hert.list[[i]][["HPD"]][,1],length=0.1,angle=90)
  arrows(1:4,hert.list[[i]][["H2"]],1:4,hert.list[[i]][["HPD"]][,2],length=0.1,angle=90)
}


