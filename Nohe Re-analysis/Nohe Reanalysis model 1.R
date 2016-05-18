#Nohe et al. Reanalysis
# The goal is to re-create Nohe et al.'s models on uncorrected correlations

options(scipen=99)

#load libraries
library(metaSEM)
require('matrixcalc')
require('lavaan')
require('semTools')
library(Matrix)
source('/Users/pdownes/Dropbox/Meta SEM/random.matrices.R') 

# Retrieve Data
library(gdata)
data <- read.xls('/Users/pdownes/Dropbox/Meta SEM/MASEM Shared Folder/vector files/Nohe et al vectors with reliabilities.xlsx',sheet = 1) #note there are 2 sheets due to 2 models

####################################
# Specify the Model
####################################
varnames <- c('W1','W2','S1','S2')
A <- mxMatrix('Full',ncol=4,nrow=4,byrow=T,
              values = c(0,.2,0,.2,
                         0,0,0,0,
                         0,.2,0,.2,
                         0,0,0,0),
              free=c(F,T,F,T,
                     F,F,F,F,
                     F,T,F,T,
                     F,F,F,F
                     ),
              labels=c(NA,"betaw1w2",NA,"betaw1s2",
                       NA,NA,NA,NA,
                       NA,"betas1w2",NA,"betas1s2",
                       NA,NA,NA,NA
              ),
              name="A")
S <- mxMatrix('Full',ncol=4,nrow=4,byrow=T,
              values = c(1,0,.2,0,
                         0,1,0,.2,
                         .2,0,1,0,
                         0,.2,0,1),
              free=c(F,F,T,F,
                     F,F,F,T,
                     T,F,F,F,
                     F,T,F,F),
              labels=c("varw1",NA,"covs1w1",NA,
                       NA,"varw2",NA,"covs2w2",
                       "covs1w1",NA,"vars1",NA,
                       NA,"covs2w2",NA,"vars2"
              ),
              name="S")
matrF <- mxMatrix(type="Iden",nrow=4,ncol=4,name="F")
exp         <- mxExpectationRAM("A","S","F", dimnames=varnames )

####################################
# conduct univariate meta-analysis (Schmidt & Hunter)
####################################
weightedmeans <- (sapply(data[,7:12],function(x) sum(x*data$N)/sum(data$N)))
vars <- sapply(7:12,function(i) {
  rs <- data[!is.na(data[,i]),i];
  ns <- data[!is.na(data[,i]),"N"];
  sum(ns*(rs-rep(weightedmeans[i-6],length(rs)))**2)/sum(ns)
}) # NOTE: I am using bare-bones Meta-Analysis on UNCORRECTED CORRELATIONS 
## (see the file; columns 7:12 are uncorrected, corrected correlations start in column 13)

aveve <- (1-weightedmeans**2)**2/(mean(data$N)-1)
sdrho <- sqrt(vars-aveve)

sapply(1:length(weightedmeans),function(i)
  {
  lb <- round(weightedmeans[i]-1.28*sdrho[i],2)
  ub <- round(weightedmeans[i]+1.28*sdrho[i],2)
  paste(lb,ub,sep="; ")
})

####################################
# conduct univariate FIMASEM
####################################
#first generate random matrices
#multivariate bootstrapping technique - specify inter-correlations to 0
sigma <- diag(sdrho**2)
#sigma <- cov(data[,7:12])
rownames(sigma) <- colnames(sigma) <- names(weightedmeans)
matrices.univ <- random.matrices.multivariate(1000,weightedmeans,sigma,varnames) #generate the matrices based on the point matrix

#run the SEM
coefs.fits.univariate.FIMASEM <- as.data.frame(t(sapply(1:1000,function(i) {
  openmxmodel <- mxModel("temp",mxData(matrices.univ[,,i],type="cov",numObs = mean(data$N)),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML());
  openmxfit <- mxRun(openmxmodel,silent=T);
  modelsummary <- summary(openmxfit);
  coefs <- mxStandardizeRAMpaths(openmxfit);
  output <- c(coefs[,8], modelsummary$CFI); 
  names(output) <- c(coefs[,2],"cfi")
  output
}))) 

## VISH AND ONES RESULTS
openmxmodel <- mxModel("temp",mxData(rho,type="cov",numObs = mean(data$N)),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML());
openmxfit <- mxRun(openmxmodel,silent=T);
modelsummary <- summary(openmxfit);
coefs <- mxStandardizeRAMpaths(openmxfit);
output <- c(coefs[,8], modelsummary$CFI); 
names(output) <- c(coefs[,2],"cfi")
modelsummary 
output

## FIMASEM RESULTS
sapply(coefs.fits.univariate.FIMASEM[1:8],mean)#MEAN BETA
sapply(coefs.fits.univariate.FIMASEM[1:8],function(x) {
  lb <- round(mean(x)-1.28*sd(x),2)
  ub <- round(mean(x)+1.28*sd(x),2)
  #lb-ub
  paste(lb,ub,sep="; ")
  })#80% CV BETA

####################################
# compose data for TSSEM and TS-FIMASEM
####################################
cormats <- list()
for (i in 1:nrow(data))
{
  temp <- diag(length(varnames));
  colnames(temp) <-  rownames(temp) <- varnames
  temp[lower.tri(temp)] <- as.numeric(data[i,c(11,9,7,8,10,12)]); # put in the lower half
  temp2 <- t(temp)
  temp[upper.tri(temp)] <- temp2[upper.tri(temp2)]
  rm(temp2)
  cormats[[i]] <- temp
}

####################################
# conduct TS-FIMASEM
####################################
step.one <- tssem1(cormats,data$N,method="REM",RE.type="Diag") #NOTE: HAD TO REMOVE 2 NPD STUDIES
rho.mult <- coef(step.one,select="fixed")
sigma.mult <- vcov(step.one,select="fixed")
sigma.mult.alt <- sigma.mult
diag(sigma.mult.alt) <- coef(step.one,select="random")
matrices.mult <- random.matrices.multivariate(500,rho.mult,sigma.mult.alt,varnames)

#run the SEM
coefs.fits.multivariate.FIMASEM <- as.data.frame(t(sapply(1:500,function(i) {
  openmxmodel <- mxModel("temp",mxData(matrices.mult[,,i],type="cov",numObs = mean(data$N)),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML());
  openmxfit <- mxRun(openmxmodel,silent=T);
  modelsummary <- summary(openmxfit);
  coefs <- mxStandardizeRAMpaths(openmxfit);
  output <- c(coefs[,8], modelsummary$CFI); 
  names(output) <- c(coefs[,2],"cfi")
  output
}))) 

## TS-FIMASEM RESULTS
sapply(coefs.fits.multivariate.FIMASEM[1:8],mean)#MEAN BETA
sapply(coefs.fits.multivariate.FIMASEM[1:8],function(x) {
  lb <- round(mean(x)-1.28*sd(x),2)
  ub <- round(mean(x)+1.28*sd(x),2)
  paste(lb,ub,sep="; ")
})#80% CV BETA

####################################'
# conduct TSSEM
####################################
step.two <- tssem2(step.one,A,S,Fmatrix=matrF)
summary(step.two)



