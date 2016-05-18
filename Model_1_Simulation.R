# Simulation Study 1 (Model 1)
# The script requires the "random.matrices.R" script, which defines the functions to generate random matrices.
# Written by Patrick Downes, Joya Yu, Kameron Carter, and Ernest O'Boyle. Correspondence to Joya Yu.

########################################
# Define Relevant Inputs
########################################
rho <- .3
sdrho <- .1
k <- 50 #number of primary studies. doesn't affect outcomes.
N <- 1000 #total sample size (i.e., sum of n from 1 to k).


########################################
# Run the simulation!
# Note: You shouldn't need to edit below this line to run the Model 1 simulations
########################################
#load libraries
library(metaSEM)
require('matrixcalc')
require('semTools')
library(Matrix)
source('Dropbox/Meta SEM/random.matrices.R') #updated as of 12-14-15

####################################
# Generate Data (Basically a Meta-Analytic Coding Sheet)
####################################
obs_sd <- sqrt(sdrho**2 + (((1-rho**2)**2)/(N-1))) #need to add average sampling variance to sdrho in order to get observed sd
data <- data.frame(rnorm(k,rho,obs_sd),rnorm(k,rho,obs_sd),rnorm(k,0,obs_sd),rnorm(k,0,obs_sd),rnorm(k,rho,obs_sd),rnorm(k,rho,obs_sd)) #generate data
names(data) <- c('xm1','xm2','xy','m1m2','m1y','m2y') #add names to data

#Reformat the data for TSSEM input - generating a list of correlation matrices
cormats <- list()
for (i in 1:nrow(data)){
  temp <- matrix(1,nrow=4,ncol=4)
  temp[lower.tri(temp)] <- as.numeric(data[i,])
  temp2 <- t(temp)
  temp[upper.tri(temp)] <- temp2[upper.tri(temp2)]
  cormats[[i]] <- temp
}

# If a matrix is not positive definte, replace it with the nearest positive definite matrix
# NOTE: For certain combinations of k, rho, and sdrho, this is not enough to handle all NPD
## issues. In these cases, you will need to generate more than k studies (e.g., 5 X k studies)
## and only use the positive definite matrices.
for(i in 1:k){
  if (!is.positive.definite(cormats[[i]])) {
    cormats[[i]] <- as.matrix(nearPD(cormats[[i]],corr=T,keepDiag=T)$mat)
    data[i,] <- as.numeric(cormats[[i]][lower.tri(cormats[[i]])])
  }
}

####################################
# Specify the Model using OpenMx specification
# See http://openmx.psyc.virginia.edu/docs/OpenMx/2.3.1/Regression_Matrix.html
####################################
varnames <- c('X','M1','M2','Y')
A <- mxMatrix('Full',ncol=4,nrow=4,byrow=T,
              values = c(0,rho,rho,0,
                         0,0,0,rho,
                         0,0,0,rho,
                         0,0,0,0),
              free=c(F,T,T,F,
                     F,F,F,T,
                     F,F,F,T,
                     F,F,F,F
              ),
              labels=c(NA,"betaxm1","betaxm2",NA,
                       NA,NA,NA,"betam1y",
                       NA,NA,NA,"betam2y",
                       NA,NA,NA,NA
              ),
              name="A")
S <- mxMatrix('Full',ncol=4,nrow=4,byrow=T,
              values = c(1,0,0,0,
                         0,1,.2,0,
                         0,.2,1,0,
                         0,0,0,1),
              free=c(F,F,F,F,
                     F,F,T,F,
                     F,T,F,F,
                     F,F,F,F),
              labels=c("varx",NA,NA,NA,
                       NA,"varm1","covm1m2",NA,
                       NA,"covm1m2","varm2",NA,
                       NA,NA,NA,"vary"
              ),
              name="S")
matrF <- mxMatrix(type="Iden",nrow=4,ncol=4,name="F")
exp         <- mxExpectationRAM("A","S","F", dimnames=varnames )


####################################
# Conduct univariate meta-analysis (Hunter & Schmidt)
####################################
weightedmeans <- sapply(data,function(x) {
  numerator <- sum(x*(N/k),na.rm=T);
  denominator <- N
  numerator/denominator
}) # a vector of all effect sizes; note this assumes effects have been individually corrected for artifacts

vars <- sapply(1:ncol(data),function(i) {
  rs <- data[!is.na(data[,i]),i];
  ns <- rep(N/k,length(rs))
  sum(ns*(rs-rep(weightedmeans[i],length(rs)))**2)/N
}) # a vector of observed variances

aveve.denom <- rep(N-1,ncol(data)) 
aveve.numerator <- sapply(weightedmeans,function(x)
  ((1-x**2)**2)
) # using the weighted mean in the numerator (instead of each study r) is more accurate (see Hunter & Schmidt)

aveve <- aveve.numerator/aveve.denom #compute average sampling error variance
sdrho <- sqrt(vars-aveve) # subtract the sampling variances of each effect from the observed variances of each effect to get var(rho))


####################################
# Conduct Univariate FIMASEM
####################################
#first generate random matrices
#multivariate bootstrapping technique, specifying intercorrelations to zero
rho <- diag(length(varnames)) #create a rho matrix
rho[lower.tri(rho)] <- weightedmeans #input the meta-results into the rho matrix
temp <- t(rho) #symmetrize the rho matrix (part 1)
rho[upper.tri(rho)] <- temp[upper.tri(temp)]#symmetrize the rho matrix (part 2)
rm(temp)
rownames(rho) <- colnames(rho) <-  varnames #ensure names match from above
sigma <- diag(sdrho**2) #the sigma matrix is essentially the variance-covariance matrix of the effect sizes; we assume covariances are all zero, and we put var(rho) along the diagonal; when random matrices are bootstrapped we'll use sigma
rownames(sigma) <- colnames(sigma) <- names(weightedmeans) #make the names match
matrices.univ <- random.matrices.multivariate(500,weightedmeans,sigma,varnames) #generate the matrices based on the rho matrix

#run the SEM across all the random matrices
coefs.fits.univariate.FIMASEM <- as.data.frame(t(sapply(1:500,function(i) {
  openmxmodel <- mxModel("temp",mxData(matrices.univ[,,i],type="cov",numObs = N/k),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML()); #specify the model and data
  openmxfit <- mxRun(openmxmodel,silent=T);#estimate the SEM
  modelsummary <- summary(openmxfit); #retrieve outputs
  coefs <- mxStandardizeRAMpaths(openmxfit); #standardize the coefficients
  output <- c(coefs[,8], modelsummary$CFI,modelsummary$Chi,modelsummary$RMSEA); #combine the path coefficients and fit statistics
  names(output) <- c(coefs[,2],"cfi","chisq",'RMSEA')
  output
}))) #returns a dataframe of SEM parameter estimates (i.e., fit indices and path coefficients)

####################################
# conduct multivariate FIMASEM
####################################
step.one <- tssem1(cormats,rep(N/k,k),method="REM",RE.type="Diag") #run TSSEM Step 1
rho.mult <- coef(step.one,select="fixed")
sigma.mult.alt <- diag(coef(step.one,select="random")) # the sigma matrix from this analysis will have Tau^2 along the diagonal and zeros in the off-diagonal. There are a number of possible approaches to this sigma matrix; we think this is most appropriate.
#sigma.mult.alt <- nearPD(sigma.mult.alt)$mat #defining the matrix this leaves the possbility for sigma to be NPD. We use nearPD() in this case.
matrices.mult <- random.matrices.multivariate(500,rho.mult,sigma.mult.alt,varnames) #generate the random matrices

# Run SEM on those random matrices using the same technique we use for V&O FIMASEM
coefs.fits.multivariate.FIMASEM <- as.data.frame(t(sapply(1:nrow(data),function(i) {
  openmxmodel <- mxModel("temp",mxData(matrices.mult[,,i],type="cov",numObs = N),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML()); #specify the model and data
  openmxfit <- mxRun(openmxmodel,silent=T);#estimate the SEM
  modelsummary <- summary(openmxfit);#retrieve outputs
  coefs <- mxStandardizeRAMpaths(openmxfit);#standardize the coefficients
  output <- c(coefs[,8], modelsummary$CFI,modelsummary$Chi,modelsummary$RMSEA); #combine the path coefficients and fit statistics
  names(output) <- c(coefs[,2],"cfi","chisq",'RMSEA')
  output
}))) #returns a dataframe of SEM parameter estimates (i.e., fit indices and path coefficients)

####################################'
# conduct second stage of TSSEM
####################################
step.two <- tssem2(step.one,A,S,Fmatrix=matrF) 

####################################'
# Get Results from all techniques
####################################
## VISH AND ONES RESULTS
openmxmodel <- mxModel("temp",mxData(rho,type="cov",numObs = N/k),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML());
openmxfit <- mxRun(openmxmodel,silent=T);
modelsummary <- summary(openmxfit);
coefs <- mxStandardizeRAMpaths(openmxfit);
output <- c(coefs[,8], modelsummary$CFI); 
names(output) <- c(coefs[,2],"cfi")

## UNIVARIATE FIMASEM RESULTS
sapply(coefs.fits.univariate.FIMASEM[1:9],mean)#MEAN BETA
sapply(coefs.fits.univariate.FIMASEM[1:9],function(x) {
  lb <- round(mean(x)-1.28*sd(x),2)
  ub <- round(mean(x)+1.28*sd(x),2)
  paste(lb,ub,sep="; ")
})#80% CV BETA

## MULTIVARIATE FIMASEM RESULTS
sapply(coefs.fits.multivariate.FIMASEM[1:9],mean)#MEAN BETA
sapply(coefs.fits.multivariate.FIMASEM[1:9],function(x) {
  lb <- round(mean(x)-1.28*sd(x),2)
  ub <- round(mean(x)+1.28*sd(x),2)
  paste(lb,ub,sep="; ")
})#80% CV BETA

## TSSEM RESULTS
summary(step.two)

cor(data)#"within-study" correlations
