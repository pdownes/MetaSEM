# Simulation Study 1 (Model 2)
# The script requires the "random.matrices.R" script, which defines the functions to generate random matrices.
# Written by Patrick Downes, Joya Yu, Kameron Carter, and Ernest O'Boyle. Correspondence to Joya Yu.


########################################
# Define Relevant Inputs
########################################
rho <- .3
sdrho <- .1
k <- 50 #number of primary studies. doesn't affect outcomes.
generated_k <- 200 #due to the more complex model, you often need to generate more than k studies and eliminate non-positive definite data in order to run TSSEM. 3-4x k seems to work okay.
N <- 10000 #total sample size (i.e., sum of n from 1 to k).


########################################
# Run the simulation!
# Note: You shouldn't need to edit below this line
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
data <- data.frame(rnorm(generated_k,0,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,0,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,rho,obs_sd),rnorm(generated_k,0,obs_sd)) #generate data
names(data) <- c('x1x2','x1x3','x1m1','x1m2','x1y1','x1y2','x2x3','x2m1','x2m2','x2y1','x2y2','x3m1','x3m2','x3y1','x3y2','m1m2','m1y1','m1y2','m2y1','m2y2','y1y2') #add names to dataframe

#Reformat the data for TSSEM input as a list of matrices
test_cormats <- list()
for (i in 1:nrow(data)){
  temp <- matrix(1,nrow=7,ncol=7)
  temp[lower.tri(temp)] <- as.numeric(data[i,])
  temp2 <- t(temp)
  temp[upper.tri(temp)] <- temp2[upper.tri(temp2)]
  test_cormats[[i]] <- temp
}

data <- data[sapply(test_cormats,is.positive.definite),] #eliminate NPD matrices
data <- data[1:50,] # use only the first 50 positive definite matrices

cormats <- list()
for (i in 1:nrow(data)){
  temp <- matrix(1,nrow=7,ncol=7)
  temp[lower.tri(temp)] <- as.numeric(data[i,])
  temp2 <- t(temp)
  temp[upper.tri(temp)] <- temp2[upper.tri(temp2)]
  cormats[[i]] <- temp
}

####################################
# Specify the Model using OpenMx specification
# See http://openmx.psyc.virginia.edu/docs/OpenMx/2.3.1/Regression_Matrix.html
####################################
varnames <- c('X1','X2','X3','M1','M2','Y1','Y2')
A <- mxMatrix('Full',ncol=7,nrow=7,byrow=T,
              values = c(0,0,0,rho,rho,0,0,
                         0,0,0,rho,rho,0,0,
                         0,0,0,rho,rho,0,0,
                         0,0,0,0,0,rho,rho,
                         0,0,0,0,0,rho,rho,
                         0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0),
              free=c(F,F,F,T,T,F,F,
                     F,F,F,T,T,F,F,
                     F,F,F,T,T,F,F,
                     F,F,F,F,F,T,T,
                     F,F,F,F,F,T,T,
                     F,F,F,F,F,F,F,
                     F,F,F,F,F,F,F
              ),
              labels=c(NA,NA,NA,"betax1m1","betax1m2",NA,NA,
                       NA,NA,NA,"betax2m1","betax2m2",NA,NA,
                       NA,NA,NA,"betax3m1","betax3m2",NA,NA,
                       NA,NA,NA,NA,NA,"betam1y1","betam1y2",
                       NA,NA,NA,NA,NA,"betam2y1","betam2y2",
                       NA,NA,NA,NA,NA,NA,NA,
                       NA,NA,NA,NA,NA,NA,NA
              ),
              name="A")
S <- mxMatrix('Full',ncol=7,nrow=7,byrow=T,
              values = c(1,0,0,0,0,0,0,
                         0,1,0,0,0,0,0,
                         0,0,1,0,0,0,0,
                         0,0,0,1,0,0,0,
                         0,0,0,0,1,0,0,
                         0,0,0,0,0,1,0,
                         0,0,0,0,0,0,1
                         ),
              free=c(F,F,F,F,F,F,F,
                     F,F,F,F,F,F,F,
                     F,F,F,F,F,F,F,
                     F,F,F,F,T,F,F,
                     F,F,F,T,F,F,F,
                     F,F,F,F,F,F,T,
                     F,F,F,F,F,T,F
                     ),
              labels=c("varx1",NA,NA,NA,NA,NA,NA,
                       NA,"varx2",NA,NA,NA,NA,NA,
                       NA,NA,"varx3",NA,NA,NA,NA,
                       NA,NA,NA,"varm1","covm1m2",NA,NA,
                       NA,NA,NA,"covm1m2","varm2",NA,NA,
                       NA,NA,NA,NA,NA,"vary1","covy1y2",
                       NA,NA,NA,NA,NA,"covy1y2","vary2"
              ),
              name="S")
matrF <- mxMatrix(type="Iden",nrow=,7,ncol=7,name="F")
exp         <- mxExpectationRAM("A","S","F", dimnames=varnames )


####################################
# conduct univariate meta (Hunter & Schmidt)
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

aveve.denom <- rep(N-1,ncol(data)) #begin computing average sampling error variance
aveve.numerator <- sapply(weightedmeans,function(x)
  ((1-x**2)**2)
) # using the weighted mean in the numerator (instead of each study r) is more accurate (see Hunter & Schmidt)

aveve <- aveve.numerator/aveve.denom
sdrho <- sqrt(vars-aveve) # subtract the sampling variances of each effect from the observed variances of each effect to get var(rho))


####################################
# conduct univariate FIMASEM
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

#sum(sapply(1:500, function(i) is.positive.definite(matrices.univ[,,i]))) # how many matrices are positive definite? if there are many npd matrices, consider using nearPD() on all matrices

#run the SEM across all the random matrices
coefs.fits.univariate.FIMASEM <- as.data.frame(t(sapply(1:500,function(i) {
  openmxmodel <- mxModel("temp",mxData(matrices.univ[,,i],type="cov",numObs = N/k),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML());
  openmxfit <- mxRun(openmxmodel,silent=T);
  modelsummary <- summary(openmxfit);
  coefs <- mxStandardizeRAMpaths(openmxfit);
  output <- c(coefs[,8], modelsummary$CFI,modelsummary$Chi,modelsummary$RMSEA,fitMeasuresMx(openmxfit)[[27]]); 
  names(output) <- c(coefs[,2],"cfi","chisq",'RMSEA','srmr')
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
  openmxmodel <- mxModel("temp",mxData(matrices.mult[,,i],type="cov",numObs = N),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML());
  openmxfit <- mxRun(openmxmodel,silent=T);
  modelsummary <- summary(openmxfit);
  coefs <- mxStandardizeRAMpaths(openmxfit);
  output <- c(coefs[,8], modelsummary$CFI,modelsummary$Chi,modelsummary$RMSEA,fitMeasuresMx(openmxfit)[[27]]); 
  names(output) <- c(coefs[,2],"cfi","chisq",'RMSEA','srmr')
  output
}))) #returns a dataframe of SEM parameter estimates (i.e., fit indices and path coefficients)

####################################'
# conduct second stage of TSSEM
####################################
step.two <- tssem2(step.one,A,S,Fmatrix=matrF)

####################################'
# Get Results
####################################
print('## VISH AND ONES RESULTS')
openmxmodel <- mxModel("temp",mxData(rho,type="cov",numObs = N/k),matrA = A,matrS = S,matrF=matrF,exp=exp,funML=mxFitFunctionML());
openmxfit <- mxRun(openmxmodel,silent=T);
modelsummary <- summary(openmxfit);
coefs <- mxStandardizeRAMpaths(openmxfit);
output <- c(coefs[,8], modelsummary$CFI,modelsummary$Chi,modelsummary$RMSEA,fitMeasuresMx(openmxfit)[[27]]) 
names(output) <- c(coefs[,2],"cfi",'chisq','RMSEA','srmr')
output

print('## UNIVARIATE FIMASEM RESULTS')
print('means')
sapply(coefs.fits.univariate.FIMASEM,mean)#MEANS
print('80% CVs')
sapply(coefs.fits.univariate.FIMASEM,function(x) {
  lb <- round(mean(x)-1.28*sd(x),2)
  ub <- round(mean(x)+1.28*sd(x),2)
  paste(lb,ub,sep="; ")
})#80% CVs
print('% SRMR < .10')
sum(coefs.fits.univariate.FIMASEM$srmr < .10)/length(coefs.fits.univariate.FIMASEM$srmr)# %srmr < .10
sum(coefs.fits.univariate.FIMASEM$cfi > .90)/length(coefs.fits.univariate.FIMASEM$cfi)# %cfi > .90

print('## MULTIVARIATE FIMASEM RESULTS')
print('means')
sapply(coefs.fits.multivariate.FIMASEM,mean)#MEANS
print('80% CVs')
sapply(coefs.fits.multivariate.FIMASEM,function(x) {
  lb <- round(mean(x)-1.28*sd(x),2)
  ub <- round(mean(x)+1.28*sd(x),2)
  paste(lb,ub,sep="; ")
})#80% CVs
print('% SRMR < .10')
sum(coefs.fits.multivariate.FIMASEM$srmr < .10)/length(coefs.fits.multivariate.FIMASEM$srmr)# %srmr < .10
sum(coefs.fits.multivariate.FIMASEM$cfi > .90)/length(coefs.fits.multivariate.FIMASEM$cfi)# %cfi > .90

print('## TSSEM RESULTS')
summary(step.two)


