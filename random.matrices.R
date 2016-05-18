random.matrices.univariate <- function(point.matrix,reps)
{
  build.matrix <- function()
  {
    sigma <- diag(point.matrix[,,2][lower.tri(point.matrix[,,2])])
    vec <- mvrnorm(1,point.matrix[,,1][lower.tri(point.matrix[,,1])],Sigma=sigma)
    mat <- diag(1,nrow(point.matrix))
    mat[lower.tri(mat)] <- vec
    trans.mat <- t(mat)
    mat[upper.tri(mat)] <- trans.mat[upper.tri(trans.mat)]
    if (!is.positive.definite((mat)))
    {
      return(build.matrix())
    }else{
      return(mat)
    }
  }
  output <- array(1,dim = c(nrow(point.matrix),ncol(point.matrix),reps))
  for (rep in 1:reps)
  { 
    output[,,rep] <- build.matrix()
  }
  rownames(output) <- colnames(output) <- rownames(point.matrix)
  return(output)
}


random.matrices.multivariate <- function(reps,estimates,sigma,varnames)
{
  temp.data <- mvrnorm(reps,estimates,Sigma=sigma)
  output <- array(NA,dim = c(length(varnames),length(varnames),reps))
  for (i in 1:reps)
  {
    temp <- diag(length(varnames));
    temp[lower.tri(temp)] <- temp.data[i,];
    temp2 <- t(temp);
    temp[upper.tri(temp)] <- temp2[upper.tri(temp2)];
     if (is.positive.definite(temp)){
       output[,,i] <- temp;
     }else{
       output[,,i] <- as.matrix(nearPD(temp,keepDiag = T)$mat);
     }
  }
  rownames(output) <- varnames
  colnames(output) <- varnames
  return(output)
}