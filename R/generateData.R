# Creates equicorrelation matrix
equicor <- function(p, rho, sblock = p){
  Sigma = matrix(0, p, p)
  Sigma[1:sblock, 1:sblock]=rho
  diag(Sigma)=1
  return(Sigma)
}

# Creates autocorrelation matrix
autocor <- function(p, rho, sblock = p){
  Sigma = matrix(0, p, p)
  for (i in 1:(sblock-1)){
    for (j in (i+1):sblock){
      Sigma[i,j]=rho^abs(i-j)
    }
  }
  Sigma = Sigma + t(Sigma)
  diag(Sigma)=1
  return(Sigma)
}

# Creates spiked matrix
spiked <- function(q1, q2){
  p = length(q1)
  Sigma = 30*tcrossprod(q1)+2*tcrossprod(q2)+ diag(p)
  return(Sigma)
}

# Create models for simulations
generateModel <- function(id, p){
  mu1 = rep(0,p)
  mu2 = c(rep(1,5), rep(-1,5), rep(0, p-10))
  if (id==1){
    # Both euqicorrelation 100
    Sigma1 = equicor(p, rho = 0.5, sblock = 100)
    Sigma2 = Sigma1
  }else if (id==2){
    # One equicorrelation 100, another autocorrelation 100
    Sigma1 = autocor(p, rho = 0.8, sblock = 100)
    Sigma2 = equicor(p, rho = 0.5, sblock = 100)
  }else if (id==3){
    # Both autocorrelation of size 10 with different parameters
    Sigma1 = autocor(p, rho = 0.5, sblock = 10)
    Sigma2 = autocor(p, rho = 0.8, sblock = 10)
  }else if (id==4){
    # One auto, one equicorrelation size 10 block
    Sigma1 = autocor(p, rho = 0.5, sblock = 10)
    Sigma2 = equicor(p, rho = 0.8, sblock = 10)
  }else if (id==5){
    # Both spiked size 10, one has reversed eigenvectors
    q1 = c(rep(1,5), rep(0, p-5))
    q2 = c(rep(0,5), rep(1,5), rep(0, p-10))
    q1 <- q1/sqrt(crossprod(q1))
    q2 <- q2/sqrt(crossprod(q2))
    Sigma1 = spiked(q1,q2)
    Sigma2 = spiked(q2,q1)
  }else if (id==6){
    # One spiked size 10, another spiked size 100
    q1 = c(rep(1,5), rep(0, p-5))
    q2 = c(rep(0,5), rep(1,5), rep(0, p-10))
    q1 <- q1/sqrt(crossprod(q1))
    q2 <- q2/sqrt(crossprod(q2))
    Sigma1 = spiked(q2,q1)
    q11 = c(c(1:100), rep(0, p-100))
    q11 <- q11/sqrt(crossprod(q11))
    q21 = c(c(100:1), rep(0, p-100))
    q21 = (diag(p) - tcrossprod(q11))%*%q21
    q21 <- q21/as.numeric(sqrt(crossprod(q21)))
    Sigma2 = spiked(q11,q21)
  }else if (id==7){
    # One spiked size 10, another equicorrelation size 10
    q1 = c(rep(1,5), rep(0, p-5))
    q2 = c(rep(0,5), rep(1,5), rep(0, p-10))
    q1 <- q1/sqrt(crossprod(q1))
    q2 <- q2/sqrt(crossprod(q2))
    Sigma1 = spiked(q1,q2)
    Sigma2 = equicor(p, rho = 0.8, sblock = 10)
  }else if (id ==8){
    # One spiked size 10, another equicorrelation size 100
    q1 = c(rep(1,5), rep(0, p-5))
    q2 = c(rep(0,5), rep(1,5), rep(0, p-10))
    q1 <- q1/sqrt(crossprod(q1))
    q2 <- q2/sqrt(crossprod(q2))
    Sigma1 = spiked(q1,q2)
    Sigma2 = equicor(p, rho = 0.3, sblock = 100)
  }else if (id == 9){
    # One spiked size 100, another equicorrelation size 100
    q11 = c(c(1:100), rep(0, p-100))
    q11 <- q11/sqrt(crossprod(q11))
    q21 = c(c(100:1), rep(0, p-100))
    q21 = (diag(p) - tcrossprod(q11))%*%q21
    q21 <- q21/as.numeric(sqrt(crossprod(q21)))
    Sigma1 = spiked(q11,q21)
    Sigma2 = equicor(p, rho = 0.3, sblock = 100)
  }else{
    stop("Unknown model id")
  }
  return(list = list(mu1 = mu1, mu2 = mu2, Sigma1 = Sigma1, Sigma2 = Sigma2))
}


# Generate data for a given model
# model - output from generateModel, list
generateData <- function(model, n1, n2, n_test){
  #training data
  ytrain<-c(rep(1,n1),rep(2,n2))
  x1=mvrnorm(n=n1,mu=model$mu1,Sigma=model$Sigma1)
  x2=mvrnorm(n=n2,mu=model$mu2,Sigma=model$Sigma2)
  xtrain<-rbind(x1,x2)
  #test data
  x1_test=mvrnorm(n=n_test,mu=model$mu1,Sigma=model$Sigma1)
  x2_test=mvrnorm(n=n_test,mu=model$mu2,Sigma=model$Sigma2)
  xtest<-rbind(x1_test,x2_test)
  ytest<-c(rep(1,n_test),rep(2,n_test))
  return(list(xtrain=xtrain, ytrain = ytrain, xtest = xtest, ytest = ytest))
}