#' Standardize the data
#'
#' Given matrix \eqn{X} and class-indicator \eqn{Y}, the following function centers \eqn{X}, and then forms \eqn{X1} and \eqn{X2} which are scaled, as well as coefficients for back-scaling.
#'
#' @param X A n by p matrix of n observations with p measurements per each.
#' @param Y The class-indicator denoting a binary label for each observation.
#' @param center wheather \eqn{X} will be centered or not. The default is "TURE".
#'
#' @return A list with components:
#'     \item{X1}{Standardized observations labeled in group 1.}
#'     \item{X2}{Standardized observations labeled in group 2.}
#'     \item{coef1}{Coefficient of \eqn{X1} for back-scaling.}
#'     \item{coef2}{Coefficient of \eqn{X2} for back-scaling.}
#'     \item{Xmean}{Column mean for the matrix \eqn{X}.}
#'
#' @example man/examples/standardizeData_eg.R
#'
#' @export
standardizeData <- function(X, Y, center = T){
  # center X
  Xmean = colMeans(X)
  if (center){
    X <- X - matrix(Xmean, nrow(X), ncol(X), byrow = T)
  }
  X1 = X[Y==1,]
  X2 = X[Y==2,]
  coef1 = sqrt(colSums(X1^2)/sum(Y==1))
  coef2 = sqrt(colSums(X2^2)/sum(Y==2))
  X1 = X1/ matrix(coef1, nrow(X1), ncol(X1), byrow = T)
  X2 = X2/ matrix(coef2, nrow(X2), ncol(X2), byrow = T)
  return(list(X1 = X1, X2 = X2, coef1 = coef1, coef2 = coef2, Xmean = Xmean))
}

#' Solve Optimization Problem (C version, fixed lambda)
#'
#' Solving group lasso using block-coordinate decent algorithm for a
#' fixed value of lambda. C code is used for coding basic functions
#' to speed up.
#'
#' @param X1 A matrix (n_1 by p) of group 1 data (scaled).
#' @param X2 A matrix (n_2 by p) of group 2 data (scaled).
#' @param lambda A fixed value of the tuning parameter.
#' @param Vinit Starting point. The default is "NULL".
#' @param eps Convergence threshold for block-coordinate decent
#' algorithm. Each block-coordinate decent algorithm loop continuoues
#' until the maximum iteration number exceeds \code{maxiter} or the
#' maximum element-wise change in \eqn{V} is less than \code{eps}.
#' Default is 1e-02.
#' @param maxiter Maximum number of iterations. Default is 10000.
#'
#' @return A list with following components:
#'        \item{V}{The projection matrix.}
#'        \item{nfeature}{The number of selected features.}
#'        \item{iter}{Iteration numbers used to converge.}
#'
#' @section Warnings:
#' Please take scaled X1 and X2 for this function. One can use
#' function "standardizeData" to do so.
#'
#' @example man/examples/solve_DAP_C
#'
#' @export
solve_DAP_C <-function(X1, X2, lambda, Vinit = NULL, eps = 1e-02, maxiter = 10000){
  p = ncol(X1)
  n1 = nrow(X1)
  n2 = nrow(X2)
  if (is.null(Vinit)){
    Vinit=matrix(0,p,2)
  }
  nitr = 0
  out =.C("solveProj_withS", as.double(X1), as.double(X2), as.double(Vinit), as.double(lambda), as.integer(p),as.integer(n1), as.integer(n2), as.double(eps), as.integer(maxiter), as.integer(nitr))
  if (out[[10]]==maxiter){
    warning(paste("Projections didn't converge in ", maxiter, " iterations, try increasing the number of iterations or decreasing the convergence level.", sep=""))
  }
  V = matrix(out[[3]], p, 2)
  nfeature = sum(rowSums(abs(V))>0)

  return(list(V=V, nfeature = nfeature, iter = out[[10]]))
}

####This function is used to solve optimization problem with  a sequence of lambda#####
#' Title
#'
#' @param X1
#' @param X2
#' @param lambda_seq
#' @param eps
#' @param m_max
#' @param feature_max
#'
#' @return
#' @export
#'
#' @examples
solve_DAP_seq <- function(X1, X2, lambda_seq, eps = 1e-2, m_max = 10000, feature_max = nrow(X1) + nrow(X2)){
  p =ncol(X1)
  n_lambda = length(lambda_seq)

  ####initilize V1_mat, V2_mat, both p by n_lambda, V0
  V1_mat = V2_mat = matrix(NA, nrow = p, ncol = n_lambda)
  V0 = matrix(0, nrow = p, ncol = 2)
  nfeature_vec = rep(NA, n_lambda)
  for (i in 1:n_lambda){
    ####use solve_proj for each lambda
    out = solve_DAP_C(X1, X2, lambda = lambda_seq[i], Vinit = V0, eps = eps, maxiter = m_max)
    ####use V0 from previous iteration as a warm start
    V0 = out$V
    V1_mat[, i] = V0[, 1]
    V2_mat[, i] = V0[, 2]
    nfeature_vec[i] = out$nfeature
    # The number of nonzero variables exceed the set threshold in feature_max
    if (out$nfeature >= feature_max){
      return (list(V1_mat = V1_mat[,1:i], V2_mat = V2_mat[,1:i], lambda_seq = lambda_seq[1:i], nfeature_vec = nfeature_vec[1:i]))
    }
  }
  return (list(V1_mat = V1_mat, V2_mat = V2_mat, lambda_seq = lambda_seq, nfeature_vec = nfeature_vec))
}

#Classify using 2-dimensional subspace formed by V
classify_DAP <- function(xtrain, ytrain, xtest, V, prior = TRUE){
  #xtrain = scale(xtrain, scale = F)
  #xtest = scale(xtest, scale = F, center = attr(xtrain, which = "scaled:center"))
  Vsvd <- svd(V)$d
  ## Only for 2 column case for now, and 2 classes
  if ((min(Vsvd)<1e-6)|(length(Vsvd)==1)){
    V <- V[,1, drop=F]
    trainproj <- xtrain%*%V
    testproj <- xtest %*%V
    Omega1 <- 1/var(trainproj[ytrain==1,])
    Omega2 <- 1/var(trainproj[ytrain==2,])
    mu1 <- mean(trainproj[ytrain==1,])
    mu2 <- mean(trainproj[ytrain==2,])
    predict_vec <- rep(2,nrow(xtest))
    difference <- diag(tcrossprod(testproj)*(Omega1 - Omega2))-2*as.vector(testproj*(Omega1*mu1 - Omega2*mu2))-mu2^2*Omega2+mu1^2*Omega1-log(Omega1)+log(Omega2)
    if (prior == TRUE){
      difference <- difference - 2*log(sum(ytrain == 1)/sum(ytrain == 2))
    }
    predict_vec[difference < 0] <- 1
    return(predict_vec)
  }else{
    # Use qda from MASS for prediction
    if (prior == T){
      out = qda(xtrain %*% V, grouping = ytrain)
    }else{
      out = qda(xtrain %*% V, grouping = ytrain, prior = c(1/2,1/2))
    }
    pred = predict(out, newdata = xtest %*% V)

    return(as.numeric(pred$class))
  }
}

# CV using projections as before
cv_DAP <-function(X, Y, lambda_seq, nfolds = 5, rho = 0, gamma1 = 0, gamma2 = 0, eps = 1e-6, m_max = 1000, myseed = 1001, prior = TRUE){

  n = length(Y)
  n_lambda = length(lambda_seq)

  error_mat = matrix(0.5, nfolds, n_lambda)
  #cor_mat = matrix(0, nfolds, n_lambda)
  nfeature_mat = matrix(NA, nfolds, n_lambda)

  ####random set split the whole data set into k folds corresponding to n1 and n2
  set.seed(myseed)
  id = rep(NA, n)
  id[Y==1]<- sample(rep(seq_len(nfolds), length.out = sum(Y==1)))
  id[Y==2]<- sample(rep(seq_len(nfolds), length.out = sum(Y==2)))

  for (nf in 1: nfolds){
    cat(nf)
    ####set training data and test data
    xtrain = X[id != nf, ]
    ytrain = Y[id != nf]
    xtest = X[id == nf, ]
    ytest = Y[id == nf]
    out_s <- standardizeData(xtrain, ytrain, center = T)

    # Calculate mean difference based on test data
    d_test = colMeans(xtest[ytest==1,])-colMeans(xtest[ytest==2,])

    ####use solve_proj_seq
    fit_tmp = solve_DAP_seq(X1 = out_s$X1, X2 = out_s$X2, lambda_seq = lambda_seq, eps = eps, m_max = m_max, feature_max = n)
    nfeature_mat[nf,1:length(fit_tmp$nfeature_vec)] = fit_tmp$nfeature_vec

    #### Calculate errors for each lambda
    for (j in 1:length(fit_tmp$lambda_seq)){
      V = cbind(diag(1/out_s$coef1)%*%fit_tmp$V1_mat[,j],diag(1/out_s$coef2)%*% fit_tmp$V2_mat[,j])
      ypred = classify_DAP(xtrain- matrix(out_s$Xmean, nrow(xtrain), ncol(xtest), byrow = T), ytrain, xtest = xtest - matrix(out_s$Xmean, nrow(xtest), ncol(xtest), byrow = T), V, prior = prior)
      error_mat[nf, j] = sum(ypred != ytest)/length(ytest)
    }
  }
  ####calculate cvm, cvse, lambda_min, lambda_1se
  cvm = colMeans(error_mat)
  index = which.min(cvm)
  lambda_min = lambda_seq[index]

  cvse = apply(error_mat, 2, sd)/sqrt(nfolds)
  se_lambda = cvse[index]
  up_b = cvm[index] + se_lambda
  lambda_1se = (lambda_seq[cvm <= up_b])[1]


  return (list(lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse, lambda_seq = lambda_seq, nfeature_mat = nfeature_mat, error_mat = error_mat))
}

# Apply DAP
apply_DAP <- function(xtrain, ytrain, xtest, ytest, lambda_seq = NULL, n_lambda = 50,  maxmin_ratio = 0.1, nfolds = 5, eps = 1e-4, m_max = 10000, myseed = 1001, prior = TRUE){

  Xmean = colMeans(xtrain)
  xtrain <- xtrain - matrix(Xmean, nrow(xtrain), ncol(xtrain), byrow = T)
  xtest <- xtest - matrix(Xmean, nrow(xtest), ncol(xtest), byrow = T)
  out_s <- standardizeData(xtrain, ytrain, center = F)

  ## generate lambda sequence
  l_max <- max(sqrt(colMeans(out_s$X1)^2 + colMeans(out_s$X2)^2))
  if (!is.null(lambda_seq)) {
    n_l= length(lambda_seq)
    if (nl < 1){
      warning(paste("There is no qualified lambda value. New values will be generated automatically. n_lambda will be set as.", n_lambda,sep = " "))
      lambda_seq = exp(seq(log(l_max * maxmin_ratio), log(l_max), length.out = n_lambda))
    }else{
      n_lambda = nl
    }
  }else {
    lambda_seq = exp(seq(log(l_max * maxmin_ratio), log(l_max), length.out = n_lambda))
  }
  lambda_seq = sort(lambda_seq, decreasing = TRUE)

  ####use cv to select the tuning parameter
  out.cv = cv_DAP(X = xtrain, Y = ytrain, lambda_seq = lambda_seq, gamma1 = 0, gamma2 = 0, nfolds = nfolds, rho = rho, eps = eps, m_max = m_max, myseed = myseed, prior = prior)

  ####solve for selected tuning parameter
  out.proj = solve_DAP_C(X1 = out_s$X1, X2 = out_s$X2, lambda = out.cv$lambda_min, eps = eps, maxiter = m_max)
  V = cbind(diag(1/out_s$coef1)%*%out.proj$V[,1],diag(1/out_s$coef2)%*% out.proj$V[,2])

  ### calculate the error and corresponding number of variables
  ####back scaling
  if (out.proj$nfeature > 0){
    ypred = classify_DAP(xtrain, ytrain, xtest, V = V, prior = prior)
    error = sum(ypred != ytest)/length(ytest)
    return(list(error = error, features = out.proj$nfeature, features_id = c(1:ncol(xtrain))[rowSums(abs(V)) != 0]))
  }else{
    return(list(error = 0.5, features = 0, features_id = NA))
  }
}
