#' Apply DAP for binary classification 
#'
#' Applies Discriminant Analysis via Projections to perform binary classification on the test dataset based on the training data. 
#'
#' @param xtrain A n x p training dataset; n observations on the rows and p features on the columns.
#' @param ytrain A n vector of training group labels, either 1 or 2.
#' @param xtest A m x p testing dataset; m observations on the rows and p features on the columns.
#' @param ytest An optional m vector of testing group labels, either 1 or 2. If supplied,
#' the function returns misclassification error rate;
#' if \code{NULL}, the function returns predicted labels for \code{xtest}. 
#' Default is \code{NULL}.
#' @param lambda_seq An optional sequence of tunning parameters lambda. Default is \code{NULL}, and the function generates its own sequence.
#' @param n_lambda Number of lambda values, the default is 50.
#' @param maxmin_ratio Smallest value for lambda, as a fraction of maximal value for which all coefficients are zero. The default is 0.1.
#' @param nfolds Number of folds for cross-validation, the default is 5.
#' @param eps Convergence threshold for the block-coordinate decent
#' algorithm based on the maximum element-wise change in \eqn{V}. The
#' default is 1e-4.
#' @param maxiter Maximum number of iterations, the default is 10000.
#' @param myseed Optional specification of random seed for generating the folds, the default value is 1001.
#' @param prior A logical indicating whether to put larger weights to the groups of larger size; the default value is \code{TRUE}.
#' 
#' @return A list of
#'        \item{error}{Misclassification error rate (if \code{ytest} is provided).}
#'        \item{ypred}{Predicted labels on the test set (if \code{ytest} is \code{NULL}).}
#'        \item{features}{Number of selected features.}
#'        \item{feature_id}{Index of selected features.}
#' @details If no feature is selected by DAP, the function will return \code{error} of 0.5 and no \code{ypred}, indicating that the classifier is no better than random guessing.
#' 
#' @example man/examples/apply_DAP_eg.R
#' 
#' @export
#' 
apply_DAP <- function(xtrain, ytrain, xtest, ytest = NULL, lambda_seq = NULL, n_lambda = 50,  maxmin_ratio = 0.1, nfolds = 5, eps = 1e-4, maxiter = 10000, myseed = 1001, prior = TRUE){
  
  Xmean = colMeans(xtrain)
  xtrain = xtrain - matrix(Xmean, nrow(xtrain), ncol(xtrain), byrow = T)
  xtest = xtest - matrix(Xmean, nrow(xtest), ncol(xtest), byrow = T)
  out_s = standardizeData(xtrain, ytrain, center = F)
  
  ## generate lambda sequence
  l_max = max(sqrt(colMeans(out_s$X1)^2 + colMeans(out_s$X2)^2))
  if (!is.null(lambda_seq)) {
    n_l = length(lambda_seq)
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
  out.cv = cv_DAP(X = xtrain, Y = ytrain, lambda_seq = lambda_seq, nfolds = nfolds, eps = eps, maxiter = maxiter, myseed = myseed, prior = prior)
  
  ####solve for selected tuning parameter
  out.proj = solve_DAP_C(X1 = out_s$X1, X2 = out_s$X2, lambda = out.cv$lambda_min, eps = eps, maxiter = maxiter)
  V = cbind(diag(1 / out_s$coef1) %*% out.proj$V[, 1], diag(1 / out_s$coef2) %*% out.proj$V[, 2])
  
  ### calculate the error and corresponding number of variables
  ####back scaling
  if (out.proj$nfeature > 0){
    ypred = classify_DAP(xtrain, ytrain, xtest, V = V, prior = prior)
    if (!is.null(ytest))
    {
      error = sum(ypred != ytest) / length(ytest)
      return(list(error = error, features = out.proj$nfeature, features_id = c(1:ncol(xtrain))[rowSums(abs(V)) != 0]))
    }else{
      return(list(ypred = ypred, features = out.proj$nfeature, features_id = c(1:ncol(xtrain))[rowSums(abs(V)) != 0]))
    }
  }else{
    return(list(error = 0.5, features = 0, features_id = NA))
  }
}
