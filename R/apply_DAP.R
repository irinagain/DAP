#' Apply DAP for binary calssification 
#'
#' Apply sparse quadratic classification rules via linear dimension
#' reduction (projection matrix). Meanwhile, variable selection via
#' group lasso will be implemented. Can deal with high-dimensional
#' binary classification problem.
#'
#' @param xtrain Total training data.
#' @param ytrain Total training label, either "1" or "2".
#' @param xtest Test data.
#' @param ytest Test data label, either "1" or "2". Once it's provided,
#' this function will return a misclassification error rate on the test set;
#' otherwise, the function will return predicted labels for the test set. 
#' Default is NULL.
#' @param lambda_seq A sequence of tunning parameter, lambda. Deafult is NULL.
#' @param n_lambda Length of lambda_seq, used for generating lambda_seq if it's
#' NULL. Default is 50.
#' @param maxmin_ratio A ratio to control the minimum lambda in lambda_seq.
#' Default is 0.1.
#' @param nfolds Set folds number for cross-validation. Default is 5.
#' @param eps Convergence threshold for block-coordinate decent
#' algorithm. Each block-coordinate decent algorithm loop continuoues
#' until the maximum iteration number exceeds \code{maxiter} or the
#' maximum element-wise change in \eqn{V} is less than \code{eps}.
#' Default is 1e-4.
#' @param m_max Maximum number of iterations. Default is 10000.
#' @param myseed Seed for random spliting the data set into traininf
#' and testing. Default seed is 1001.
#' @param prior If "TRUE", the proportions for the training set will
#' be used to adjust the classification rule. Default is "TRUE".
#' 
#' @return A list as below.
#'        \item{error}{ Misclassification error rate if ytest is provided.}
#'        \item{ypred}{predicted label on the test set if ytest is NULL.}
#'        \item{features}{Number of selected features.}
#'        \item{feature_id}{Index of selected features.}
#' @details If no feature is selected by DAP, no matter ytest is NULL or 
#' provided, the function will return error = 0.5 and no ypred. In this case,
#' the classifier is no better than randomly guessing.
#' 
#' @example man/examples/apply_DAP_eg.R
#' 
#' @export
#' 
apply_DAP <- function(xtrain, ytrain, xtest, ytest = NULL, lambda_seq = NULL, n_lambda = 50,  maxmin_ratio = 0.1, nfolds = 5, eps = 1e-4, m_max = 10000, myseed = 1001, prior = TRUE){
  
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
  out.cv = cv_DAP(X = xtrain, Y = ytrain, lambda_seq = lambda_seq, nfolds = nfolds, eps = eps, m_max = m_max, myseed = myseed, prior = prior)
  
  ####solve for selected tuning parameter
  out.proj = solve_DAP_C(X1 = out_s$X1, X2 = out_s$X2, lambda = out.cv$lambda_min, eps = eps, maxiter = m_max)
  V = cbind(diag(1/out_s$coef1)%*%out.proj$V[,1],diag(1/out_s$coef2)%*% out.proj$V[,2])
  
  ### calculate the error and corresponding number of variables
  ####back scaling
  if (out.proj$nfeature > 0){
    ypred = classify_DAP(xtrain, ytrain, xtest, V = V, prior = prior)
    if (!is.null(ytest))
    {
      error = sum(ypred != ytest)/length(ytest)
      return(list(error = error, features = out.proj$nfeature, features_id = c(1:ncol(xtrain))[rowSums(abs(V)) != 0]))
    }else{
      return(list(ypred = ypred, features = out.proj$nfeature, features_id = c(1:ncol(xtrain))[rowSums(abs(V)) != 0]))
    }
  }else{
    return(list(error = 0.5, features = 0, features_id = NA))
  }
}
