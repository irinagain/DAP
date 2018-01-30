#' Cross-validation for DAP
#'
#' Chooses optimal tuning parameter lambda for DAP based on the k-fold cross-validation to minimize the misclassification error rate
#'
#' @param X A n x p training dataset; n observations on the rows and p features on the columns.
#' @param Y A n vector of training group labels, either 1 or 2.
#' @param lambda_seq A sequence of tuning parameters to choose from.
#' @param nfolds Number of folds for cross-validation, the default is 5.
#' @param eps Convergence threshold for the block-coordinate decent algorithm based on the maximum element-wise change in \eqn{V}. The default is 1e-4.
#' @param maxiter Maximum number of iterations, the default is 10000.
#' @param myseed Optional specification of random seed for generating the folds, the default value is 1001.
#' @param prior A logical indicating whether to put larger weights to the groups of larger size; the default value is \code{TRUE}.
#'
#' @return A list of
#'         \item{lambda_seq}{The sequence of tuning parameters used.}
#'        \item{cvm}{The mean cross-validated error rate - a vector of length \code{length(lambda_seq)}}
#'        \item{cvse}{The estimated standard error vector corresponding to \code{cvm}.}
#'        \item{lambda_min}{Value of tuning parameter corresponding to the minimal error in \code{cvm}.}
#'        \item{lambda_1se}{The largest value of tuning parameter such that the correspondig error is within 1 standard error of the minimal error in \code{cvm}.}
#'        \item{nfeature_mat}{A \code{nfolds} x \code{length(lambda_seq)} matrix of the number of selected features.}
#'        \item{error_mat}{A \code{nfolds} x \code{length(lambda_seq)} matrix of the error rates.}
#'
#' @example man/examples/cv_DAP_eg.R
#'
#' @export 
#' 
cv_DAP <- function(X, Y, lambda_seq, nfolds = 5, eps = 1e-4, maxiter = 1000, myseed = 1001, prior = TRUE){
  
  n = length(Y)
  n_lambda = length(lambda_seq)
  
  error_mat = matrix(0.5, nfolds, n_lambda)
  nfeature_mat = matrix(NA, nfolds, n_lambda)
  
  ####random set split the whole data set into k folds corresponding to n1 and n2
  set.seed(myseed)
  id = rep(NA, n)
  id[Y==1]= sample(rep(seq_len(nfolds), length.out = sum(Y == 1)))
  id[Y==2]= sample(rep(seq_len(nfolds), length.out = sum(Y == 2)))
  
  for (nf in 1: nfolds){
    cat(nf)
    ####set training data and test data
    xtrain = X[id != nf, ]
    ytrain = Y[id != nf]
    xtest = X[id == nf, ]
    ytest = Y[id == nf]
    out_s = standardizeData(xtrain, ytrain, center = TRUE)
    
    # Calculate mean difference based on test data
    d_test = colMeans(xtest[ytest == 1, ])-colMeans(xtest[ytest == 2, ])
    
    ####use solve_proj_seq
    fit_tmp = solve_DAP_seq(X1 = out_s$X1, X2 = out_s$X2, lambda_seq = lambda_seq, eps = eps, maxiter = maxiter, feature_max = n)
    nfeature_mat[nf, 1:length(fit_tmp$nfeature_vec)] = fit_tmp$nfeature_vec
    
    #### Calculate errors for each lambda
    for (j in 1:length(fit_tmp$lambda_seq)){
      V = cbind(diag(1 / out_s$coef1) %*% fit_tmp$V1_mat[, j], diag(1 / out_s$coef2) %*% fit_tmp$V2_mat[, j])
      ypred = classify_DAP(xtrain - matrix(out_s$Xmean, nrow(xtrain), ncol(xtest), byrow = TRUE), ytrain, xtest = xtest - matrix(out_s$Xmean, nrow(xtest), ncol(xtest), byrow = TRUE), V, prior = prior)
      error_mat[nf, j] = sum(ypred != ytest) / length(ytest)
    }
  }
  ####calculate cvm, cvse, lambda_min, lambda_1se
  cvm = colMeans(error_mat)
  index = which.min(cvm)
  lambda_min = lambda_seq[index]
  
  cvse = apply(error_mat, 2, stats::sd) / sqrt(nfolds)
  se_lambda = cvse[index]
  up_b = cvm[index] + se_lambda
  lambda_1se = (lambda_seq[cvm <= up_b])[1]
  
  
  return (list(lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse, lambda_seq = lambda_seq, nfeature_mat = nfeature_mat, error_mat = error_mat))
}