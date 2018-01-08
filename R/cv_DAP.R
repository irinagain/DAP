#' Cross-validation for DAP
#'
#' Does k-fold cross-validation for DAP.
#'
#' @param X Training data set. No need for standardization.
#' @param Y Training labels, either "1" or "2".
#' @param lambda_seq A sequence of tunning parameter, lambda.
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
#' @return A list with following component:
#'        \item{lambda_min}{Value of \code{lambda} corresponding to
#'        the minimum \code{cvm}.}
#'        \item{lambda_1se}{The largest value of \code{lambda} such
#'        that the error is within 1 standard error of the minimum
#'        \code{cvm}.}
#'        \item{cvm}{The mean of k-fold cross-validation error
#'        with respect to the sequence of \code{lambda}.}
#'        \item{cvse}{The estimate of standard error of \code{cvm}.}
#'        \item{lambda_seq}{The sequence of \code{lambda} used in the
#'        fits.}
#'        \item{nfeature_mat}{The matrix of the number of selected
#'         features with respect to \code{lambda} sequence.}
#'        \item{error_mat}{The matrix of the fitting errors
#'         with respect to \code{lambda} sequence.}
#'
#' @example man/examples/cv_DAP_eg.R
#'
#' @export 
#' 
cv_DAP <-function(X, Y, lambda_seq, nfolds = 5, eps = 1e-4, m_max = 1000, myseed = 1001, prior = TRUE){
  
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