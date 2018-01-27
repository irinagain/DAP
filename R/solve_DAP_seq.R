#' Solves DAP optimization problem for a given sequence of lambda values
#'
#' Uses block-coordinate descent algorithm with warm initializations, starts with the maximal supplied lambda value.
#'
#' @param X1 A n1 x p matrix of group 1 data (scaled).
#' @param X2 A n2 x p matrix of group 2 data (scaled).
#' @param lambda_seq A supplied sequence of tunning parameters.
#' @param eps Convergence threshold for the block-coordinate decent algorithm based on the maximum element-wise change in \eqn{V}. The default is 1e-4.
#' @param maxiter Maximum number of iterations, the default is 10000.
#' @param feature_max An upper bound on the number of nonzero features in the solution; the default value is the total sample size. The algorithm trims the supplied \code{lambda_seq} to eliminate solutions that exceed \code{feature_max}.
#'
#' @return A list of
#'    \item{lambda_seq}{A sequence of considered lambda values.}
#'        \item{V1_mat}{A p x m matrix with columns corresponding to the 1st projection vector V1 found at each lambda from \code{lambda_seq}.}
#'        \item{V2_mat}{A p x m matrix with columns corresponding to the 2nd projection vector V2 found at each lambda from \code{lambda_seq}.}
#'       
#'        \item{nfeature_vec}{A sequence of corresponding number of selected features for each value in \code{lambda_seq}.}
#'
#' @example man/examples/solve_DAP_seq_eg.R
#'
#' @export
#'
solve_DAP_seq <- function(X1, X2, lambda_seq, eps = 1e-4, maxiter = 10000, feature_max = nrow(X1) + nrow(X2)){
  p = ncol(X1)
  n_lambda = length(lambda_seq)
  
  ####initilize V1_mat, V2_mat, both p by n_lambda, V0
  V1_mat = V2_mat = matrix(NA, nrow = p, ncol = n_lambda)
  V0 = matrix(0, nrow = p, ncol = 2)
  nfeature_vec = rep(NA, n_lambda)
  for (i in 1:n_lambda){
    ####use solve_proj for each lambda
    out = solve_DAP_C(X1, X2, lambda = lambda_seq[i], Vinit = V0, eps = eps, maxiter = maxiter)
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
  return(list(V1_mat = V1_mat, V2_mat = V2_mat, lambda_seq = lambda_seq, nfeature_vec = nfeature_vec))
}
