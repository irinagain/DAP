#' Solves DAP optimization problem for a given sequence of lambda values
#'
#' Uses block-coordinate descent algorithm with warm initializations, starts with the maximal supplied lambda value.
#'
#' @param X1 A n1 x p matrix of group 1 data (scaled).
#' @param X2 A n2 x p matrix of group 2 data (scaled).
#' @param lambda_seq A supplied sequence of tunning parameters.
#' @param eps Convergence threshold for the block-coordinate decent algorithm based on the maximum element-wise change in \eqn{V}. The default is 1e-4.
#' @param maxiter Maximum number of iterations, the default is 10000.
#' @param feature_max The maximum number of features that can be
#' selected. Default is the total sample size. Once the maximum is
#' reached, the function will return the results and larger lambda
#' values won't be applied.
#'
#' @return A list of
#'        \item{V1_mat}{A matrix of the first projection vector V1
#'        corresponding to the sequence of lambda.}
#'        \item{V2_mat}{A matrix of the second projection vector V2
#'        corresponding to the sequence of lambda.}
#'        \item{lambda_seq}{The sequence of lambda that has been
#'        applied to the data.}
#'        \item{nfeature_vec}{A sequence of number of selected
#'        features.}
#'
#' @example man/examples/solve_DAP_seq_eg.R
#'
#' @export
#'
solve_DAP_seq <- function(X1, X2, lambda_seq, eps = 1e-4, maxiter = 10000, feature_max = nrow(X1) + nrow(X2)){
  p =ncol(X1)
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
  return (list(V1_mat = V1_mat, V2_mat = V2_mat, lambda_seq = lambda_seq, nfeature_vec = nfeature_vec))
}
