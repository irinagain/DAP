####This function is used to solve optimization problem with  a sequence of lambda#####
#' Solve Optimization Problem for a sequence of lambda
#'
#' Solving group lasso problem using block-coordinate algorithm for a
#' sequence of lambda.
#'
#' @param X1 A matrix (\code{n_1} by p) of group 1 data (scaled).
#' @param X2 A matrix (\code{n_2} by p) of group 2 data (scaled).
#' @param lambda_seq A sequence of tunning parameter lambda.
#' @param eps Convergence threshold for block-coordinate decent
#' algorithm. Each block-coordinate decent algorithm loop continuoues
#' until the maximum iteration number exceeds \code{maxiter} or the
#' maximum element-wise change in \eqn{V} is less than \code{eps}.
#' Default is 1e-4.
#' @param m_max Maximum number of iterations. Default is 10000.
#' @param feature_max The maximum number of features that can be
#' selected. Default is the total sample size. Once the maximum is
#' reached, the function will return the reasults and larger lambda
#' values won't be applied.
#'
#' @return A list with the following components.
#'        \item{V1_mat}{A matrix of the first projection vector V1
#'        correstponding to the sequence of lambda.}
#'        \item{V2_mat}{A matrix of the second projection vector V2
#'        correstponding to the sequence of lambda.}
#'        \item{lambda_seq}{The sequence of lambda that has been
#'        applied to the data.}
#'        \item{nfeature_vec}{A sequence of number of selected
#'        features.}
#'
#' @example man/examples/solve_DAP_seq_eg.R
#'
#' @export
#'
solve_DAP_seq <- function(X1, X2, lambda_seq, eps = 1e-4, m_max = 10000, feature_max = nrow(X1) + nrow(X2)){
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
