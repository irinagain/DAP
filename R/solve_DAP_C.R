#' Solve Optimization Problem (C version, fixed lambda)
#'
#' Solving group lasso using block-coordinate descent algorithm for a
#' fixed value of lambda. C code is used for coding basic functions
#' to speed up.
#'
#' @useDynLib DAP solveProj_withS
#' 
#' @param X1 A \code{n_1} x \code{p} training data set for group 1 (scaled).
#' @param X2 A \code{n_2} x \code{p} training data set for group 2 (scaled).
#' @param lambda A fixed value of the tuning parameter.
#' @param Vinit Starting point, the default is "NULL".
#' @param eps Convergence threshold for block-coordinate descent
#' algorithm. Each block-coordinate descent algorithm loop continuously
#' until the maximum iteration number exceeds \code{maxiter} or the
#' maximum element-wise change in \code{V} is less than \code{eps}.
#' Default is 1e-4.
#' @param maxiter Maximum number of iterations, the default is 10000.
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
#' @example man/examples/solve_DAP_C_eg.R
#'
#' @export
solve_DAP_C <-function(X1, X2, lambda, Vinit = NULL, eps = 1e-4, maxiter = 10000){
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