#' Solves DAP optimization problem for one given lambda value 
#'
#' Uses block-coordinate descent algorithm to solve DAP problem.
#'
#' @useDynLib DAP solveProj_withS
#' 
#' @param X1 A n1 x p matrix of group 1 data
#' (scaled).
#' @param X2 A n2 x p matrix of group 2 data
#' (scaled).
#' @param lambda A value of the tuning parameter lambda.
#' @param Vinit Optional starting point, the default is NULL, and the algorithm starts with the matrix of zeros.
#' @param eps Convergence threshold for the block-coordinate decent algorithm based on the maximum element-wise change in \eqn{V}. The default is 1e-4.
#' @param maxiter Maximum number of iterations, the default is 10000.
#'
#' @return A list of
#'        \item{V}{The projection matrix.}
#'        \item{nfeature}{The number of selected features.}
#'        \item{iter}{Number of iterations until convergence.}
#'
#' @section Warnings:
#' Please use scaled \code{X1} and \code{X2} for this function, they can be obtained using \code{standardizeData} to do so.
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