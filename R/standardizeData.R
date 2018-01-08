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