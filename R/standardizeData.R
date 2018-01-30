#' Divides the features matrix into two standardized submatrices
#'
#' Given matrix \code{X} with corresponding class labels in \code{Y}, the function column-centers \code{X}, divides it into two submatrices corresponding to each class, and scales the columns of each submatrix to have eucledean norm equal to one.
#'
#' @param X A n x p training dataset; n observations on the rows and p features on the columns.
#' @param Y A n vector of training group labels, either 1 or 2.
#' @param center A logical indicating whether \code{X} should be centered, the default is TRUE.
#'
#' @return A list of
#'     \item{X1}{A n1 x p standardized matrix with observations from group 1.}
#'     \item{X2}{A n2 x p standardized matrix with observations from group 2.}
#'     \item{coef1}{Back-scaling coefficients for \code{X1}.}
#'     \item{coef2}{Back-scaling coefficients for \code{X2}.}
#'     \item{Xmean}{Column means of the matrix \code{X} before centering.}
#'
#' @example man/examples/standardizeData_eg.R
#'
#' @export
standardizeData <- function(X, Y, center = TRUE){
  # center X
  Xmean = colMeans(X)
  if (center){
    X = X - matrix(Xmean, nrow(X), ncol(X), byrow = TRUE)
  }
  X1 = X[Y==1,]
  X2 = X[Y==2,]
  coef1 = sqrt(colSums(X1^2)/sum(Y==1))
  coef2 = sqrt(colSums(X2^2)/sum(Y==2))
  X1 = X1/ matrix(coef1, nrow(X1), ncol(X1), byrow = TRUE)
  X2 = X2/ matrix(coef2, nrow(X2), ncol(X2), byrow = TRUE)
  return(list(X1 = X1, X2 = X2, coef1 = coef1, coef2 = coef2, Xmean = Xmean))
}