#' Classification via DAP
#'
#' Classify observations in the test set using the supplied matrix V and the training data.
#'
#' @param xtrain A n x p training dataset; n observations on the rows and p features on the columns.
#' @param ytrain A n vector of training group labels, either 1 or 2.
#' @param xtest A m x p testing dataset; m observations on the rows and p features on the columns.
#' @param V A p x 2 projection matrix.
#' @param prior A logical indicating whether to put larger weights to the groups of larger size; the default value is \code{TRUE}.
#'
#' @return Predicted class labels for the test data.
#'
#' @example man/examples/classify_DAP_eg.R
#'
#' @export
#'
classify_DAP <- function(xtrain, ytrain, xtest, V, prior = TRUE){
  Vsvd = base::svd(V)$d
  ## Only for 2 column case for now, and 2 classes
  if ((min(Vsvd) < 1e-6) | (length(Vsvd) == 1)){
    V = V[ , 1, drop = FALSE]
    trainproj = xtrain %*% V
    testproj = xtest %*% V
    Omega1 = 1 / stats::var(trainproj[ytrain == 1, ])
    Omega2 = 1 / stats::var(trainproj[ytrain == 2, ])
    mu1 = mean(trainproj[ytrain == 1, ])
    mu2 = mean(trainproj[ytrain == 2, ])
    predict_vec = rep(2, nrow(xtest))
    difference = diag(tcrossprod(testproj) * (Omega1 - Omega2)) - 2 * as.vector(testproj * (Omega1 * mu1 - Omega2 * mu2)) - mu2^2 * Omega2 + mu1^2 * Omega1 - log(Omega1) + log(Omega2)
    if (prior == TRUE){
      difference = difference - 2*log(sum(ytrain == 1)/sum(ytrain == 2))
    }
    predict_vec[difference < 0] = 1
    return(predict_vec)
  }else{
    # Use qda from MASS for prediction
    if (prior == TRUE){
      out = MASS::qda(xtrain %*% V, grouping = ytrain)
    }else{
      out = MASS::qda(xtrain %*% V, grouping = ytrain, prior = c(1/2,1/2))
    }
    pred = stats::predict(out, newdata = xtest %*% V)
    
    return(as.numeric(pred$class))
  }
}
