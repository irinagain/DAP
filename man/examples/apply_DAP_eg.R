## This is an example for apply_DAP

## Generate data
n1 = n2 = 50
n_test = 50
p = 100
mu1 = rep(0, p)
mu2 = rep(3, p)
Sigma1 = diag(p)
Sigma2 = 0.5* diag(p)

## Build training data and test data
x1 = MASS::mvrnorm(n = n1, mu = mu1, Sigma = Sigma1)
x2 = MASS::mvrnorm(n = n2, mu = mu2, Sigma = Sigma2)
y1 = rep(1, n1)
y2 = rep(2, n2)
xtrain <- rbind(x1, x2)
x1_test = MASS::mvrnorm(n = n_test, mu = mu1, Sigma = Sigma1)
x2_test = MASS::mvrnorm(n = n_test, mu = mu2, Sigma = Sigma2)
xtest <- rbind(x1_test, x2_test)
ytrain <- c(rep(1, n1), rep(2, n2))
ytest <- c(rep(1, n_test), rep(2, n_test))

## Apply DAP

# Given ytest, the function will return a miclassification error rate.
ClassificationError = apply_DAP(xtrain, ytrain, xtest, ytest)

# Without ytest, the function will return predictions.
Ypredict = apply_DAP(xtrain, ytrain, xtest)
