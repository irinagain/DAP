# An example for the function standardizeData

## Generate data
n1 = n2 = 50
n_test = 50
p = 100
mu1 = rep(0, p)
mu2 = rep(3, p)
Sigma1 = diag(p)
Sigma2 = 0.5* diag(p)

## Build training data
x1 = MASS::mvrnorm(n = n1, mu = mu1, Sigma = Sigma1)
x2 = MASS::mvrnorm(n = n2, mu = mu2, Sigma = Sigma2)
y1 = rep(1, n1)
y2 = rep(2, n2)
xtrain = rbind(x1, x2)
ytrain = c(rep(1, n1), rep(2, n2))

## Standardize data
out_s = standardizeData(xtrain, ytrain, center = F)

