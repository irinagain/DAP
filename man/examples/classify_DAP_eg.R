## This is an example for classify_DAP

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

# Standardize the data
out_s = standardizeData(xtrain, ytrain, center = F)

## Find V
out.proj = solve_DAP_C(X1 = out_s$X1, X2 = out_s$X2, lambda = 0.3)
V = cbind(diag(1/out_s$coef1)%*%out.proj$V[,1],diag(1/out_s$coef2)%*% out.proj$V[,2])

# Predict y using classify_DAP
ypred = classify_DAP(xtrain, ytrain, xtest, V = V)

