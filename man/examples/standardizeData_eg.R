# An example for the function standardizeData
n <- 200
p <- 100
mu = rep(3, p)
Sigma = diag(p)
X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
Y <- c(rep(1, n/2), rep(2, n/2))
outcome = standardizeData(X = X, Y = Y, center = TRUE)
