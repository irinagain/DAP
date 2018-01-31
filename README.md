# DAP: Discriminant Analysis via Projections
The R package `DAP` provides tools for high-dimensional binary classification in the case of unequal covariance matrices. It  implements methods from the following paper:
* [Sparse quadratic classification rules via linear dimension reduction](https://arxiv.org/abs/1711.04817) by Gaynanova and Wang (2017).

## Installation

To install the latest version from Github, use
```s
library(devtools)
devtools::install_github("irinagain/DAP")
```
## Usage

```s
library(DAP)
library(MASS)

# Example 

## Specify model parameters
p = 100
mu1 = rep(0, p)
mu2 = c(rep(3, 10), rep(0, p-10))
Sigma1 = diag(p)
Sigma2 = 0.5*diag(p)

## Build training data and test data
n_train = 50
n_test = 50
x1 = MASS::mvrnorm(n = n_train, mu = mu1, Sigma = Sigma1)
x2 = MASS::mvrnorm(n = n_train, mu = mu2, Sigma = Sigma2)
xtrain = rbind(x1, x2)
x1_test = MASS::mvrnorm(n = n_test, mu = mu1, Sigma = Sigma1)
x2_test = MASS::mvrnorm(n = n_test, mu = mu2, Sigma = Sigma2)
xtest = rbind(x1_test, x2_test)
ytrain = c(rep(1, n_train), rep(2, n_train))
ytest = c(rep(1, n_test), rep(2, n_test))

## Apply DAP
# Given ytest, the function returns the miclassification error rate.
ClassificationError = apply_DAP(xtrain, ytrain, xtest, ytest)

# Without ytest, the function returns predicted labels.
Ypredict = apply_DAP(xtrain, ytrain, xtest)

```

<!--## List of functions

### Exported

#### DAP
The function `apply_DAP` applies the proposed method DAP

#### Cross Validation
The function `cv_DAP` selects tuning parameter using 5-fold cross validation as default.

#### Standardize original data
The function `standardizeData` scales the data matrix by centering and scaling each variable.

#### solve optimization problem with  a single lambda (using C code)
The function `solve_DAP_C` applies block-coordinate algorithm to solve the optimization problem with a single lambda.

#### solve optimization problem with  a sequence of lambda
The function `solve_DAP_seq` solves the optimization problem for a sequence of lambda values via `solve_DAP_C`.

#### implement new classification rule using 2-dimensional space formed by DAP
The function `classify_DAP` performs the classification based on found solution.

-->

## License
This package is free and open source software, licensed under GPL (>=2).
