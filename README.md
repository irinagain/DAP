# DAP
The R package **DAP** provides tools for Discriminant analysis via projections.

## Installation

You can install the **DAP** version from
[Github](https://github.com/irinagain/DAP).
.
```s
install.packages("DAP")
devtools::install_github("irinagain/DAP")
```
## Usage

```s
library(DAP)

# Example 1

# Example 2
```

## List of functions

### Exported

#### DAP
The function apply_DAP applies the proposed method DAP

#### Cross Validation
The function cv_DAP selects tuning parameter using 5-fold cross validation as default.

#### Standardize original data
The function standardizeData can scale the data matrix \eqn{X} by centering each column and scaling the diagonal element to be 1.

#### solve optimization problem with  a single lambda (using C code)
The function solve_DAP_C applies block-coordinate algorithm to solve the optimization probalem with a single lambda.

#### solve optimization problem with  a sequence of lambda
The function solve_DAP_seq solve the optimization with a sequence of lambda using solve_DAP_C.

#### implement new classification rule using 2-dimensional space formed by DAP
The function classify_DAP apply the optimal solution of DAP to perform the classification.

### Not Exported
* update_hdrda
* predict.hdrda
* hdrda_cv_variables
* quadform
* apply_hdrda_full

## License

This package is free and open source software, licensed under GPL (>=2).
