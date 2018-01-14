# DAP
The R package `DAP` provides tools for Discriminant analysis via projections, which is based on the following paper:
* [Sparse quadratic classification rules via linear dimension reduction](https://arxiv.org/abs/1711.04817) by Gaynanova and Wang (2017).

## Installation

You can install `DAP` from
[Github](https://github.com/irinagain/DAP).

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

### Not Exported
* update_hdrda
* predict.hdrda
* hdrda_cv_variables
* quadform
* apply_hdrda_full

## License

This package is free and open source software, licensed under GPL (>=2).
