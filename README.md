# DAP
The R package **DAP** provides tools for Discriminant analysis via projections.

## Installation

You can install the **DAP** version from
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
The function apply_DAP applies the proposed method DAP

#### RDAP
The function apply_RDAP applies RDA to the subset selected by DAP. 

#### Cross Validation
The function cv_DAP selects tuning parameter using 5-fold cross validation as default.

### Not Exported

*standardizeData
*solve_DAP_R 
*solve_DAP_C
*solve_DAP_seq
*classify_DAP
*update_hdrda
*predict.hdrda
*hdrda_cv_variables
*quadform
*apply_hdrda_full

## License

This package is free and open source software, licensed under GPL (>=2).