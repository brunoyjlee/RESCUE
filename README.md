# RESCUE: Recovery of idiosyncratic expression patterns

## Overview
We developed a new computational method that can recover the spatial expression patterns missed by standard single-cell inference on spatial transcriptomics (ST) data. RESCUE is based on making a fundamental correction to the class of statistical models currently used for integrating ST data and reference sc/snRNA-seq data. RESCUE partitions observed gene expression into one component explainable by reference sc/snRNA-seq data, which we refer to as **“canonical”** expression, and another component that is not, which we refer to as **“idiosyncratic”** expression. Consequently, our purified canonical ST expression component can achieve more accurate single-cell inference when used with cell-type deconvolution. Our recovered idiosyncratic ST expression component gives us the new ability to specifically target biology that is distinct from our reference data.

<img width="1400" height="749" alt="Image" src="https://github.com/user-attachments/assets/9ba4cbe4-2a69-4ce3-ba37-84655117f623" />

## Dependencies
- R version >= 4.4.0.
- Dependent R packages: RcppML, nnls, parallel, stats, Matrix, tidyverse, dplyr, data.table, ggplot2
```ruby
# install devtools if necessary
install.packages('devtools')

# install the RESCUE package
devtools::install_github("brunoyjlee/RESCUE")

# load package
library(RESCUE)
```

How to use `RESCUE`
-------------------
- Cichlid Fish Visium + snRNA-seq: Details in [Tutorial](https://github.com/brunoyjlee/RESCUE/blob/main/docs/fish-example.md).
