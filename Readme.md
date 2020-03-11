Causal Inference Task View
================
Imke Mayer, Julie Josse and Stefan Wager
11 March, 2020

<!-- README.md is generated from README.Rmd. Please edit that file -->
The Causal Inference task view source is in the file `CausalInference.ctv` that can be transformed into an HTML file using the **R** package [ctv](https://CRAN.R-project.org/package=ctv)

``` r
library(ctv)
setwd(paste(getwd(), "source", sep = "/"))
source_ctv <- read.ctv("CausalInference.ctv")
```

``` r
ctv2html(source_ctv)
```

In addition, CTV can be checked before submission to CRAN with:

``` r
check_ctv_packages("source/CausalInference.ctv")
```
