Causal Inference Task View
================
Imke Mayer, Pan Zhao and Julie Josse
07 June, 2022

<!-- README.md is generated from README.Rmd. Please edit that file -->

The Causal Inference task view source is in the file
`CausalInference.ctv` that can be transformed into an HTML file using
the **R** package [ctv](https://CRAN.R-project.org/package=ctv)

``` r
library(ctv)
setwd(paste(getwd(), "source", sep = "/"))
# ctv:::ctv_xml_to_rmd("CausalInference.ctv")
source_ctv <- "CausalInference.md"
ctv2html(source_ctv)
```

In addition, CTV can be checked before submission to CRAN with:

``` r
check_ctv_packages("source/CausalInference.md")
```
