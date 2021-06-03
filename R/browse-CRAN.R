# packages
library(tidytext)
library(tidyverse)
library(XML)

# keywords for missing data
mda_words <- c("causal",
               "causal inference",
               "treatment effect",
               "treatment effects",
               "treatment",
               "ATE",
               "HTE",
               "policy",
               "policy learning",
               "causal analysis",
               "intervention",
               "do-calculus",
               "back-door criterion",
               "propensity",
               "propensity score",
               "propensity scores",
               "potential outcome",
               "potential outcomes",
               "treatment",
               "RCT",
               "structural equation model",
               "structural equation models",
               "SEM",
               "causal model",
               "causal models",
               "causal structure",
               "instrumental variable",
               "instrumental variables",
               "IV",
               "IPW",
               "IPTW",
               "IPSW")

# extract already listed packages
cur_page <- htmlTreeParse("../source/CausalInference.ctv", useInternalNodes = TRUE)
listed_packages <- unique(xpathSApply(cur_page, "//pkg", xmlValue))

# browse CRAN
cran_db <- tools::CRAN_package_db()

cran_tbl <- tibble::as_tibble(cran_db[,-65])

cran_tbl_short <- cran_tbl %>% select(Package, Description)

tidy_desc <- cran_tbl_short %>% unnest_tokens(word, Description)

data("stop_words")

cleaned_desc <- tidy_desc %>% anti_join(stop_words)

causal_pkgs <- cleaned_desc %>% 
  group_by(Package) %>%
  filter(word %in% mda_words) %>%
  #filter(!(Package %in% listed_packages)) %>%
  ungroup()


n_distinct(causal_pkgs$Package)
#> [1] 482

write.table(unique(causal_pkgs$Package), file = "../data/causal_packages.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

table(causal_pkgs$word)

# causal intervention       policy   propensity    treatment 
#    226           49           53           82          491

sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.15.3
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] XML_3.98-1.20   forcats_0.4.0   stringr_1.4.0   dplyr_0.8.3     purrr_0.3.3     readr_1.3.1    
# [7] tidyr_1.0.0     tibble_2.1.3    ggplot2_3.2.1   tidyverse_1.2.1 tidytext_0.2.1 
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.3        cellranger_1.1.0  pillar_1.4.2      compiler_3.5.1    tokenizers_0.2.1 
# [6] tools_3.5.1       zeallot_0.1.0     lubridate_1.7.4   jsonlite_1.6      lifecycle_0.1.0  
# [11] nlme_3.1-140      gtable_0.3.0      lattice_0.20-38   pkgconfig_2.0.3   rlang_0.4.2      
# [16] Matrix_1.2-17     cli_1.1.0         rstudioapi_0.10   yaml_2.2.0        haven_2.1.0      
# [21] withr_2.1.2       xml2_1.2.0        httr_1.4.0        janeaustenr_0.1.5 hms_0.4.2        
# [26] generics_0.0.2    vctrs_0.2.0       grid_3.5.1        tidyselect_0.2.5  glue_1.3.1       
# [31] R6_2.4.1          readxl_1.3.1      modelr_0.1.4      magrittr_1.5      backports_1.1.5  
# [36] scales_1.1.0      SnowballC_0.6.0   rvest_0.3.4       assertthat_0.2.1  colorspace_1.4-1 
# [41] stringi_1.4.3     lazyeval_0.2.2    munsell_0.5.0     broom_0.5.2       crayon_1.3.4    
