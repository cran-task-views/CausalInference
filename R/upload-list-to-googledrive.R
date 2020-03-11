library(tidyverse)
library(readr)

pkg_name_text <- read_lines(here::here("data/causal_packages.txt"))

causal_tv <- tibble(
  pkg_name = pkg_name_text,
  reviewer = "imke",
  url = glue::glue("https://CRAN.R-project.org/package={pkg_name_text}"),
  updated_in_12m = "",
  date_updated = "",
  vignette = "",
  comments = "",
  decision = ""
)

causal_tv

write_csv(causal_tv, path = here::here("data/causal_tv.csv"))

# # upload to google drive
# library(googledrive)
# googlesheets::gs_upload(file = here::here("data/causael_tv.csv"),
#                         "causal-inference-taskview")
