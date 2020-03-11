# download current gdrive

library(googledrive)

# googledrive::drive_download(file = "causal-inference-taskview",
#                             path = "data/causal-inference-taskview.csv")


library(tidyverse)
causal_tv <- read_csv(here::here("data", "causal-inference-taskview.csv"))

causal_tv %>% 
  filter(reviewer == "imke",
         decision == TRUE,
         updated_in_12m == TRUE,
         vignette == TRUE) %>%
  select(pkg_name,
         comments,
         everything()) %>%
  View()
