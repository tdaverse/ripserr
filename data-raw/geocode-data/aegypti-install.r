library(readr)
library(magrittr)
library(dplyr)
library(stringr)

read_csv("ignore/aegypti_annotated.csv") %>%
  rename_all(tolower) %>%
  rename_all(~ str_remove(., "^geocoded_")) %>%
  filter(year == 2013) %>%
  mutate(polygon_admin = ifelse(
    polygon_admin == -999,
    NA_real_, polygon_admin
  )) %>%
  mutate(country_code = toupper(country_code)) %>%
  select(-year, -country, -country_id, -gaul_ad0,
         -country_name, -country_code) ->
  aegypti

# convert state names to ASCII
aegypti %>%
  mutate(state_name = iconv(state_name, from = "", to = "ASCII//TRANSLIT")) %>%
  mutate(state_name = str_remove_all(state_name, "'|~|\\^")) ->
  aegypti

usethis::use_data(aegypti)
