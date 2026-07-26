### Elk harvest data management
### Last updated: June 3, 2026
### xprockox@gmail.com

###########################################################################################
# packages

library(tidyverse)

###########################################################################################
# data import

dat <- read.csv('data/ElkHarvest.csv')

###########################################################################################
# data management

annual_harvest <- dat %>%
  pivot_longer(
    cols = starts_with("Y"),
    names_to = "year",
    values_to = "n_harvested"
  ) %>%
  mutate(
    year = as.integer(sub("Y", "", year)),
    age_class = case_when(
      Age == 1 ~ "age_1",
      Age >= 2 & Age <= 13 ~ "age_2_13",
      Age >= 14 ~ "age_14_plus"
    )
  ) %>%
  group_by(year, age_class) %>%
  summarise(n = sum(n_harvested, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = age_class,
    values_from = n
  ) %>%
  mutate(
    total_female_harvested = age_1 + age_2_13 + age_14_plus
  ) %>%
  select(year, total_female_harvested, age_1, age_2_13, age_14_plus)

annual_harvest

###########################################################################################
# data writing

stop('The following will overwrite data.\n
     Are you sure you would like to proceed?')

write.csv(annual_harvest, 'data/covariates/annual_elk_harvest.csv')

###########################################################################################