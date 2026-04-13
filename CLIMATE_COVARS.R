### Climate covariate extraction from PRISM data 
### 25m spatial resolution, monthly temporal resolution
### Last updated: Apr. 10, 2026
### xprockox@gmail.com

################################################################################
############---------------- Load packages ----------------#####################
################################################################################

library(prism)
library(terra)
library(sf)
library(lubridate)
library(dplyr)
library(stringr)

################################################################################
############---------------- User settings ----------------#####################
################################################################################

prism_set_dl_dir("data/prism_data")

start_year <- 1994
end_year <- 2023

nr_shp <- "data/Northern_Range_Bound/Northern_Range_Bound.shp"

################################################################################
############---------------- Load study area ----------------###################
################################################################################

nr <- vect(nr_shp)

################################################################################
############---------------- Find local PRISM archives ---------------##########
################################################################################

all_archives <- prism_archive_ls()

ppt_archives <- all_archives[
  str_detect(all_archives, "prism_ppt_us_25m_")
]

tmean_archives <- all_archives[
  str_detect(all_archives, "prism_tmean_us_25m_")
]

# keep only desired years
ppt_archives <- ppt_archives[
  str_detect(ppt_archives, paste0("(", paste(start_year:end_year, collapse = "|"), ")"))
]

tmean_archives <- tmean_archives[
  str_detect(tmean_archives, paste0("(", paste(start_year:end_year, collapse = "|"), ")"))
]

################################################################################
############---------------- Build raster stacks ----------------###############
################################################################################

ppt_stack_r <- prism_stack(ppt_archives)
tmean_stack_r <- prism_stack(tmean_archives)

ppt_stack <- rast(ppt_stack_r)
tmean_stack <- rast(tmean_stack_r)   # convert PRISM tmean to deg C

################################################################################
############---------------- Align projection + crop ----------------###########
################################################################################

nr <- project(nr, crs(ppt_stack))

ppt_nr <- crop(ppt_stack, nr) %>% 
  mask(nr)

tmean_nr <- crop(tmean_stack, nr) %>% 
  mask(nr)

################################################################################
############---------------- Extract monthly area means ---------------#########
################################################################################

ppt_monthly <- global(ppt_nr, mean, na.rm = TRUE)[, 1]
tmean_monthly <- global(tmean_nr, mean, na.rm = TRUE)[, 1]

################################################################################
############---------------- Build monthly dataframe ---------------############
################################################################################

dates <- seq(
  as.Date(sprintf("%d-01-01", start_year)),
  as.Date(sprintf("%d-12-01", end_year)),
  by = "month"
)

climate_monthly <- data.frame(
  date = dates,
  year = year(dates),
  month = month(dates),
  ppt_mm = as.numeric(ppt_monthly),
  tmean_c = as.numeric(tmean_monthly)
) %>%
  arrange(date)

print(climate_monthly)

################################################################################
############---------------- Build simple annual summaries ---------------######
################################################################################

################################################################################
############---------------- Build simple annual summaries ---------------######
################################################################################

climate_annual <- climate_monthly %>%
  mutate(
    winter_year = ifelse(month %in% c(12), year + 1, year)
  ) %>%
  group_by(year) %>%
  summarise(
    annual_ppt_mm = sum(ppt_mm, na.rm = TRUE),
    summer_ppt_mm = sum(ppt_mm[month %in% c(6, 7, 8)], na.rm = TRUE),
    summer_tmean_c = mean(tmean_c[month %in% c(6, 7, 8)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    climate_monthly %>%
      mutate(
        winter_year = ifelse(month %in% c(12), year + 1, year)
      ) %>%
      filter(month %in% c(12, 1, 2)) %>%
      group_by(winter_year) %>%
      summarise(
        winter_ppt_mm = sum(ppt_mm, na.rm = TRUE),
        winter_tmean_c = mean(tmean_c, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(year = winter_year),
    by = "year"
  ) %>%
  arrange(year)

print(climate_annual)

################################################################################
############---------------- Export outputs ----------------####################
################################################################################

stop("WARNING: The following lines will overwrite data. Are you sure you would like to proceed?")

write.csv(
  climate_monthly,
  "data/covariates/prism_monthly_precip_tmean.csv",
  row.names = FALSE
)

write.csv(
  climate_annual,
  "data/covariates/prism_annual_precip_tmean.csv",
  row.names = FALSE
)
