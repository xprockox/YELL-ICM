### Grizzly Bear Data Management Script
### Last updated: Apr. 24, 2026
### michael.procko@usu.edu
###
### --- Description: 
###
### The grizzly bear data being used as of Apr. 24, 2026 were scraped from the
### grizzly abundance plot in Gould et al. 2024 (Grizzly IPM paper for GYE) using 
### this Web Plot Digitizer tool: https://automeris.io/
###
### Accordingly, the data contain errant decimals and multiple rows per year. 
### This script collapses that data by using the abundance estimate linked to the 
### row where the year value is closest to a whole number. It also rounds the
### abundance estimate to a whole number
##
#

##################################################################################
# packages
library(tidyverse)

##################################################################################
# data loading and basic cleaning
griz <- read.csv('data/covariates/grizzly_abundances.csv')

colnames(griz) <- c('N','year')

##################################################################################
# data management
griz_yearly <- griz %>%
  mutate(
    year_round = round(year),
    dist_to_whole = abs(year - year_round)
  ) %>%
  group_by(year_round) %>%
  slice_min(order_by = dist_to_whole, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    year = year_round,
    N
  ) %>%
  arrange(year) %>%
  mutate(N = round(N))

griz_yearly

##################################################################################
# data writing 
stop("Warning: the following line with overwrite data. Are you sure you would like to proceed?")
write.csv(griz_yearly, 'grizzly_abundances_cleaned.csv')
