### Integrated Community Model (ICM)
### Abundance validation plots
### Last updated: July 26, 2026
### xprockox@gmail.com

################################################################################
#########------------------- Load packages ---------------------################
################################################################################

library(tidyverse)
library(MCMCvis)
library(coda)
library(stringr)

################################################################################
############------------------ Load results --------------------################
################################################################################

# load model results
load("data/outputs/ICM_parallel_output_2026-07-24.RData")

################################################################################
############--------------- Load helper functions --------------################
################################################################################

source("ICM_RESULTS_HELPER_FUNCTIONS.R")

################################################################################
#########----------------------- Settings ----------------------################
################################################################################

# Set FALSE to omit a validation plot
plot_elk <- TRUE
plot_wolves <- TRUE
plot_grizzlies <- TRUE
plot_bison <- TRUE
plot_cougars <- TRUE

################################################################################
#########--------------------- Plot labels ---------------------################
################################################################################

elk_stage_map <- c(
  elk_N_1y = "Yearling",
  elk_N_ya = "Young Adult",
  elk_N_oa = "Old Adult",
  elk_N_female = "Total Females"
)

wolf_stage_map <- c(
  wolf_N_p = "Pups",
  wolf_N_a = "Adults",
  wolf_N_tot = "Total Wolves"
)

################################################################################
#########--------------- Clean posterior samples ---------------################
################################################################################

cleaned <- clean_mcmc_samples(icm_samples)

icm_clean <- cleaned$icm_clean
post_mat <- cleaned$post_mat
bad_summary <- cleaned$bad_summary

# Review parameters removed because they contained nonfinite values
bad_summary

################################################################################
#########---------------- Posterior summaries ------------------################
################################################################################

elk_N_summ <- extract_abundance_summary(
  mcmc_obj = icm_clean,
  params = c(
    "elk_N_1y",
    "elk_N_ya",
    "elk_N_oa",
    "elk_N_female"
  ),
  years = community_years,
  stage_map = elk_stage_map
)

wolf_N_summ <- extract_abundance_summary(
  mcmc_obj = icm_clean,
  params = c(
    "wolf_N_p",
    "wolf_N_a",
    "wolf_N_tot"
  ),
  years = community_years,
  stage_map = wolf_stage_map
)

griz_N_summ <- extract_abundance_summary(
  mcmc_obj = icm_clean,
  params = "griz_N",
  years = community_years,
  stage_map = c(
    griz_N = "Grizzly abundance"
  )
)

bison_N_summ <- extract_abundance_summary(
  mcmc_obj = icm_clean,
  params = "bison_N",
  years = community_years,
  stage_map = c(
    bison_N = "Northern Range bison abundance"
  )
)

cougar_N_summ <- extract_abundance_summary(
  mcmc_obj = icm_clean,
  params = "cougar_N",
  years = community_years,
  stage_map = c(
    cougar_N = "Cougar abundance"
  )
)

################################################################################
#########------------------ Observed datasets ------------------################
################################################################################

# Ensure total female elk abundance is available
if (!"n_female" %in% names(elk_dat_n)) {
  elk_dat_n <- elk_dat_n %>%
    mutate(
      n_female = n_cow + n_calf / 2
    )
}

elk_obs_long <- make_observed_long(
  data = elk_dat_n,
  year_col = "year",
  cols_map = c(
    n_female = "Total Females"
  )
)

wolf_obs_long <- make_observed_long(
  data = wolf_pop,
  year_col = "seasonal.year",
  cols_map = c(
    dec_pups = "Pups",
    dec_adults = "Adults",
    total_abundance = "Total Wolves"
  )
)

griz_obs_long <- make_observed_long(
  data = grizzly,
  year_col = "year",
  cols_map = c(
    griz_N = "Grizzly abundance"
  )
)

bison_obs_long <- make_observed_long(
  data = bison,
  year_col = "year",
  cols_map = c(
    NR_Bison = "Northern Range bison abundance"
  )
)

cougar_obs_long <- make_observed_long(
  data = cougars,
  year_col = "year",
  cols_map = c(
    cougar_N = "Cougar abundance"
  )
)

################################################################################
#########----------------- Validation plots --------------------################
################################################################################

if (plot_elk) {
  
  elk_validation_plot <- plot_validation(
    summary_df = elk_N_summ,
    observed_df = elk_obs_long,
    title = "Elk posterior abundance estimates with observed data"
  )
  
} else {
  
  elk_validation_plot <- NULL
}


if (plot_wolves) {
  
  wolf_validation_plot <- plot_validation(
    summary_df = wolf_N_summ,
    observed_df = wolf_obs_long,
    title = "Wolf posterior abundance estimates with observed data"
  )
  
} else {
  
  wolf_validation_plot <- NULL
}


if (plot_grizzlies) {
  
  griz_validation_plot <- plot_validation(
    summary_df = griz_N_summ,
    observed_df = griz_obs_long,
    title = "Grizzly posterior abundance estimates with observed data"
  )
  
} else {
  
  griz_validation_plot <- NULL
}


if (plot_bison) {
  
  bison_validation_plot <- plot_validation(
    summary_df = bison_N_summ,
    observed_df = bison_obs_long,
    title = "Bison posterior abundance estimates with observed data"
  )
  
} else {
  
  bison_validation_plot <- NULL
}


if (plot_cougars) {
  
  cougar_validation_plot <- plot_validation(
    summary_df = cougar_N_summ,
    observed_df = cougar_obs_long,
    title = "Cougar posterior abundance estimates with observed data"
  )
  
} else {
  
  cougar_validation_plot <- NULL
}

################################################################################
#########--------------------- View plots -----------------------###############
################################################################################

elk_validation_plot
wolf_validation_plot
griz_validation_plot
bison_validation_plot
cougar_validation_plot

################################################################################
################################################################################
################################################################################