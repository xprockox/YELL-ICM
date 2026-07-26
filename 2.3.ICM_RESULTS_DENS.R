### Integrated Community Model (ICM)
### Coefficient density plots
### Last updated: July 26, 2026
### xprockox@gmail.com

################################################################################
#########------------------- Load packages ---------------------################
################################################################################

library(tidyverse)
library(coda)
library(stringr)

################################################################################
############------------------ Load results --------------------################
################################################################################

# load model results
load("data/outputs/ICM_parallel_output_2026-07-24.RData")

################################################################################
#########---------------------- Settings -----------------------################
################################################################################

# Coefficient selectors:
# - "all" = all listed coefficients
# - "intercepts" = only intercepts
# - "non_intercepts" = everything except intercepts
# - character vector = exact pretty labels to retain
# - "regex:<pattern>" = regex matched against pretty labels
#
# Examples:
# dens_coefs_elk <- "all"
# dens_coefs_elk <- "non_intercepts"
# dens_coefs_elk <- "regex:(wolf|grizzly|cougar)"
# dens_coefs_elk <- c("Calf survival wolf effect", "YA survival wolf effect")

dens_coefs_elk <- "non_intercepts"
dens_coefs_ya <- "non_intercepts"
dens_coefs_wolf <- "all"

################################################################################
############--------------- Load helper functions --------------################
################################################################################

# load model results
source("ICM_RESULTS_HELPER_FUNCTIONS.R")

################################################################################
#########---------------- Coefficient labels -------------------################
################################################################################

pretty_coef_labels <- c(
  # elk calf survival
  beta0_calfSurv = "Calf survival intercept",
  beta1_calfSurv_wolfN = "Calf survival wolf effect",
  beta2_calfSurv_winterSeverity = "Calf survival winter severity effect",
  beta3_calfSurv_grizN = "Calf survival grizzly effect",
  beta4_calfSurv_cougarN = "Calf survival cougar effect",
  beta5_calfSurv_elkN = "Calf survival density dependence effect",
  
  # elk young-adult survival
  beta0_yaSurv = "YA survival intercept",
  beta1_yaSurv_wolfN = "YA survival wolf effect",
  beta2_yaSurv_winterSeverity = "YA survival winter severity effect",
  beta4_yaSurv_cougarN = "YA survival cougar effect",
  beta5_yaSurv_elkN = "YA survival density dependence effect",
  beta6_yaSurv_annualNpp = "YA survival annual NPP effect",
  beta7_yaSurv_browndown = "YA survival browndown effect",
  beta8_yaSurv_pdsi = "YA survival PDSI effect",
  
  # elk old-adult survival
  beta0_oaSurv = "OA survival intercept",
  beta1_oaSurv_wolfN = "OA survival wolf effect",
  beta2_oaSurv_winterSeverity = "OA survival winter severity effect",
  beta4_oaSurv_cougarN = "OA survival cougar effect",
  beta5_oaSurv_elkN = "OA survival density dependence effect",
  beta6_oaSurv_annualNpp = "OA survival annual NPP effect",
  beta7_oaSurv_browndown = "OA survival browndown effect",
  beta8_oaSurv_pdsi = "OA survival PDSI effect",
  
  # wolf pup survival
  beta0_wpupSurv = "Wolf pup survival intercept",
  beta1_wpupSurv_elkN = "Wolf pup survival elk effect",
  beta2_wpupSurv_bisonN = "Wolf pup survival bison effect",
  beta3_wpupSurv_wolfN = "Wolf pup survival density dependence effect",
  
  # wolf adult survival
  beta0_wadSurv = "Wolf adult survival intercept",
  beta1_wadSurv_elkN = "Wolf adult survival elk effect",
  beta2_wadSurv_bisonN = "Wolf adult survival bison effect",
  beta3_wadSurv_wolfN = "Wolf adult survival density dependence effect",
  
  # species-specific population-process coefficients
  beta1_griz_elkCalves = "Grizzly elk calves effect",
  beta1_bison_cull = "Bison cull/harvest effect",
  beta1_cougar_elk = "Cougar elk effect"
)

################################################################################
#########------------- Regression coefficient groups ------------###############
################################################################################

elk_model_specs <- list(
  "Calf survival" = list(
    intercept = "beta0_calfSurv",
    terms = c(
      wolf_N_tot = "beta1_calfSurv_wolfN",
      winter_severity = "beta2_calfSurv_winterSeverity",
      griz_N = "beta3_calfSurv_grizN",
      cougar_N = "beta4_calfSurv_cougarN",
      elk_N_female = "beta5_calfSurv_elkN"
    )
  ),
  
  "Young adult survival" = list(
    intercept = "beta0_yaSurv",
    terms = c(
      wolf_N_tot = "beta1_yaSurv_wolfN",
      winter_severity = "beta2_yaSurv_winterSeverity",
      cougar_N = "beta4_yaSurv_cougarN",
      elk_N_female = "beta5_yaSurv_elkN",
      annual_npp = "beta6_yaSurv_annualNpp",
      browndown_onset_greenness_min = "beta7_yaSurv_browndown",
      summer_avg_pdsi = "beta8_yaSurv_pdsi"
    )
  ),
  
  "Old adult survival" = list(
    intercept = "beta0_oaSurv",
    terms = c(
      wolf_N_tot = "beta1_oaSurv_wolfN",
      winter_severity = "beta2_oaSurv_winterSeverity",
      cougar_N = "beta4_oaSurv_cougarN",
      elk_N_female = "beta5_oaSurv_elkN",
      annual_npp = "beta6_oaSurv_annualNpp",
      browndown_onset_greenness_min = "beta7_oaSurv_browndown",
      summer_avg_pdsi = "beta8_oaSurv_pdsi"
    )
  )
)

wolf_model_specs <- list(
  "Wolf pup survival" = list(
    intercept = "beta0_wpupSurv",
    terms = c(
      elk_N_female = "beta1_wpupSurv_elkN",
      NR_Bison = "beta2_wpupSurv_bisonN",
      wolf_N_tot = "beta3_wpupSurv_wolfN"
    )
  ),
  
  "Wolf adult survival" = list(
    intercept = "beta0_wadSurv",
    terms = c(
      elk_N_female = "beta1_wadSurv_elkN",
      NR_Bison = "beta2_wadSurv_bisonN",
      wolf_N_tot = "beta3_wadSurv_wolfN"
    )
  )
)

################################################################################
#########--------------- Clean posterior samples ---------------################
################################################################################

cleaned <- clean_mcmc_samples(icm_samples)

icm_clean <- cleaned$icm_clean
post_mat <- cleaned$post_mat
bad_summary <- cleaned$bad_summary

# Review any posterior columns removed because they contained nonfinite values
bad_summary

################################################################################
#########------------- Select coefficient groups -----------------##############
################################################################################

# All elk-survival regression coefficients listed in elk_model_specs
elk_coef_names <- get_all_model_coef_names(
  model_specs = elk_model_specs,
  post_mat = post_mat,
  include_intercepts = TRUE
)

# All young-adult elk-survival regression coefficients
ya_coef_names <- get_stage_coef_names(
  model_specs = elk_model_specs,
  stage_name = "Young adult survival",
  post_mat = post_mat,
  include_intercept = TRUE
)

# All wolf-survival regression coefficients listed in wolf_model_specs
wolf_coef_names <- get_all_model_coef_names(
  model_specs = wolf_model_specs,
  post_mat = post_mat,
  include_intercepts = TRUE
)

# Inspect selected parameter names
elk_coef_names
ya_coef_names
wolf_coef_names

################################################################################
#########---------------- Elk density plots ---------------------###############
################################################################################

# All elk survival stages
elk_coef_plot <- make_coef_density_plot_stacked(
  post_mat = post_mat,
  coef_names = elk_coef_names,
  label_map = pretty_coef_labels,
  keep = dens_coefs_elk,
  fill_color = "#236192",
  title = "Posterior distributions of elk survival coefficients",
  effect_order = c(
    "Intercept",
    "Wolf abundance",
    "Grizzly abundance",
    "Cougar abundance",
    "Density dependence",
    "Winter severity",
    "Annual NPP",
    "Browndown",
    "PDSI"
  ),
  model_order = c(
    "Calf survival",
    "YA survival",
    "OA survival"
  ),
  scale_factor = 0.35
)

# Young-adult elk survival only
ya_coef_plot <- make_coef_density_plot_stacked(
  post_mat = post_mat,
  coef_names = ya_coef_names,
  label_map = pretty_coef_labels,
  keep = dens_coefs_ya,
  fill_color = "#236192",
  title = "Posterior distributions of young-adult elk survival coefficients",
  effect_order = c(
    "Intercept",
    "Wolf abundance",
    "Cougar abundance",
    "Density dependence",
    "Winter severity",
    "Annual NPP",
    "Browndown",
    "PDSI"
  ),
  model_order = "YA survival",
  scale_factor = 0.35
)

################################################################################
#########---------------- Wolf density plots --------------------###############
################################################################################

wolf_coef_plot <- make_coef_density_plot_stacked(
  post_mat = post_mat,
  coef_names = wolf_coef_names,
  label_map = pretty_coef_labels,
  keep = dens_coefs_wolf,
  fill_color = "#6F263D",
  title = "Posterior distributions of wolf survival coefficients",
  effect_order = c(
    "Intercept",
    "Elk abundance",
    "Bison abundance",
    "Density dependence"
  ),
  model_order = c(
    "Wolf pup survival",
    "Wolf adult survival"
  ),
  scale_factor = 0.35
)

################################################################################
#########----------- Population-process density plots ------------##############
################################################################################

griz_coef_plot <- make_single_coef_density_plot(
  post_mat = post_mat,
  parameter = "beta1_griz_elkCalves",
  fill_color = "#4B7F52",
  title = "Posterior distribution of elk calf effect on grizzly abundance"
)

bison_coef_plot <- make_single_coef_density_plot(
  post_mat = post_mat,
  parameter = "beta1_bison_cull",
  fill_color = "#8C510A",
  title = "Posterior distribution of bison cull/harvest effect on bison abundance"
)

cougar_coef_plot <- make_single_coef_density_plot(
  post_mat = post_mat,
  parameter = "beta1_cougar_elk",
  fill_color = "#7A5C99",
  title = "Posterior distribution of elk effect on cougar abundance"
)

################################################################################
#########--------------------- View plots -----------------------###############
################################################################################

elk_coef_plot
ya_coef_plot
wolf_coef_plot

griz_coef_plot
bison_coef_plot
cougar_coef_plot

################################################################################
################################################################################
################################################################################