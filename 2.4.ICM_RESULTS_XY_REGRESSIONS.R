### Integrated Community Model (ICM)
### Survival-abundance regression line plots
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

# Years whose survival estimates were modeled using the regression structure
regression_years <- setdiff(
  community_years,
  drop_regression_years
)

# Number of points used to draw each regression curve
curve_length <- 200

################################################################################
#########--------------- Clean posterior samples ---------------################
################################################################################

cleaned <- clean_mcmc_samples(icm_samples)

icm_clean <- cleaned$icm_clean
post_mat <- cleaned$post_mat
bad_summary <- cleaned$bad_summary

# Review posterior columns removed because they contained nonfinite values
bad_summary

################################################################################
#########--------------- Survival model definitions -------------##############
################################################################################

elk_model_specs <- list(
  
  "Calf survival" = list(
    survival_parameter = "elk_s_c",
    intercept = "beta0_calfSurv",
    terms = c(
      wolf_N_tot = "beta1_calfSurv_wolfN",
      griz_N = "beta3_calfSurv_grizN",
      cougar_N = "beta4_calfSurv_cougarN",
      elk_N_female = "beta5_calfSurv_elkN"
    )
  ),
  
  "Young adult survival" = list(
    survival_parameter = "elk_s_ya",
    intercept = "beta0_yaSurv",
    terms = c(
      wolf_N_tot = "beta1_yaSurv_wolfN",
      cougar_N = "beta4_yaSurv_cougarN",
      elk_N_female = "beta5_yaSurv_elkN"
    )
  ),
  
  "Old adult survival" = list(
    survival_parameter = "elk_s_oa",
    intercept = "beta0_oaSurv",
    terms = c(
      wolf_N_tot = "beta1_oaSurv_wolfN",
      cougar_N = "beta4_oaSurv_cougarN",
      elk_N_female = "beta5_oaSurv_elkN"
    )
  )
)


wolf_model_specs <- list(
  
  "Pup survival" = list(
    survival_parameter = "wolf_s_p",
    intercept = "beta0_wpupSurv",
    terms = c(
      elk_N_female = "beta1_wpupSurv_elkN",
      bison_N = "beta2_wpupSurv_bisonN",
      wolf_N_tot = "beta3_wpupSurv_wolfN"
    )
  ),
  
  "Adult survival" = list(
    survival_parameter = "wolf_s_a",
    intercept = "beta0_wadSurv",
    terms = c(
      elk_N_female = "beta1_wadSurv_elkN",
      bison_N = "beta2_wadSurv_bisonN",
      wolf_N_tot = "beta3_wadSurv_wolfN"
    )
  )
)

################################################################################
#########--------------- Abundance predictor definitions --------###############
################################################################################

predictor_specs <- list(
  
  wolf_N_tot = list(
    parameter = "wolf_N_tot",
    mean = icm_constants$wolf_tot_mean,
    sd = icm_constants$wolf_tot_sd,
    label = "Wolf abundance"
  ),
  
  griz_N = list(
    parameter = "griz_N",
    mean = icm_constants$griz_N_mean,
    sd = icm_constants$griz_N_sd,
    label = "Grizzly abundance"
  ),
  
  cougar_N = list(
    parameter = "cougar_N",
    mean = icm_constants$cougar_N_mean,
    sd = icm_constants$cougar_N_sd,
    label = "Cougar abundance"
  ),
  
  elk_N_female = list(
    parameter = "elk_N_female",
    mean = icm_constants$elk_N_female_mean,
    sd = icm_constants$elk_N_female_sd,
    label = "Female elk abundance"
  ),
  
  bison_N = list(
    parameter = "bison_N",
    mean = icm_constants$bison_N_mean,
    sd = icm_constants$bison_N_sd,
    label = "Northern Range bison abundance"
  )
)


################################################################################
#########----------------- Elk line-plot data -------------------###############
################################################################################

elk_wolf_points <- make_survival_abundance_points(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "wolf_N_tot",
  predictor_specs = predictor_specs,
  community_years = community_years,
  regression_years = regression_years
)

elk_wolf_curve <- make_survival_abundance_curve(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "wolf_N_tot",
  predictor_specs = predictor_specs,
  point_df = elk_wolf_points,
  length_out = curve_length
)


elk_grizzly_points <- make_survival_abundance_points(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "griz_N",
  predictor_specs = predictor_specs,
  community_years = community_years,
  regression_years = regression_years
)

elk_grizzly_curve <- make_survival_abundance_curve(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "griz_N",
  predictor_specs = predictor_specs,
  point_df = elk_grizzly_points,
  length_out = curve_length
)


elk_cougar_points <- make_survival_abundance_points(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "cougar_N",
  predictor_specs = predictor_specs,
  community_years = community_years,
  regression_years = regression_years
)

elk_cougar_curve <- make_survival_abundance_curve(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "cougar_N",
  predictor_specs = predictor_specs,
  point_df = elk_cougar_points,
  length_out = curve_length
)


elk_density_points <- make_survival_abundance_points(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "elk_N_female",
  predictor_specs = predictor_specs,
  community_years = community_years,
  regression_years = regression_years
)

elk_density_curve <- make_survival_abundance_curve(
  post_mat = post_mat,
  model_specs = elk_model_specs,
  predictor_name = "elk_N_female",
  predictor_specs = predictor_specs,
  point_df = elk_density_points,
  length_out = curve_length
)

################################################################################
#########---------------- Wolf line-plot data -------------------###############
################################################################################

wolf_elk_points <- make_survival_abundance_points(
  post_mat = post_mat,
  model_specs = wolf_model_specs,
  predictor_name = "elk_N_female",
  predictor_specs = predictor_specs,
  community_years = community_years,
  regression_years = regression_years
)

wolf_elk_curve <- make_survival_abundance_curve(
  post_mat = post_mat,
  model_specs = wolf_model_specs,
  predictor_name = "elk_N_female",
  predictor_specs = predictor_specs,
  point_df = wolf_elk_points,
  length_out = curve_length
)


wolf_bison_points <- make_survival_abundance_points(
  post_mat = post_mat,
  model_specs = wolf_model_specs,
  predictor_name = "bison_N",
  predictor_specs = predictor_specs,
  community_years = community_years,
  regression_years = regression_years
)

wolf_bison_curve <- make_survival_abundance_curve(
  post_mat = post_mat,
  model_specs = wolf_model_specs,
  predictor_name = "bison_N",
  predictor_specs = predictor_specs,
  point_df = wolf_bison_points,
  length_out = curve_length
)


wolf_density_points <- make_survival_abundance_points(
  post_mat = post_mat,
  model_specs = wolf_model_specs,
  predictor_name = "wolf_N_tot",
  predictor_specs = predictor_specs,
  community_years = community_years,
  regression_years = regression_years
)

wolf_density_curve <- make_survival_abundance_curve(
  post_mat = post_mat,
  model_specs = wolf_model_specs,
  predictor_name = "wolf_N_tot",
  predictor_specs = predictor_specs,
  point_df = wolf_density_points,
  length_out = curve_length
)

################################################################################
#########--------------------- Elk plots ------------------------###############
################################################################################

elk_wolf_plot <- plot_survival_abundance(
  point_df = elk_wolf_points,
  curve_df = elk_wolf_curve,
  x_label = predictor_specs$wolf_N_tot$label,
  title = "Estimated relationship between wolf abundance and elk survival"
)

elk_grizzly_plot <- plot_survival_abundance(
  point_df = elk_grizzly_points,
  curve_df = elk_grizzly_curve,
  x_label = predictor_specs$griz_N$label,
  title = "Estimated relationship between grizzly abundance and elk calf survival"
)

elk_cougar_plot <- plot_survival_abundance(
  point_df = elk_cougar_points,
  curve_df = elk_cougar_curve,
  x_label = predictor_specs$cougar_N$label,
  title = "Estimated relationship between cougar abundance and elk survival"
)

elk_density_plot <- plot_survival_abundance(
  point_df = elk_density_points,
  curve_df = elk_density_curve,
  x_label = predictor_specs$elk_N_female$label,
  title = "Estimated relationship between elk abundance and elk survival"
)

################################################################################
#########-------------------- Wolf plots ------------------------###############
################################################################################

wolf_elk_plot <- plot_survival_abundance(
  point_df = wolf_elk_points,
  curve_df = wolf_elk_curve,
  x_label = predictor_specs$elk_N_female$label,
  title = "Estimated relationship between elk abundance and wolf survival"
)

wolf_bison_plot <- plot_survival_abundance(
  point_df = wolf_bison_points,
  curve_df = wolf_bison_curve,
  x_label = predictor_specs$bison_N$label,
  title = "Estimated relationship between bison abundance and wolf survival"
)

wolf_density_plot <- plot_survival_abundance(
  point_df = wolf_density_points,
  curve_df = wolf_density_curve,
  x_label = predictor_specs$wolf_N_tot$label,
  title = "Estimated relationship between wolf abundance and wolf survival"
)

################################################################################
#########--------------------- View plots -----------------------###############
################################################################################

# Elk survival-abundance pairings
elk_wolf_plot
elk_grizzly_plot
elk_cougar_plot
elk_density_plot

# Wolf survival-abundance pairings
wolf_elk_plot
wolf_bison_plot
wolf_density_plot

################################################################################
################################################################################
################################################################################