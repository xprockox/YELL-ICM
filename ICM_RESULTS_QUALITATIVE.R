### Integrated Community Model (ICM)
### Qualitative results summaries
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
#########--------------- Clean posterior samples ---------------################
################################################################################

cleaned <- clean_mcmc_samples(icm_samples)

icm_clean <- cleaned$icm_clean
post_mat <- cleaned$post_mat
bad_summary <- cleaned$bad_summary

# Review posterior columns removed because they contained nonfinite values
bad_summary

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
  
  # population-process coefficients
  beta1_griz_elkCalves = "Grizzly elk calves effect",
  beta1_bison_cull = "Bison cull/harvest effect",
  beta1_cougar_elk = "Cougar elk effect"
)

################################################################################
#########---------------- Evidence classification --------------################
################################################################################

coef_names <- names(pretty_coef_labels)[
  names(pretty_coef_labels) %in% colnames(post_mat)
]

coef_evidence <- bind_rows(
  lapply(
    coef_names,
    function(parameter) {
      
      classify_evidence(
        post_mat[, parameter]
      ) %>%
        mutate(
          parameter = parameter,
          label = pretty_coef_labels[[parameter]],
          .before = 1
        )
    }
  )
) %>%
  mutate(
    evidence = factor(
      evidence,
      levels = c(
        "Strong",
        "Moderate",
        "Weak",
        "Little/none"
      )
    )
  ) %>%
  arrange(
    evidence,
    desc(abs(median))
  )

strong_effects <- coef_evidence %>%
  filter(evidence == "Strong")

moderate_effects <- coef_evidence %>%
  filter(evidence == "Moderate")

weak_effects <- coef_evidence %>%
  filter(evidence == "Weak")

little_effects <- coef_evidence %>%
  filter(evidence == "Little/none")

# Inspect complete coefficient evidence table
coef_evidence

################################################################################
#########---------------- Coefficient summaries ----------------################
################################################################################

cat(
  "\n",
  "--------------------------------------------------------\n",
  "COEFFICIENT EVIDENCE SUMMARY\n",
  "--------------------------------------------------------\n\n",
  
  "-----------------------------------\n",
  "Strong evidence\n",
  "95% credible interval excludes zero\n",
  "-----------------------------------\n",
  format_evidence_lines(
    strong_effects,
    interval = "95"
  ),
  "\n\n",
  
  "-----------------------------------\n",
  "Moderate evidence\n",
  "80% credible interval excludes zero\n",
  "-----------------------------------\n",
  format_evidence_lines(
    moderate_effects,
    interval = "80"
  ),
  "\n\n",
  
  "-----------------------------------\n",
  "Weak evidence\n",
  "50% credible interval excludes zero\n",
  "-----------------------------------\n",
  format_evidence_lines(
    weak_effects,
    interval = "50"
  ),
  "\n\n",
  
  "-----------------------------------\n",
  "Little or no evidence\n",
  "50% credible interval includes zero\n",
  "-----------------------------------\n",
  format_evidence_lines(
    little_effects,
    interval = "50"
  ),
  "\n"
)

################################################################################
#########---------------- Abundance summaries ------------------################
################################################################################

elk_stage_map <- c(
  elk_N_1y = "Yearling",
  elk_N_ya = "Young adult",
  elk_N_oa = "Old adult",
  elk_N_female = "Total females"
)

wolf_stage_map <- c(
  wolf_N_p = "Pups",
  wolf_N_a = "Adults",
  wolf_N_tot = "Total wolves"
)

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

################################################################################
#########---------------- Elk population change ----------------################
################################################################################

elk_total_df <- elk_N_summ %>%
  filter(stage == "Total females") %>%
  arrange(year) %>%
  drop_na(mean)

elk_yearling_df <- elk_N_summ %>%
  filter(stage == "Yearling") %>%
  arrange(year) %>%
  drop_na(mean)

elk_ya_df <- elk_N_summ %>%
  filter(stage == "Young adult") %>%
  arrange(year) %>%
  drop_na(mean)

elk_oa_df <- elk_N_summ %>%
  filter(stage == "Old adult") %>%
  arrange(year) %>%
  drop_na(mean)

elk_first_year <- min(
  elk_total_df$year,
  na.rm = TRUE
)

elk_last_year <- max(
  elk_total_df$year,
  na.rm = TRUE
)

elk_first_total <- elk_total_df %>%
  filter(year == elk_first_year) %>%
  pull(mean)

elk_last_total <- elk_total_df %>%
  filter(year == elk_last_year) %>%
  pull(mean)

elk_pct_change <- 100 * (
  elk_last_total - elk_first_total
) / elk_first_total

elk_lambda_geom <- (
  elk_last_total / elk_first_total
)^(1 / (elk_last_year - elk_first_year))

elk_annual_pct_change <- 100 * (
  elk_lambda_geom - 1
)

# Annual posterior-mean stage distribution
elk_stage_distribution <- elk_yearling_df %>%
  select(
    year,
    yearling = mean
  ) %>%
  left_join(
    elk_ya_df %>%
      select(
        year,
        young_adult = mean
      ),
    by = "year"
  ) %>%
  left_join(
    elk_oa_df %>%
      select(
        year,
        old_adult = mean
      ),
    by = "year"
  ) %>%
  mutate(
    stage_total =
      yearling +
      young_adult +
      old_adult,
    
    prop_yearling =
      yearling / stage_total,
    
    prop_young_adult =
      young_adult / stage_total,
    
    prop_old_adult =
      old_adult / stage_total
  )

elk_first_stage_distribution <- elk_stage_distribution %>%
  filter(year == elk_first_year)

elk_last_stage_distribution <- elk_stage_distribution %>%
  filter(year == elk_last_year)

################################################################################
#########--------------- Wolf population change ----------------################
################################################################################

wolf_total_df <- wolf_N_summ %>%
  filter(stage == "Total wolves") %>%
  arrange(year) %>%
  drop_na(mean)

wolf_pup_df <- wolf_N_summ %>%
  filter(stage == "Pups") %>%
  arrange(year) %>%
  drop_na(mean)

wolf_adult_df <- wolf_N_summ %>%
  filter(stage == "Adults") %>%
  arrange(year) %>%
  drop_na(mean)

wolf_first_year <- min(
  wolf_total_df$year,
  na.rm = TRUE
)

wolf_last_year <- max(
  wolf_total_df$year,
  na.rm = TRUE
)

wolf_first_total <- wolf_total_df %>%
  filter(year == wolf_first_year) %>%
  pull(mean)

wolf_last_total <- wolf_total_df %>%
  filter(year == wolf_last_year) %>%
  pull(mean)

wolf_pct_change <- 100 * (
  wolf_last_total - wolf_first_total
) / wolf_first_total

wolf_lambda_geom <- (
  wolf_last_total / wolf_first_total
)^(1 / (wolf_last_year - wolf_first_year))

wolf_annual_pct_change <- 100 * (
  wolf_lambda_geom - 1
)

# Annual posterior-mean stage distribution
wolf_stage_distribution <- wolf_pup_df %>%
  select(
    year,
    pups = mean
  ) %>%
  left_join(
    wolf_adult_df %>%
      select(
        year,
        adults = mean
      ),
    by = "year"
  ) %>%
  mutate(
    stage_total =
      pups +
      adults,
    
    prop_pups =
      pups / stage_total,
    
    prop_adults =
      adults / stage_total
  )

wolf_first_stage_distribution <- wolf_stage_distribution %>%
  filter(year == wolf_first_year)

wolf_last_stage_distribution <- wolf_stage_distribution %>%
  filter(year == wolf_last_year)

################################################################################
#########------------- Population change summaries -------------################
################################################################################

cat(
  "\n",
  "----------------------------\n",
  "------  ELK SUMMARY  -------\n",
  "----------------------------\n\n",
  
  "Population change\n",
  "-----------------\n",
  "Years: ",
  elk_first_year,
  " to ",
  elk_last_year,
  "\n",
  
  "First-year total females: ",
  round(elk_first_total),
  "\n",
  
  "Last-year total females: ",
  round(elk_last_total),
  "\n",
  
  "Total percent change: ",
  round(elk_pct_change, 1),
  "%\n",
  
  "Geometric lambda: ",
  round(elk_lambda_geom, 3),
  "\n",
  
  "Annual percent change: ",
  round(elk_annual_pct_change, 1),
  "%\n\n",
  
  "Stage distribution\n",
  "------------------\n",
  
  "Yearlings: ",
  elk_first_year,
  " = ",
  round(
    100 * elk_first_stage_distribution$prop_yearling,
    1
  ),
  "%; ",
  elk_last_year,
  " = ",
  round(
    100 * elk_last_stage_distribution$prop_yearling,
    1
  ),
  "%\n",
  
  "Young adults: ",
  elk_first_year,
  " = ",
  round(
    100 * elk_first_stage_distribution$prop_young_adult,
    1
  ),
  "%; ",
  elk_last_year,
  " = ",
  round(
    100 * elk_last_stage_distribution$prop_young_adult,
    1
  ),
  "%\n",
  
  "Old adults: ",
  elk_first_year,
  " = ",
  round(
    100 * elk_first_stage_distribution$prop_old_adult,
    1
  ),
  "%; ",
  elk_last_year,
  " = ",
  round(
    100 * elk_last_stage_distribution$prop_old_adult,
    1
  ),
  "%\n\n",
  
  "----------------------------\n",
  "------  WOLF SUMMARY  ------\n",
  "----------------------------\n\n",
  
  "Population change\n",
  "-----------------\n",
  "Years: ",
  wolf_first_year,
  " to ",
  wolf_last_year,
  "\n",
  
  "First-year total wolves: ",
  round(wolf_first_total),
  "\n",
  
  "Last-year total wolves: ",
  round(wolf_last_total),
  "\n",
  
  "Total percent change: ",
  round(wolf_pct_change, 1),
  "%\n",
  
  "Geometric lambda: ",
  round(wolf_lambda_geom, 3),
  "\n",
  
  "Annual percent change: ",
  round(wolf_annual_pct_change, 1),
  "%\n\n",
  
  "Stage distribution\n",
  "------------------\n",
  
  "Pups: ",
  wolf_first_year,
  " = ",
  round(
    100 * wolf_first_stage_distribution$prop_pups,
    1
  ),
  "%; ",
  wolf_last_year,
  " = ",
  round(
    100 * wolf_last_stage_distribution$prop_pups,
    1
  ),
  "%\n",
  
  "Adults: ",
  wolf_first_year,
  " = ",
  round(
    100 * wolf_first_stage_distribution$prop_adults,
    1
  ),
  "%; ",
  wolf_last_year,
  " = ",
  round(
    100 * wolf_last_stage_distribution$prop_adults,
    1
  ),
  "%\n"
)

################################################################################
#########---------------- Optional inspection ------------------################
################################################################################

# Complete evidence-classification table
coef_evidence

# Annual stage distributions
elk_stage_distribution
wolf_stage_distribution

# First- and last-year stage distributions
elk_first_stage_distribution
elk_last_stage_distribution

wolf_first_stage_distribution
wolf_last_stage_distribution

################################################################################
################################################################################
################################################################################