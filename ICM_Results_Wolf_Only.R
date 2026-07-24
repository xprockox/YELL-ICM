### Integrated Community Model (ICM)
### Results exploration
### Last updated: July 24, 2026
### xprockox@gmail.com

################################################################################
#########------------------- Load packages ---------------------################
################################################################################

library(tidyverse)
library(MCMCvis)
library(coda)
library(cowplot)
library(popbio)
library(stringr)
library(nimble)

################################################################################
############------------------ Load data ---------------------##################
################################################################################

# load model results
load("data/outputs/ICM_parallel_output_2026-07-23_wolfOnly.RData")

# import elk data
elk_dat_n <- read.csv("data/elk_abundanceEstimates_stages.csv")
elk_dat_n$n_female <- elk_dat_n$n_cow + (elk_dat_n$n_calf / 2)
elk_dat_n <- elk_dat_n %>%
  filter(year %in% community_years) %>%
  arrange(match(year, community_years))

# select either "NR" or "full" for northern range vs. total YELL population
wolf_range <- "full"

# import wolf data
wolf_pop <- switch(
  wolf_range,
  NR = read.csv("data/wolf_nr_pop.csv"),
  full = read.csv("data/wolf_full_park_pop.csv"),
  stop("wolf_range must be 'NR' or 'full'")
)
wolf_pop <- wolf_pop %>%
  filter(seasonal.year %in% community_years) %>%
  arrange(match(seasonal.year, community_years))

regression_years <- setdiff(community_years, drop_regression_years)

################################################################################
#########---------------------- Settings -----------------------################
################################################################################

# coefficient density plot selectors
# - "all" = all coefficients in that group
# - "intercepts" = only intercepts
# - "non_intercepts" = everything except intercepts
# - character vector = exact pretty labels to keep
# - regex:<pattern> = regex matched against pretty labels e.g., "regex:YA"

dens_coefs_elk <- "regex:(wolf|grizzly|cougar)"
dens_coefs_wolf <- "all"

################################################################################
#########---------------- Data cleaning (labels) ----------------###############
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

elk_vrate_map <- c(
  elk_s_c = "Calf survival (s_c)",
  elk_s_ya = "Young Adult survival (s_ya)",
  elk_s_oa = "Old Adult survival (s_oa)",
  elk_p_13 = "Young→Old transition (p_13)",
  elk_f_ya = "Fecundity (young) (f_ya)",
  elk_f_oa = "Fecundity (old) (f_oa)"
)

wolf_vrate_map <- c(
  wolf_s_p = "Pup survival (s_p)",
  wolf_s_a = "Adult survival (s_a)",
  wolf_f = "Fecundity (f)"
)

pretty_coef_labels <- c(
  beta0_calfSurv = "Calf survival intercept",
  beta1_calfSurv_wolfN = "Calf survival wolf effect",
  beta2_calfSurv_winterSeverity = "Calf survival winter severity effect",
  beta3_calfSurv_grizN = "Calf survival grizzly effect",
  beta4_calfSurv_cougarN = 'Calf survival cougar effect',
  beta5_calfSurv_elkN = "Calf survival density dependence effect",
  
  beta0_yaSurv = "YA survival intercept",
  beta1_yaSurv_wolfN = "YA survival wolf effect",
  beta2_yaSurv_winterSeverity = "YA survival winter severity effect",
  beta4_yaSurv_cougarN = 'YA survival cougar effect',
  beta5_yaSurv_elkN = "YA survival density dependence effect",
  beta6_yaSurv_annualNpp = "YA survival annual NPP effect",
  beta7_yaSurv_browndown = "YA survival browndown effect",
  beta8_yaSurv_pdsi = "YA survival PDSI effect",
  
  beta0_oaSurv = "OA survival intercept",
  beta1_oaSurv_wolfN = "OA survival wolf effect",
  beta2_oaSurv_winterSeverity = "OA survival winter severity effect",
  beta4_oaSurv_cougarN = 'OA survival cougar effect',
  beta5_oaSurv_elkN = "OA survival density dependence effect",
  beta6_oaSurv_annualNpp = "OA survival annual NPP effect",
  beta7_oaSurv_browndown = "OA survival browndown effect",
  beta8_oaSurv_pdsi = "OA survival PDSI effect",
  
  beta0_wpupSurv = "Wolf pup survival intercept",
  beta1_wpupSurv_elkN = "Wolf pup survival elk effect",
  beta2_wpupSurv_bisonN = "Wolf pup survival bison effect",
  beta3_wpupSurv_wolfN = "Wolf pup survival density dependence effect",
  
  beta0_wadSurv = "Wolf adult survival intercept",
  beta1_wadSurv_elkN = "Wolf adult survival elk effect",
  beta2_wadSurv_bisonN = "Wolf adult survival bison effect",
  beta3_wadSurv_wolfN = "Wolf adult survival density dependence effect",
  
  beta1_griz_elkCalves = "Grizzly elk calves effect",
  beta1_bison_cull = 'Bison cull/harvest effect',
  beta1_cougar_elk = 'Cougar elk effect'
)

################################################################################
#########------------------ Helper functions --------------------###############
################################################################################

clean_mcmc_samples <- function(chain_samples) {
  post_mat_raw <- as.matrix(mcmc.list(chain_samples))
  
  bad_cols <- which(apply(post_mat_raw, 2, function(x) any(!is.finite(x))))
  bad_summary <- tibble(
    parameter = colnames(post_mat_raw)[bad_cols],
    n_bad = sapply(bad_cols, function(i) sum(!is.finite(post_mat_raw[, i]))),
    n_total = nrow(post_mat_raw)
  )
  
  good_cols <- colnames(post_mat_raw)[apply(post_mat_raw, 2, function(x) all(is.finite(x)))]
  
  icm_clean <- mcmc.list(lapply(chain_samples, function(ch) {
    mcmc(as.matrix(ch)[, good_cols, drop = FALSE])
  }))
  
  list(
    icm_clean = icm_clean,
    post_mat = as.matrix(icm_clean),
    bad_summary = bad_summary
  )
}

extract_indexed_summary <- function(mcmc_obj, params, years, label_map = NULL,
                                    fecundity_params = NULL) {
  out <- MCMCsummary(mcmc_obj, params = params) %>%
    as.data.frame() %>%
    rownames_to_column("param") %>%
    rename(mean = mean, low = `2.5%`, high = `97.5%`) %>%
    mutate(
      param_base = str_extract(param, "^[^\\[]+"),
      year_index = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
    )
  
  if (!is.null(label_map)) {
    out <- out %>%
      mutate(param_base = recode(param_base, !!!label_map))
  }
  
  if (is.null(fecundity_params)) {
    out <- out %>%
      mutate(year = years[year_index])
  } else {
    out <- out %>%
      mutate(
        year = case_when(
          param_base %in% fecundity_params ~ years[-length(years)][year_index],
          TRUE ~ years[year_index]
        )
      )
  }
  
  out
}

extract_abundance_summary <- function(mcmc_obj, params, years, stage_map) {
  MCMCsummary(mcmc_obj, params = params) %>%
    as.data.frame() %>%
    rownames_to_column("param") %>%
    mutate(
      stage = str_extract(param, "^[^\\[]+"),
      t = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
    ) %>%
    transmute(
      stage = recode(stage, !!!stage_map),
      year = years[t],
      mean = mean,
      low = `2.5%`,
      high = `97.5%`
    )
}

get_rhat_issues <- function(summary_df, threshold = 1.1) {
  tibble(
    param = rownames(summary_df),
    Rhat = summary_df[, "Rhat"]
  ) %>%
    filter(Rhat > threshold)
}

extract_param_draws <- function(post_mat, prefix) {
  cols <- grep(paste0("^", prefix, "\\["), colnames(post_mat), value = TRUE)
  post_mat[, cols, drop = FALSE]
}

summarize_indexed_draws <- function(post_mat, prefix, years, value_name) {
  draws <- extract_param_draws(post_mat, prefix)
  
  tibble(
    year = years,
    !!value_name := apply(draws, 2, mean),
    low = apply(draws, 2, quantile, probs = 0.025, na.rm = TRUE),
    high = apply(draws, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

make_observed_long <- function(data, year_col, cols_map) {
  data %>%
    select(all_of(c(year_col, names(cols_map)))) %>%
    rename(year = all_of(year_col)) %>%
    pivot_longer(-year, names_to = "stage", values_to = "value") %>%
    mutate(stage = recode(stage, !!!cols_map))
}

plot_validation <- function(summary_df, observed_df, title, obs_stage = NULL) {
  p <- ggplot(summary_df, aes(x = year, y = mean, group = stage)) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = stage), alpha = 0.2) +
    geom_line(linewidth = 1) +
    facet_wrap(~stage, scales = "free_y") +
    theme_bw() +
    labs(
      x = "Year",
      y = "Abundance",
      title = title,
      subtitle = "Ribbon = 95% credible interval, line = posterior mean, red = observed"
    ) +
    theme(legend.position = "none")
  
  if (is.null(obs_stage)) {
    p <- p +
      geom_point(data = observed_df, aes(y = value), color = "red", size = 2) +
      geom_line(data = observed_df, aes(y = value), color = "red", linetype = 2)
  } else {
    p <- p +
      geom_point(
        data = observed_df %>% filter(stage == obs_stage),
        aes(y = value),
        color = "red",
        size = 2
      ) +
      geom_line(
        data = observed_df %>% filter(stage == obs_stage),
        aes(y = value),
        color = "red",
        linetype = 2
      )
  }
  
  p
}

plot_vital_rates <- function(vrate_df, title) {
  ggplot(vrate_df, aes(x = year, y = mean)) +
    geom_ribbon(aes(ymin = low, ymax = high, fill = param_base), alpha = 0.2) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~param_base, scales = "free_y") +
    theme_minimal() +
    labs(
      x = "Year",
      y = "Estimated value",
      title = title
    ) +
    theme(legend.position = "none")
}

resolve_coef_keep <- function(coef_names, label_map, keep = "all") {
  pretty_names <- unname(label_map[coef_names])
  pretty_names <- pretty_names[!is.na(pretty_names)]
  
  if (length(keep) == 1) {
    if (identical(keep, "all")) {
      return(pretty_names)
    }
    
    if (identical(keep, "intercepts")) {
      return(pretty_names[str_detect(pretty_names, regex("intercept", ignore_case = TRUE))])
    }
    
    if (identical(keep, "non_intercepts")) {
      return(pretty_names[!str_detect(pretty_names, regex("intercept", ignore_case = TRUE))])
    }
    
    if (str_detect(keep, "^regex:")) {
      pattern <- str_remove(keep, "^regex:")
      return(pretty_names[str_detect(pretty_names, regex(pattern, ignore_case = TRUE))])
    }
  }
  
  pretty_names[pretty_names %in% keep]
}

################################################################################
#########--------------- Clean posterior objects ---------------################
################################################################################

cleaned <- clean_mcmc_samples(icm_samples)
icm_clean <- cleaned$icm_clean
post_mat <- cleaned$post_mat
bad_summary <- cleaned$bad_summary

icm_summary <- MCMCsummary(icm_clean)
post_sum <- as.data.frame(icm_summary) %>%
  rownames_to_column("param")

bad_summary
round(icm_summary, 2)

################################################################################
#########---------------------- Core summaries ----------------------###########
################################################################################

elk_N_summ <- extract_abundance_summary(
  icm_clean,
  params = c("elk_N_1y", "elk_N_ya", "elk_N_oa", "elk_N_female"),
  years = community_years,
  stage_map = elk_stage_map
)

wolf_N_summ <- extract_abundance_summary(
  icm_clean,
  params = c("wolf_N_p", "wolf_N_a", "wolf_N_tot"),
  years = community_years,
  stage_map = wolf_stage_map
)

griz_N_summ <- extract_abundance_summary(
  icm_clean,
  params = "griz_N",
  years = community_years,
  stage_map = c(griz_N = "Grizzly abundance")
)

bison_N_summ <- extract_abundance_summary(
  icm_clean,
  params = "bison_N",
  years = community_years,
  stage_map = c(bison_N = "NR bison abundance")
)

cougar_N_summ <- extract_abundance_summary(
  icm_clean,
  params = "cougar_N",
  years = community_years,
  stage_map = c(cougar_N = "Cougar abundance")
)

elk_vrates2 <- extract_indexed_summary(
  icm_clean,
  params = c("elk_s_c", "elk_s_ya", "elk_s_oa", "elk_p_13", "elk_f_ya", "elk_f_oa"),
  years = community_years,
  label_map = elk_vrate_map,
  fecundity_params = c("Fecundity (young) (f_ya)", "Fecundity (old) (f_oa)")
)

wolf_vrates2 <- extract_indexed_summary(
  icm_clean,
  params = c("wolf_s_p", "wolf_s_a", "wolf_f"),
  years = community_years,
  label_map = wolf_vrate_map,
  fecundity_params = c("Fecundity (f)")
)

griz_logN_summ <- extract_indexed_summary(
  icm_clean,
  params = "griz_logN",
  years = community_years
)

bison_logN_summ <- extract_indexed_summary(
  icm_clean,
  params = "bison_logN",
  years = community_years
)

cougar_logN_summ <- extract_indexed_summary(
  icm_clean,
  params = "cougar_logN",
  years = community_years
) %>%
  filter(year %in% community_years[-1])

griz_mu_summ <- extract_indexed_summary(
  icm_clean,
  params = "griz_mu",
  years = community_years
) %>%
  filter(year %in% community_years[-1])

bison_mu_summ <- extract_indexed_summary(
  icm_clean,
  params = "bison_mu",
  years = community_years
) %>%
  filter(year %in% community_years[-1])

cougar_mu_summ <- extract_indexed_summary(
  icm_clean,
  params = "cougar_mu",
  years = community_years
) %>%
  filter(year %in% community_years[-1])

elk_N_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = c("elk_N_1y", "elk_N_ya", "elk_N_oa", "elk_N_female")
))

elk_vrate_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = c("elk_s_c", "elk_s_ya", "elk_s_oa", "elk_p_13", "elk_f_ya", "elk_f_oa")
))

wolf_N_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = c("wolf_N_p", "wolf_N_a", "wolf_N_tot")
))

wolf_vrate_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = c("wolf_s_p", "wolf_s_a", "wolf_f")
))

griz_coef_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = "beta1_griz_elkCalves"
))

bison_coef_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = "beta1_bison_cull"
))

cougar_coef_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = 'beta1_cougar_elk'
))

reg_coef_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = c(
    "beta0_calfSurv", 
    "beta1_calfSurv_wolfN", 
    "beta2_calfSurv_winterSeverity", 
    "beta4_calfSurv_cougarN",
    "beta3_calfSurv_grizN", 
    "beta5_calfSurv_elkN",
    
    "beta0_yaSurv", 
    "beta1_yaSurv_wolfN",
    "beta2_yaSurv_winterSeverity", 
    "beta4_yaSurv_cougarN",
    "beta5_yaSurv_elkN",
    "beta6_yaSurv_annualNpp", 
    "beta7_yaSurv_browndown", 
    "beta8_yaSurv_pdsi",
    
    "beta0_oaSurv", 
    "beta1_oaSurv_wolfN", 
    "beta2_oaSurv_winterSeverity", 
    "beta4_oaSurv_cougarN",
    "beta5_oaSurv_elkN",
    "beta6_oaSurv_annualNpp", 
    "beta7_oaSurv_browndown", 
    "beta8_oaSurv_pdsi",
    
    "beta0_wpupSurv", 
    "beta1_wpupSurv_elkN", 
    "beta2_wpupSurv_bisonN", 
    "beta3_wpupSurv_wolfN",
    
    "beta0_wadSurv", 
    "beta1_wadSurv_elkN", 
    "beta2_wadSurv_bisonN", 
    "beta3_wadSurv_wolfN"
  )
))

elk_N_rhat
elk_vrate_rhat
wolf_N_rhat
wolf_vrate_rhat
griz_coef_rhat
bison_coef_rhat
cougar_coef_rhat
reg_coef_rhat

################################################################################
##########----------------- Validation / time series -----------------##########
################################################################################

elk_obs_long <- make_observed_long(
  elk_dat_n,
  year_col = "year",
  cols_map = c(
    n_calf = "Yearling",
    n_cow_youngadult = "Young Adult",
    n_cow_oldadult = "Old Adult",
    n_female = "Total Females"
  )
) %>%
  filter(stage %in% c("Yearling", "Young Adult", "Old Adult", "Total Females"))

wolf_obs_long <- make_observed_long(
  wolf_pop,
  year_col = "seasonal.year",
  cols_map = c(
    dec_pups = "Pups",
    dec_adults = "Adults",
    total_abundance = "Total Wolves"
  )
)

bison_obs_long <- make_observed_long(
  bison,
  year_col = "year",
  cols_map = c(
    NR_Bison = "NR bison abundance"
  )
)

griz_obs_long <- make_observed_long(
  grizzly,
  year_col = "year",
  cols_map = c(
    griz_N = "Grizzly abundance"
  )
)

cougar_obs_long <- make_observed_long(
  cougars,
  year_col = 'year',
  cols_map = c(
    cougar_N = "Cougar abundance"
  )
)

elk_validation_plot <- plot_validation(
  elk_N_summ,
  elk_obs_long,
  title = "Elk posterior population estimates with validation data",
  obs_stage = "Total Females"
)

wolf_validation_plot <- plot_validation(
  wolf_N_summ,
  wolf_obs_long,
  title = "Wolf posterior population estimates with validation data"
)

griz_validation_plot <- plot_validation(
  griz_N_summ,
  griz_obs_long,
  title = "Grizzly posterior abundance estimates with observed data"
)

bison_validation_plot <- plot_validation(
  bison_N_summ,
  bison_obs_long,
  title = "Bison posterior abundance estimates with observed data"
)

cougar_validation_plot <- plot_validation(
  cougar_N_summ,
  cougar_obs_long,
  title = 'Cougar posterior abundance estimates with observed data'
)  

elk_vrate_plot <- plot_vital_rates(
  elk_vrates2,
  title = "Elk posterior time-varying vital rates (95% credible intervals)"
)

wolf_vrate_plot <- plot_vital_rates(
  wolf_vrates2,
  title = "Wolf posterior time-varying vital rates (95% credible intervals)"
)

# abundance validation
elk_validation_plot
wolf_validation_plot
griz_validation_plot
bison_validation_plot
cougar_validation_plot

elk_vrate_plot
wolf_vrate_plot

################################################################################
##########----------------- Plotting data builders -----------------############
################################################################################

# elk survival points
elk_surv_pts <- bind_rows(
  elk_vrates2 %>%
    filter(param_base == "Calf survival (s_c)", year %in% regression_years) %>%
    transmute(year, elk_surv = mean, elk_low = low, elk_high = high, stage = "Calf survival"),
  elk_vrates2 %>%
    filter(param_base == "Young Adult survival (s_ya)", year %in% regression_years) %>%
    transmute(year, elk_surv = mean, elk_low = low, elk_high = high, stage = "Young adult survival"),
  elk_vrates2 %>%
    filter(param_base == "Old Adult survival (s_oa)", year %in% regression_years) %>%
    transmute(year, elk_surv = mean, elk_low = low, elk_high = high, stage = "Old adult survival")
)

# wolf survival points
wolf_surv_pts <- bind_rows(
  wolf_vrates2 %>%
    filter(param_base == "Pup survival (s_p)", year %in% regression_years) %>%
    transmute(year, wolf_surv = mean, wolf_low = low, wolf_high = high, stage = "Pup survival"),
  wolf_vrates2 %>%
    filter(param_base == "Adult survival (s_a)", year %in% regression_years) %>%
    transmute(year, wolf_surv = mean, wolf_low = low, wolf_high = high, stage = "Adult survival")
)

# wolf abundance points used in elk plots
wolf_pts_elk <- summarize_indexed_draws(post_mat, "wolf_N_tot", community_years, "wolf_N_tot") %>%
  filter(year %in% regression_years) %>%
  rename(wolf_low = low, wolf_high = high)

# wolf abundance points used in wolf plots
wolf_pts_wolf <- summarize_indexed_draws(post_mat, "wolf_N_tot", community_years, "wolf_N_tot") %>%
  filter(year %in% regression_years) %>%
  rename(wolfN_low = low, wolfN_high = high)

# elk abundance points
elk_abund_pts <- summarize_indexed_draws(post_mat, "elk_N_female", community_years, "elk_N_female") %>%
  filter(year %in% regression_years) %>%
  rename(elk_low_x = low, elk_high_x = high)

# grizzly abundance points
griz_pts <- summarize_indexed_draws(post_mat, "griz_N", community_years, "griz_N") %>%
  filter(year %in% regression_years) %>%
  rename(griz_low = low, griz_high = high)

# elk calves points for grizzly state-space
elk_calf_griz_pts <- summarize_indexed_draws(
  post_mat, "elk_calves_born", community_years[-length(community_years)], "elk_calves_born"
) %>%
  filter(year %in% regression_years) %>%
  rename(x_low = low, x_high = high)

# bison culled points for bison state-space
bison_cull_pts <- covars %>%
  filter(year %in% regression_years) %>%
  select(year, total_cull_harvest)

# cougar abundance points
cougar_pts <- summarize_indexed_draws(post_mat, "cougar_N", community_years, "cougar_N") %>%
  filter(year %in% regression_years) %>%
  rename(cougar_low = low, cougar_high = high)

# elk abundance points for cougar state-space
elk_cougar_pts <- summarize_indexed_draws(
  post_mat, "elk_N_female", community_years, "elk_N_female"
) %>%
  filter(year %in% community_years[-length(community_years)]) %>%
  rename(x_low = low, x_high = high)

# create plotting dataframe for elk survival
plot_df_elk <- elk_surv_pts %>%
  left_join(wolf_pts_elk, by = "year") %>%
  left_join(elk_abund_pts, by = "year") %>%
  left_join(griz_pts, by = "year") %>%
  left_join(cougar_pts, by = "year") %>%
  left_join(
    covars %>%
      filter(year %in% regression_years) %>%
      select(
        year,
        winter_severity,
        annual_npp,
        browndown_onset_greenness_min,
        summer_avg_pdsi
      ),
    by = "year"
  )

# create plotting dataframe for wolf survival
wolf_plot_df <- wolf_surv_pts %>%
  left_join(
    summarize_indexed_draws(post_mat, "elk_N_female", community_years, "elk_N_female") %>%
      filter(year %in% regression_years) %>%
      rename(elk_low = low, elk_high = high),
    by = "year"
  ) %>%
  left_join(wolf_pts_wolf, by = "year") %>%
  left_join(
    covars %>%
      filter(year %in% regression_years) %>%
      select(year, NR_Bison),
    by = "year"
  )

################################################################################
##########-------------- Dynamic regression plot system ---------------#########
################################################################################

elk_model_specs <- list(
  "Calf survival" = list(
    intercept = "beta0_calfSurv",
    terms = c(
      wolf_N_tot = "beta1_calfSurv_wolfN",
      winter_severity = "beta2_calfSurv_winterSeverity",
      elk_N_female = "beta5_calfSurv_elkN"
    )
  ),
  
  "Young adult survival" = list(
    intercept = "beta0_yaSurv",
    terms = c(
      wolf_N_tot = "beta1_yaSurv_wolfN",
      winter_severity = "beta2_yaSurv_winterSeverity",
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
      elk_N_female = "beta5_oaSurv_elkN",
      annual_npp = "beta6_oaSurv_annualNpp",
      browndown_onset_greenness_min = "beta7_oaSurv_browndown",
      summer_avg_pdsi = "beta8_oaSurv_pdsi"
    )
  )
)

wolf_model_specs <- list(
  "Pup survival" = list(
    intercept = "beta0_wpupSurv",
    terms = c(
      elk_N_female = "beta1_wpupSurv_elkN",
      NR_Bison = "beta2_wpupSurv_bisonN",
      wolf_N_tot = "beta3_wpupSurv_wolfN"
    )
  ),
  "Adult survival" = list(
    intercept = "beta0_wadSurv",
    terms = c(
      elk_N_female = "beta1_wadSurv_elkN",
      NR_Bison = "beta2_wadSurv_bisonN",
      wolf_N_tot = "beta3_wadSurv_wolfN"
    )
  )
)

predictor_specs <- list(
  wolf_N_tot = list(
    raw_col = "wolf_N_tot",
    raw_df = plot_df_elk,
    std_fun = function(x) (x - icm_constants$wolf_tot_mean) / icm_constants$wolf_tot_sd,
    held_value = 0,
    xmin = "wolf_low",
    xmax = "wolf_high",
    label = "Wolf abundance"
  ),
  winter_severity = list(
    raw_col = "winter_severity",
    raw_df = covars %>% filter(year %in% regression_years),
    std_fun = function(x) x,
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
    label = "Winter severity"
  ),
  griz_N = list(
    raw_col = "griz_N",
    raw_df = griz_pts,
    std_fun = function(x) (x - icm_constants$griz_N_mean) / icm_constants$griz_N_sd,
    held_value = 0,
    xmin = "griz_low",
    xmax = "griz_high",
    label = "Grizzly abundance"
  ),
  elk_N_female = list(
    raw_col = "elk_N_female",
    raw_df = plot_df_elk,
    std_fun = function(x) (x - icm_constants$elk_N_female_mean) / icm_constants$elk_N_female_sd,
    held_value = 0,
    xmin = "elk_low_x",
    xmax = "elk_high_x",
    label = "Elk abundance"
  ),
  annual_npp = list(
    raw_col = "annual_npp",
    raw_df = covars %>% filter(year %in% regression_years),
    std_fun = function(x) (x - mean(covars$annual_npp, na.rm = TRUE)) / sd(covars$annual_npp, na.rm = TRUE),
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
    label = "Annual NPP"
  ),
  browndown_onset_greenness_min = list(
    raw_col = "browndown_onset_greenness_min",
    raw_df = covars %>% filter(year %in% regression_years),
    std_fun = function(x) (x - mean(covars$browndown_onset_greenness_min, na.rm = TRUE)) / sd(covars$browndown_onset_greenness_min, na.rm = TRUE),
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
    label = "Browndown"
  ),
  summer_avg_pdsi = list(
    raw_col = "summer_avg_pdsi",
    raw_df = covars %>% filter(year %in% regression_years),
    std_fun = function(x) (x - mean(covars$summer_avg_pdsi, na.rm = TRUE)) / sd(covars$summer_avg_pdsi, na.rm = TRUE),
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
    label = "Summer PDSI"
  ),
  NR_Bison = list(
    raw_col = "NR_Bison",
    raw_df = covars %>% filter(year %in% regression_years),
    std_fun = function(x) (x - mean(covars$NR_Bison, na.rm = TRUE)) / sd(covars$NR_Bison, na.rm = TRUE),
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
    label = "Bison abundance"
  ),
  cougar_N = list(
    raw_col = "cougar_N",
    raw_df = cougar_pts,
    std_fun = function(x) (x - icm_constants$cougar_N_mean) / icm_constants$cougar_N_sd,
    held_value = 0,
    xmin = "cougar_low",
    xmax = "cougar_high",
    label = "Cougar abundance"
  )
)

make_effect_curve <- function(post_mat, model_specs, focal_predictor, predictor_specs,
                              length_out = 200) {
  spec <- predictor_specs[[focal_predictor]]
  
  x_raw <- seq(
    min(spec$raw_df[[spec$raw_col]], na.rm = TRUE),
    max(spec$raw_df[[spec$raw_col]], na.rm = TRUE),
    length.out = length_out
  )
  x_std <- spec$std_fun(x_raw)
  
  out <- lapply(names(model_specs), function(stage_name) {
    stage_spec <- model_specs[[stage_name]]
    
    if (!focal_predictor %in% names(stage_spec$terms)) {
      return(NULL)
    }
    
    eta <- matrix(post_mat[, stage_spec$intercept], nrow = nrow(post_mat), ncol = length_out)
    
    for (pred_name in names(stage_spec$terms)) {
      beta_name <- stage_spec$terms[[pred_name]]
      pred_vals <- if (pred_name == focal_predictor) x_std else rep(predictor_specs[[pred_name]]$held_value, length_out)
      eta <- eta + outer(post_mat[, beta_name], pred_vals)
    }
    
    mu <- plogis(eta)
    
    tibble(
      x = x_raw,
      surv = apply(mu, 2, mean),
      low = apply(mu, 2, quantile, probs = 0.025, na.rm = TRUE),
      high = apply(mu, 2, quantile, probs = 0.975, na.rm = TRUE),
      stage = stage_name
    )
  })
  
  bind_rows(out)
}

plot_effect <- function(plot_df, curve_df, x_var, y_var, ymin_var, ymax_var,
                        xlow_var = NULL, xhigh_var = NULL,
                        x_label, y_label, title, subtitle,
                        fill_color, line_color) {
  p <- ggplot(plot_df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_ribbon(
      data = curve_df,
      aes(x = x, ymin = low, ymax = high),
      inherit.aes = FALSE,
      fill = fill_color,
      alpha = 0.25
    ) +
    geom_line(
      data = curve_df,
      aes(x = x, y = surv),
      inherit.aes = FALSE,
      color = line_color,
      linewidth = 1
    ) +
    geom_errorbar(
      aes(ymin = .data[[ymin_var]], ymax = .data[[ymax_var]]),
      width = 0
    ) +
    facet_wrap(~stage, scales = "free_y") +
    theme_classic() +
    labs(
      x = x_label,
      y = y_label,
      title = title,
      subtitle = subtitle
    )
  
  if (!is.null(xlow_var) && !is.null(xhigh_var)) {
    p <- p + geom_errorbarh(
      aes(xmin = .data[[xlow_var]], xmax = .data[[xhigh_var]]),
      height = 0,
      alpha = 0.6
    )
  }
  
  p +
    geom_point(shape = 21, fill = "gold", color = "black", size = 2)
}

################################################################################
############----------- State-space effect curve plots -------------############
################################################################################

make_bison_cull_effect_curve <- function(post_mat, cull_vals, ref_val = NULL) {
  beta_draws <- post_mat[, "beta1_bison_cull", drop = TRUE]
  
  if (is.null(ref_val)) {
    ref_val <- median(cull_vals, na.rm = TRUE)
  }
  
  x_trans <- log1p(cull_vals)
  ref_trans <- log1p(ref_val)
  
  effect_mat <- outer(beta_draws, x_trans - ref_trans)
  
  tibble(
    x = cull_vals,
    effect = apply(effect_mat, 2, mean),
    low = apply(effect_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    high = apply(effect_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

make_griz_elkCalf_effect_curve <- function(post_mat, calf_vals, ref_val = NULL) {
  beta_draws <- post_mat[, "beta1_griz_elkCalves", drop = TRUE]
  
  if (is.null(ref_val)) {
    ref_val <- median(calf_vals, na.rm = TRUE)
  }
  
  x_trans <- log(calf_vals)
  ref_trans <- log(ref_val)
  
  effect_mat <- outer(beta_draws, x_trans - ref_trans)
  
  tibble(
    x = calf_vals,
    effect = apply(effect_mat, 2, mean),
    low = apply(effect_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    high = apply(effect_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

make_cougar_elk_effect_curve <- function(post_mat, elk_vals, ref_val = NULL) {
  beta_draws <- post_mat[, "beta1_cougar_elk", drop = TRUE]
  
  if (is.null(ref_val)) {
    ref_val <- median(elk_vals, na.rm = TRUE)
  }
  
  x_trans <- log(elk_vals + 1e-6)
  ref_trans <- log(ref_val + 1e-6)
  
  effect_mat <- outer(beta_draws, x_trans - ref_trans)
  
  tibble(
    x = elk_vals,
    effect = apply(effect_mat, 2, mean),
    low = apply(effect_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    high = apply(effect_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

################################################################################
############------------ Reproducing regression plots --------------############
################################################################################

line_df_wolf_elk <- make_effect_curve(post_mat, elk_model_specs, "wolf_N_tot", predictor_specs)
# line_df_griz_elk <- make_effect_curve(post_mat, elk_model_specs, "griz_N", predictor_specs)
line_df_elkN_elk <- make_effect_curve(post_mat, elk_model_specs, "elk_N_female", predictor_specs)
line_df_winter_elk <- make_effect_curve(post_mat, elk_model_specs, "winter_severity", predictor_specs)
line_df_npp_elk <- make_effect_curve(post_mat, elk_model_specs, "annual_npp", predictor_specs)
line_df_browndown_elk <- make_effect_curve(post_mat, elk_model_specs, "browndown_onset_greenness_min", predictor_specs)
line_df_pdsi_elk <- make_effect_curve(post_mat, elk_model_specs, "summer_avg_pdsi", predictor_specs)

line_df_elk_wolf <- make_effect_curve(post_mat, wolf_model_specs, "elk_N_female", predictor_specs)
line_df_bison_wolf <- make_effect_curve(post_mat, wolf_model_specs, "NR_Bison", predictor_specs)
line_df_wolf_wolf <- make_effect_curve(post_mat, wolf_model_specs, "wolf_N_tot", predictor_specs)
# line_df_cougar_elk <- make_effect_curve(post_mat, elk_model_specs, "cougar_N", predictor_specs)

wolf_plot_elk <- plot_effect(
  plot_df = plot_df_elk,
  curve_df = line_df_wolf_elk,
  x_var = "wolf_N_tot",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  xlow_var = "wolf_low",
  xhigh_var = "wolf_high",
  x_label = "Wolf abundance",
  y_label = "Elk survival",
  title = "Estimated effect of wolf abundance on elk survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#6F263D",
  line_color = "#6F263D"
)

winter_plot_elk <- plot_effect(
  plot_df = plot_df_elk,
  curve_df = line_df_winter_elk,
  x_var = "winter_severity",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  x_label = "Winter severity",
  y_label = "Elk survival",
  title = "Estimated effect of winter severity on elk survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#236192",
  line_color = "#236192"
)

npp_plot_elk <- plot_effect(
  plot_df = plot_df_elk %>% filter(stage %in% c("Young adult survival", 'Old adult survival')),
  curve_df = line_df_npp_elk,
  x_var = "annual_npp",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  x_label = "NPP",
  y_label = "Elk survival",
  title = "Estimated effect of annual vegetation productivity (NPP) on elk survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#236192",
  line_color = "#236192"
)

browndown_plot_elk <- plot_effect(
  plot_df = plot_df_elk %>% filter(stage %in% c("Young adult survival", 'Old adult survival')),
  curve_df = line_df_browndown_elk,
  x_var = "browndown_onset_greenness_min",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  x_label = "Browndown Date",
  y_label = "Elk survival",
  title = "Estimated effect of browndown timing on elk survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#236192",
  line_color = "#236192"
)

pdsi_plot_elk <- plot_effect(
  plot_df = plot_df_elk %>% filter(stage %in% c("Young adult survival", 'Old adult survival')),
  curve_df = line_df_pdsi_elk,
  x_var = "summer_avg_pdsi",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  x_label = "PDSI",
  y_label = "Elk survival",
  title = "Estimated effect of drought (PDSI) on elk survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#236192",
  line_color = "#236192"
)

# griz_plot_elk <- plot_effect(
#   plot_df = plot_df_elk %>% filter(stage == "Calf survival"),
#   curve_df = line_df_griz_elk,
#   x_var = "griz_N",
#   y_var = "elk_surv",
#   ymin_var = "elk_low",
#   ymax_var = "elk_high",
#   x_label = "Grizzly abundance",
#   y_label = "Calf survival",
#   title = "Estimated effect of grizzly abundance on elk calf survival",
#   subtitle = "Other predictors held constant at their means",
#   fill_color = "#4B7F52",
#   line_color = "#4B7F52"
# )
# 
# cougar_plot_elk <- plot_effect(
#   plot_df = plot_df_elk,
#   curve_df = line_df_cougar_elk,
#   x_var = "cougar_N",
#   y_var = "elk_surv",
#   ymin_var = "elk_low",
#   ymax_var = "elk_high",
#   # xlow_var = "cougar_low",
#   # xhigh_var = "cougar_high",
#   x_label = "Cougar abundance",
#   y_label = "Elk survival",
#   title = "Estimated effect of cougar abundance on elk survival",
#   subtitle = "Other predictors held constant at their means",
#   fill_color = "#7A5C99",
#   line_color = "#7A5C99"
# )

elkN_plot_elk <- plot_effect(
  plot_df = plot_df_elk,
  curve_df = line_df_elkN_elk,
  x_var = "elk_N_female",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  # xlow_var = "elk_low_x",
  # xhigh_var = "elk_high_x",
  x_label = "Elk abundance",
  y_label = "Elk survival",
  title = "Estimated effect of density dependence on elk survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#8C510A",
  line_color = "#8C510A"
)

wolf_main_plot <- plot_effect(
  plot_df = wolf_plot_df,
  curve_df = line_df_elk_wolf,
  x_var = "elk_N_female",
  y_var = "wolf_surv",
  ymin_var = "wolf_low",
  ymax_var = "wolf_high",
  xlow_var = "elk_low",
  xhigh_var = "elk_high",
  x_label = "Elk female abundance",
  y_label = "Wolf survival",
  title = "Estimated effect of elk abundance on wolf survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#236192",
  line_color = "#236192"
)

bison_plot_wolf <- plot_effect(
  plot_df = wolf_plot_df,
  curve_df = line_df_bison_wolf,
  x_var = "NR_Bison",
  y_var = "wolf_surv",
  ymin_var = "wolf_low",
  ymax_var = "wolf_high",
  x_label = "Bison abundance",
  y_label = "Wolf survival",
  title = "Estimated effect of bison abundance on wolf survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#6F263D",
  line_color = "#6F263D"
)

wolfN_plot_wolf <- plot_effect(
  plot_df = wolf_plot_df,
  curve_df = line_df_wolf_wolf,
  x_var = "wolf_N_tot",
  y_var = "wolf_surv",
  ymin_var = "wolf_low",
  ymax_var = "wolf_high",
  xlow_var = "wolfN_low",
  xhigh_var = "wolfN_high",
  x_label = "Wolf abundance",
  y_label = "Wolf survival",
  title = "Estimated effect of density dependence on wolf survival",
  subtitle = "Other predictors held constant at their means",
  fill_color = "#4B7F52",
  line_color = "#4B7F52"
)

bison_cull_curve <- make_bison_cull_effect_curve(
  post_mat,
  seq(
    min(bison_cull_pts$total_cull_harvest, na.rm = TRUE),
    max(bison_cull_pts$total_cull_harvest, na.rm = TRUE),
    length.out = 200
  )
)

bison_cull_plot <- ggplot() +
  geom_ribbon(
    data = bison_cull_curve,
    aes(x = x, ymin = low, ymax = high),
    fill = "#8C510A",
    alpha = 0.25
  ) +
  geom_line(
    data = bison_cull_curve,
    aes(x = x, y = effect),
    color = "#8C510A",
    linewidth = 1
  ) +
  geom_point(
    data = bison_cull_pts,
    aes(x = total_cull_harvest, y = 0),
    shape = 21,
    fill = "gold",
    color = "black",
    size = 2
  ) +
  theme_classic() +
  labs(
    x = "Bison culled / harvested",
    y = "Change in expected log bison abundance",
    title = "Estimated effect of bison cull/harvest on bison population process",
    subtitle = "Curve shows change relative to the median observed cull/harvest level"
  )

griz_elkCalf_curve <- make_griz_elkCalf_effect_curve(
  post_mat,
  seq(
    min(elk_calf_griz_pts$elk_calves_born, na.rm = TRUE),
    max(elk_calf_griz_pts$elk_calves_born, na.rm = TRUE),
    length.out = 200
  )
)

griz_elkCalf_plot <- ggplot() +
  geom_ribbon(
    data = griz_elkCalf_curve,
    aes(x = x, ymin = low, ymax = high),
    fill = "#4B7F52",
    alpha = 0.25
  ) +
  geom_line(
    data = griz_elkCalf_curve,
    aes(x = x, y = effect),
    color = "#4B7F52",
    linewidth = 1
  ) +
  geom_point(
    data = elk_calf_griz_pts,
    aes(x = elk_calves_born, y = 0),
    shape = 21,
    fill = "gold",
    color = "black",
    size = 2
  ) +
  theme_classic() +
  labs(
    x = "Elk calves born",
    y = "Change in expected log grizzly abundance",
    title = "Estimated effect of elk calf abundance on grizzly population process",
    subtitle = "Curve shows change relative to the median observed elk calf abundance"
  )

cougar_elk_curve <- make_cougar_elk_effect_curve(
  post_mat,
  seq(
    min(elk_cougar_pts$elk_N_female, na.rm = TRUE),
    max(elk_cougar_pts$elk_N_female, na.rm = TRUE),
    length.out = 200
  )
)

cougar_elk_plot <- ggplot() +
  geom_ribbon(
    data = cougar_elk_curve,
    aes(x = x, ymin = low, ymax = high),
    fill = "#7A5C99",
    alpha = 0.25
  ) +
  geom_line(
    data = cougar_elk_curve,
    aes(x = x, y = effect),
    color = "#7A5C99",
    linewidth = 1
  ) +
  geom_point(
    data = elk_cougar_pts,
    aes(x = elk_N_female, y = 0),
    shape = 21,
    fill = "gold",
    color = "black",
    size = 2
  ) +
  theme_classic() +
  labs(
    x = "Elk female abundance",
    y = "Change in expected log cougar abundance",
    title = "Estimated effect of elk abundance on cougar population process",
    subtitle = "Curve shows change relative to the median observed elk abundance"
  )

# elk regressions
wolf_plot_elk
winter_plot_elk
# griz_plot_elk
# cougar_plot_elk
npp_plot_elk
browndown_plot_elk
pdsi_plot_elk
elkN_plot_elk
# yaHarvest_plot_elk
# oaHarvest_plot_elk

# wolf regressions
wolf_main_plot
bison_plot_wolf
wolfN_plot_wolf

# grizzlies
griz_elkCalf_plot

# bison
bison_cull_plot

# cougar
cougar_elk_plot

################################################################################
##########----------------- Coefficient density plots -----------------#########
################################################################################

make_coef_density_plot_stacked <- function(post_mat, coef_names, label_map,
                                           keep = "all",
                                           fill_color,
                                           title,
                                           effect_order = c(
                                             "Intercept",
                                             "Wolf abundance",
                                             "Grizzly abundance",
                                             "Cougar abundance",
                                             "Elk abundance",
                                             "Bison abundance",
                                             "Density dependence",
                                             "Winter severity",
                                             "Annual NPP",
                                             "Browndown",
                                             "PDSI"
                                           ),
                                           model_order = c("Calf survival", "YA survival", "OA survival",
                                                           "Wolf pup survival", "Wolf adult survival"),
                                           scale_factor = 0.35) {
  
  keep_resolved <- resolve_coef_keep(
    coef_names = coef_names,
    label_map = label_map,
    keep = keep
  )
  
  df <- as_tibble(post_mat[, coef_names, drop = FALSE]) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(
      pretty = recode(parameter, !!!label_map)
    ) %>%
    filter(pretty %in% keep_resolved) %>%
    mutate(
      model_group = case_when(
        str_detect(pretty, "^Calf survival") ~ "Calf survival",
        str_detect(pretty, "^YA survival") ~ "YA survival",
        str_detect(pretty, "^OA survival") ~ "OA survival",
        str_detect(pretty, "^Wolf pup survival") ~ "Wolf pup survival",
        str_detect(pretty, "^Wolf adult survival") ~ "Wolf adult survival",
        TRUE ~ "Other"
      ),
      effect_type = case_when(
        str_detect(pretty, regex("intercept", ignore_case = TRUE)) ~ "Intercept",
        str_detect(pretty, regex("grizzly effect", ignore_case = TRUE)) ~ "Grizzly abundance",
        str_detect(pretty, regex("bison effect", ignore_case = TRUE)) ~ "Bison abundance",
        str_detect(pretty, regex("cougar effect", ignore_case = TRUE)) ~ "Cougar abundance",
        str_detect(pretty, regex("annual NPP effect", ignore_case = TRUE)) ~ "Annual NPP",
        str_detect(pretty, regex("browndown effect", ignore_case = TRUE)) ~ "Browndown",
        str_detect(pretty, regex("PDSI effect", ignore_case = TRUE)) ~ "PDSI",
        str_detect(pretty, regex("winter severity effect", ignore_case = TRUE)) ~ "Winter severity",
        str_detect(pretty, regex("density dependence effect", ignore_case = TRUE)) ~ "Density dependence",
        str_detect(pretty, regex("wolf effect", ignore_case = TRUE)) ~ "Wolf abundance",
        str_detect(pretty, regex("elk effect", ignore_case = TRUE)) ~ "Elk abundance",
        TRUE ~ pretty
      )
    ) %>%
    mutate(
      effect_type = factor(effect_type, levels = effect_order),
      model_group = factor(model_group, levels = model_order)
    ) %>%
    drop_na(model_group, effect_type)
  
  dens_df <- df %>%
    group_by(model_group, effect_type) %>%
    group_modify(~ {
      d <- density(.x$value, na.rm = TRUE)
      tibble(x = d$x, density = d$y)
    }) %>%
    ungroup() %>%
    mutate(
      effect_num = as.numeric(effect_type),
      y = effect_num + density * scale_factor
    )
  
  label_df <- df %>%
    distinct(model_group, effect_type) %>%
    mutate(effect_num = as.numeric(effect_type))
  
  ggplot() +
    geom_ribbon(
      data = dens_df,
      aes(x = x, ymin = effect_num, ymax = y, group = interaction(model_group, effect_type)),
      fill = fill_color,
      alpha = 0.45
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = y, group = interaction(model_group, effect_type)),
      linewidth = 0.5
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(
      data = label_df,
      aes(yintercept = effect_num),
      color = "grey80",
      linewidth = 0.3
    ) +
    facet_wrap(~model_group, nrow = 1, scales = "free_x") +
    scale_y_continuous(
      breaks = label_df$effect_num,
      labels = label_df$effect_type,
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    theme_bw() +
    labs(
      x = "Posterior value",
      y = NULL,
      title = title
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}

elk_coef_plot <- make_coef_density_plot_stacked(
  post_mat = post_mat,
  coef_names = c(
    # "beta0_calfSurv", 
    "beta1_calfSurv_wolfN",
    # "beta2_calfSurv_winterSeverity", 
    # "beta3_calfSurv_grizN",
    # "beta4_calfSurv_cougarN",
    # "beta5_calfSurv_elkN",
    
    # "beta0_yaSurv", 
    "beta1_yaSurv_wolfN", 
    # "beta2_yaSurv_winterSeverity",
    # "beta4_yaSurv_cougarN",
    # "beta5_yaSurv_elkN",
    # "beta6_yaSurv_annualNpp", 
    # "beta7_yaSurv_browndown", 
    # "beta8_yaSurv_pdsi",
    
    # "beta0_oaSurv", 
    "beta1_oaSurv_wolfN"#, 
    # "beta2_oaSurv_winterSeverity",
    # "beta4_oaSurv_cougarN",
    # "beta5_oaSurv_elkN",
    # "beta6_oaSurv_annualNpp", 
    # "beta7_oaSurv_browndown", 
    # "beta8_oaSurv_pdsi"
  ),
  label_map = pretty_coef_labels,
  keep = dens_coefs_elk,
  fill_color = "#236192",
  title = "Posterior distributions of elk regression coefficients"
)

wolf_coef_plot <- make_coef_density_plot_stacked(
  post_mat = post_mat,
  coef_names = c(
    "beta0_wpupSurv", 
    "beta1_wpupSurv_elkN", 
    "beta2_wpupSurv_bisonN", 
    "beta3_wpupSurv_wolfN",
    
    "beta0_wadSurv", 
    "beta1_wadSurv_elkN", 
    "beta2_wadSurv_bisonN",
    "beta3_wadSurv_wolfN"
  ),
  label_map = pretty_coef_labels,
  keep = dens_coefs_wolf,
  fill_color = "#6F263D",
  title = "Posterior distributions of wolf regression coefficients",
  effect_order = c(
    "Intercept",
    "Elk abundance",
    "Bison abundance",
    "Density dependence"
  ),
  model_order = c("Wolf pup survival", "Wolf adult survival")
)

griz_coef_plot <- tibble(value = post_mat[, "beta1_griz_elkCalves"]) %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "#4B7F52", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distribution of elk calves effect on grizzly abundance"
  )

bison_coef_plot <- tibble(value = post_mat[, "beta1_bison_cull"]) %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "#8C510A", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distribution of bison cull/harvest effect on bison abundance"
  )

cougar_coef_plot <- tibble(value = post_mat[, "beta1_cougar_elk"]) %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "#7A5C99", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distribution of elk effect on cougar abundance"
  )

# view plots
elk_coef_plot
wolf_coef_plot
griz_coef_plot
bison_coef_plot
cougar_coef_plot

################################################################################
#########----------------- Evidence classification ----------------#############
################################################################################

classify_evidence <- function(draws) {
  ci95 <- quantile(draws, probs = c(0.025, 0.975), na.rm = TRUE)
  ci80 <- quantile(draws, probs = c(0.10, 0.90), na.rm = TRUE)
  ci50 <- quantile(draws, probs = c(0.25, 0.75), na.rm = TRUE)
  
  excl95 <- ci95[1] > 0 | ci95[2] < 0
  excl80 <- ci80[1] > 0 | ci80[2] < 0
  excl50 <- ci50[1] > 0 | ci50[2] < 0
  
  tibble(
    mean = mean(draws, na.rm = TRUE),
    median = median(draws, na.rm = TRUE),
    ci95_low = ci95[1],
    ci95_high = ci95[2],
    ci80_low = ci80[1],
    ci80_high = ci80[2],
    ci50_low = ci50[1],
    ci50_high = ci50[2],
    direction = case_when(
      median(draws, na.rm = TRUE) > 0 ~ "Positive",
      median(draws, na.rm = TRUE) < 0 ~ "Negative",
      TRUE ~ "Neutral"
    ),
    evidence = case_when(
      excl95 ~ "Strong",
      !excl95 & excl80 ~ "Moderate",
      !excl80 & excl50 ~ "Weak",
      TRUE ~ "Little/none"
    )
  )
}

coef_names <- names(pretty_coef_labels)[names(pretty_coef_labels) %in% colnames(post_mat)]

coef_evidence <- bind_rows(lapply(coef_names, function(param) {
  classify_evidence(post_mat[, param]) %>%
    mutate(parameter = param, label = pretty_coef_labels[[param]])
})) %>%
  select(parameter, label, everything()) %>%
  arrange(factor(evidence, levels = c("Strong", "Moderate", "Weak", "Little/none")), desc(abs(median)))

strong_effects <- coef_evidence %>% filter(evidence == "Strong")
moderate_effects <- coef_evidence %>% filter(evidence == "Moderate")
weak_effects <- coef_evidence %>% filter(evidence == "Weak")
little_effects <- coef_evidence %>% filter(evidence == "Little/none")

griz_coef_evidence <- classify_evidence(post_mat[, "beta1_griz_elkCalves"]) %>%
  mutate(
    parameter = "beta1_griz_elkCalves",
    label = "Grizzly elk calves effect"
  )

bison_coef_evidence <- classify_evidence(post_mat[, "beta1_bison_cull"]) %>%
  mutate(
    parameter = "beta1_bison_cull",
    label = "Bison cull/harvest effect"
  )

cougar_coef_evidence <- classify_evidence(post_mat[, "beta1_cougar_elk"]) %>%
  mutate(
    parameter = "beta1_cougar_elk",
    label = "Cougar elk effect"
  )

################################################################################
##########----------------- Elasticity analysis -----------------###############
################################################################################

elk_rates_wide <- elk_vrates2 %>%
  filter(param_base %in% c(
    "Calf survival (s_c)",
    "Young Adult survival (s_ya)",
    "Old Adult survival (s_oa)",
    "Young→Old transition (p_13)",
    "Fecundity (young) (f_ya)",
    "Fecundity (old) (f_oa)"
  )) %>%
  select(year, param_base, mean) %>%
  pivot_wider(names_from = param_base, values_from = mean) %>%
  drop_na()

elk_elasticity_df <- bind_rows(lapply(seq_len(nrow(elk_rates_wide)), function(i) {
  s_c <- as.numeric(elk_rates_wide$`Calf survival (s_c)`[i])
  s_ya <- as.numeric(elk_rates_wide$`Young Adult survival (s_ya)`[i])
  s_oa <- as.numeric(elk_rates_wide$`Old Adult survival (s_oa)`[i])
  p_13 <- as.numeric(elk_rates_wide$`Young→Old transition (p_13)`[i])
  f_ya <- as.numeric(elk_rates_wide$`Fecundity (young) (f_ya)`[i])
  f_oa <- as.numeric(elk_rates_wide$`Fecundity (old) (f_oa)`[i])
  
  A <- matrix(
    c(
      0, f_ya * s_c, f_oa * s_c,
      s_ya, s_ya * (1 - p_13), 0,
      0, s_ya * p_13, s_oa
    ),
    nrow = 3,
    byrow = TRUE
  )
  
  E <- popbio::elasticity(A)
  
  tibble(
    year = elk_rates_wide$year[i],
    lambda = popbio::lambda(A),
    a12 = E[1, 2],
    a13 = E[1, 3],
    a21 = E[2, 1],
    a22 = E[2, 2],
    a32 = E[3, 2],
    a33 = E[3, 3]
  )
}))

elk_elasticity_long <- elk_elasticity_df %>%
  pivot_longer(c(a12, a13, a21, a22, a32, a33), names_to = "transition", values_to = "elasticity") %>%
  mutate(
    transition = factor(
      transition,
      levels = c("a12", "a13", "a21", "a22", "a32", "a33"),
      labels = c(
        "Yearling recruitment from young adults",
        "Yearling recruitment from old adults",
        "Yearling → young adult",
        "Young adult survival (non-transition)",
        "Young adult → old adult",
        "Old adult survival (non-transition)"
      )
    )
  )

lambda_plot <- ggplot(elk_elasticity_df, aes(x = year, y = lambda)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    x = "Year",
    y = expression(lambda),
    title = "Annual elk population growth rate from projection matrix",
    subtitle = "Dashed line at lambda = 1 indicates stability"
  )

elasticity_bar_plot <- elk_elasticity_long %>%
  group_by(transition) %>%
  summarise(mean_elasticity = mean(elasticity, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_elasticity)) %>%
  ggplot(aes(x = reorder(transition, mean_elasticity), y = mean_elasticity)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean elasticity",
    title = "Mean elasticity of elk population growth to matrix transitions"
  )

elasticity_combo <- plot_grid(
  lambda_plot,
  elasticity_bar_plot,
  ncol = 1,
  rel_heights = c(1, 1.1),
  labels = c("A", "B")
)

################################################################################
############---------------------- Final plots ----------------------###########
################################################################################

# abundance validation
elk_validation_plot
wolf_validation_plot
griz_validation_plot
bison_validation_plot
cougar_validation_plot

# vital rates
elk_vrate_plot
wolf_vrate_plot

# elk regressions
wolf_plot_elk
winter_plot_elk
# griz_plot_elk
# cougar_plot_elk
elkN_plot_elk
npp_plot_elk
browndown_plot_elk
pdsi_plot_elk
elk_coef_plot

# wolf regressions
wolf_main_plot
bison_plot_wolf
wolfN_plot_wolf
wolf_coef_plot

# grizzlies
griz_elkCalf_plot
griz_coef_plot

# bison
bison_cull_plot
bison_coef_plot

# cougars
cougar_coef_plot

# elasticity
lambda_plot
elasticity_bar_plot
elasticity_combo

################################################################################
##########-------------- Population change summaries --------------#############
################################################################################

# calculate summary statistics for elk
elk_total_df <- elk_N_summ %>% 
  filter(stage == "Total Females") %>% 
  arrange(year) %>%
  drop_na()

elk_yearling_df <- elk_N_summ %>% 
  filter(stage == "Yearling") %>% 
  arrange(year) %>%
  drop_na()

elk_ya_df <- elk_N_summ %>% 
  filter(stage == "Young Adult") %>% 
  arrange(year) %>%
  drop_na()

elk_oa_df <- elk_N_summ %>% 
  filter(stage == "Old Adult") %>% 
  arrange(year) %>%
  drop_na()

elk_first_year <- min(elk_total_df$year, na.rm = TRUE)
elk_last_year <- max(elk_total_df$year, na.rm = TRUE)

elk_first_total <- elk_total_df$mean[elk_total_df$year == elk_first_year]
elk_last_total <- elk_total_df$mean[elk_total_df$year == elk_last_year]

elk_pct_change <- 100 * (elk_last_total - elk_first_total) / elk_first_total
elk_lambda_geom <- (elk_last_total / elk_first_total)^(1 / (elk_last_year - elk_first_year))
elk_annual_pct_change <- 100 * (elk_lambda_geom - 1)

elk_stage_df <- elk_yearling_df %>%
  select(year, yearling_mean = mean) %>%
  left_join(elk_ya_df %>% select(year, ya_mean = mean), by = "year") %>%
  left_join(elk_oa_df %>% select(year, oa_mean = mean), by = "year") %>%
  mutate(
    total_mean = yearling_mean + ya_mean + oa_mean,
    prop_yearling = yearling_mean / total_mean,
    prop_ya = ya_mean / total_mean,
    prop_oa = oa_mean / total_mean
  )

elk_first_prop_yearling <- elk_stage_df$prop_yearling[elk_stage_df$year == elk_first_year]
elk_last_prop_yearling <- elk_stage_df$prop_yearling[elk_stage_df$year == elk_last_year]

elk_first_prop_ya <- elk_stage_df$prop_ya[elk_stage_df$year == elk_first_year]
elk_last_prop_ya <- elk_stage_df$prop_ya[elk_stage_df$year == elk_last_year]

elk_first_prop_oa <- elk_stage_df$prop_oa[elk_stage_df$year == elk_first_year]
elk_last_prop_oa <- elk_stage_df$prop_oa[elk_stage_df$year == elk_last_year]

# calculate summary statistics for wolves
wolf_total_df <- wolf_N_summ %>% 
  filter(stage == "Total Wolves") %>% 
  arrange(year) %>%
  drop_na()

wolf_pup_df <- wolf_N_summ %>% 
  filter(stage == "Pups") %>% 
  arrange(year) %>%
  drop_na()

wolf_adult_df <- wolf_N_summ %>% 
  filter(stage == "Adults") %>% 
  arrange(year) %>%
  drop_na()

wolf_first_year <- min(wolf_total_df$year, na.rm=TRUE)
wolf_last_year <- max(wolf_total_df$year, na.rm=TRUE)

wolf_first_total <- wolf_total_df$mean[wolf_total_df$year == wolf_first_year]
wolf_last_total <- wolf_total_df$mean[wolf_total_df$year == wolf_last_year]

wolf_pct_change <- 100 * (wolf_last_total - wolf_first_total) / wolf_first_total
wolf_lambda_geom <- (wolf_last_total / wolf_first_total)^(1 / (wolf_last_year - wolf_first_year))
wolf_annual_pct_change <- 100 * (wolf_lambda_geom - 1)

wolf_stage_df <- wolf_pup_df %>%
  select(year, pup_mean = mean) %>%
  left_join(wolf_adult_df %>% select(year, adult_mean = mean), by = "year") %>%
  mutate(
    total_mean = pup_mean + adult_mean,
    prop_pup = pup_mean / total_mean,
    prop_adult = adult_mean / total_mean
  )

wolf_first_prop_pup <- wolf_stage_df$prop_pup[wolf_stage_df$year == wolf_first_year]
wolf_last_prop_pup <- wolf_stage_df$prop_pup[wolf_stage_df$year == wolf_last_year]
wolf_first_prop_adult <- wolf_stage_df$prop_adult[wolf_stage_df$year == wolf_first_year]
wolf_last_prop_adult <- wolf_stage_df$prop_adult[wolf_stage_df$year == wolf_last_year]

# summarize both elk and wolves
cat(
  " ----------------------------\n",
  "------  ELK SUMMARY  -------\n",
  "----------------------------\n\n",
  
  "---------------------\n",
  " Population growth\n",
  "---------------------\n",
  paste0("Years: ", elk_first_year, " to ", elk_last_year, "\n"),
  paste0("First-year total females: ", round(elk_first_total, 0), "\n"),
  paste0("Last-year total females: ", round(elk_last_total, 0), "\n"),
  paste0("Percent change: ", round(elk_pct_change, 1), "%\n"),
  paste0("Geometric lambda: ", round(elk_lambda_geom, 3), "\n"),
  paste0("Annual percent change: ", round(elk_annual_pct_change, 1), "%\n\n"),
  "---------------------\n",
  "Stage structure\n",
  "---------------------\n",
  paste0(
    "Yearling proportion: ",
    elk_first_year, " = ", round(100 * elk_first_prop_yearling, 1),
    "%; ",
    elk_last_year, " = ", round(100 * elk_last_prop_yearling, 1), "%\n"
  ),
  paste0(
    "Young adult proportion: ",
    elk_first_year, " = ", round(100 * elk_first_prop_ya, 1),
    "%; ",
    elk_last_year, " = ", round(100 * elk_last_prop_ya, 1), "%\n"
  ),
  paste0(
    "Old adult proportion: ",
    elk_first_year, " = ", round(100 * elk_first_prop_oa, 1),
    "%; ",
    elk_last_year, " = ", round(100 * elk_last_prop_oa, 1), "%\n\n"
  ),
  "----------------------------\n",
  "------  WOLF SUMMARY  ------\n",
  "----------------------------\n\n",
  
  "---------------------\n",
  " Population growth\n",
  "---------------------\n",
  paste0("Years: ", wolf_first_year, " to ", wolf_last_year, "\n"),
  paste0("First-year total wolves: ", round(wolf_first_total, 0), "\n"),
  paste0("Last-year total wolves: ", round(wolf_last_total, 0), "\n"),
  paste0("Percent change: ", round(wolf_pct_change, 1), "%\n"),
  paste0("Geometric lambda: ", round(wolf_lambda_geom, 3), "\n"),
  paste0("Annual percent change: ", round(wolf_annual_pct_change, 1), "%\n\n"),
  
  "---------------------\n",
  "Stage structure\n",
  "---------------------\n",
  paste0(
    "Pup proportion: ",
    wolf_first_year, " = ", round(100 * wolf_first_prop_pup, 1),
    "%; ",
    wolf_last_year, " = ", round(100 * wolf_last_prop_pup, 1), "%\n"
  ),
  paste0(
    "Adult proportion: ",
    wolf_first_year, " = ", round(100 * wolf_first_prop_adult, 1),
    "%; ",
    wolf_last_year, " = ", round(100 * wolf_last_prop_adult, 1), "%\n"
  )
)

################################################################################
##########----------------- Coefficient summaries -----------------#############
################################################################################

cat(
  " --------------------------------------------------------\n",
  "COEFFICIENT SUMMARY\n",
  "--------------------------------------------------------\n\n",
  
  "-----------------------------------\n",
  "Strong evidence (95% CI excludes 0)\n",
  "-----------------------------------\n",
  paste0(
    "- ",
    strong_effects %>%
      filter(!grepl("intercept", label, ignore.case = TRUE)) %>%
      transmute(
        text = paste0(
          label, ": median = ", round(median, 3),
          ", 95% CI [", round(ci95_low, 3), ", ", round(ci95_high, 3), "]",
          ", direction = ", direction
        )
      ) %>%
      pull(text),
    collapse = "\n"
  ),
  "\n\n",
  
  "-----------------------------------\n",
  "Moderate evidence (95% CI includes 0, 80% CI excludes 0)\n",
  "-----------------------------------\n",
  paste0(
    "- ",
    moderate_effects %>%
      filter(!grepl("intercept", label, ignore.case = TRUE)) %>%
      transmute(
        text = paste0(
          label, ": median = ", round(median, 3),
          ", 80% CI [", round(ci80_low, 3), ", ", round(ci80_high, 3), "]",
          ", direction = ", direction
        )
      ) %>%
      pull(text),
    collapse = "\n"
  ),
  "\n\n",
  
  "-----------------------------------\n",
  "Weak evidence (80% CI includes 0, 50% CI excludes 0)\n",
  "-----------------------------------\n",
  paste0(
    "- ",
    weak_effects %>%
      filter(!grepl("intercept", label, ignore.case = TRUE)) %>%
      transmute(
        text = paste0(
          label, ": median = ", round(median, 3),
          ", 50% CI [", round(ci50_low, 3), ", ", round(ci50_high, 3), "]",
          ", direction = ", direction
        )
      ) %>%
      pull(text),
    collapse = "\n"
  ),
  "\n\n",
  
  "-----------------------------------\n",
  "Little / no evidence (50% CI includes 0)\n",
  "-----------------------------------\n",
  paste0(
    "- ",
    little_effects %>%
      filter(!grepl("intercept", label, ignore.case = TRUE)) %>%
      transmute(
        text = paste0(
          label, ": median = ", round(median, 3),
          ", 50% CI [", round(ci50_low, 3), ", ", round(ci50_high, 3), "]",
          ", direction = ", direction
        )
      ) %>%
      pull(text),
    collapse = "\n"
  ),
  "\n"
)

cat(
  " --------------------------------------------------------\n",
  "GRIZZLY PROCESS COEFFICIENT SUMMARY\n",
  "--------------------------------------------------------\n\n",
  paste0(
    "Grizzly elk calves effect: median = ",
    round(griz_coef_evidence$median, 3),
    ", 95% CI [",
    round(griz_coef_evidence$ci95_low, 3), ", ",
    round(griz_coef_evidence$ci95_high, 3), "]",
    ", 80% CI [",
    round(griz_coef_evidence$ci80_low, 3), ", ",
    round(griz_coef_evidence$ci80_high, 3), "]",
    ", direction = ",
    griz_coef_evidence$direction,
    ", evidence = ",
    griz_coef_evidence$evidence,
    "\n"
  )
)

cat(
  " --------------------------------------------------------\n",
  "BISON PROCESS COEFFICIENT SUMMARY\n",
  "--------------------------------------------------------\n\n",
  paste0(
    "Bison cull/harvest effect: median = ",
    round(bison_coef_evidence$median, 3),
    ", 95% CI [",
    round(bison_coef_evidence$ci95_low, 3), ", ",
    round(bison_coef_evidence$ci95_high, 3), "]",
    ", 80% CI [",
    round(bison_coef_evidence$ci80_low, 3), ", ",
    round(bison_coef_evidence$ci80_high, 3), "]",
    ", direction = ",
    bison_coef_evidence$direction,
    ", evidence = ",
    bison_coef_evidence$evidence,
    "\n"
  )
)

cat(
  " --------------------------------------------------------\n",
  "COUGAR PROCESS COEFFICIENT SUMMARY\n",
  "--------------------------------------------------------\n\n",
  paste0(
    "Cougar elk effect: median = ",
    round(cougar_coef_evidence$median, 3),
    ", 95% CI [",
    round(cougar_coef_evidence$ci95_low, 3), ", ",
    round(cougar_coef_evidence$ci95_high, 3), "]",
    ", 80% CI [",
    round(cougar_coef_evidence$ci80_low, 3), ", ",
    round(cougar_coef_evidence$ci80_high, 3), "]",
    ", direction = ",
    cougar_coef_evidence$direction,
    ", evidence = ",
    cougar_coef_evidence$evidence,
    "\n"
  )
)


################################################################################
##########------------------- WAIC calculation  -------------------#############
################################################################################

icm_model <- nimbleModel(
  code = icm_code,
  constants = icm_constants,
  data = icm_data,
  inits = waic_inits
)

C_icm_model <- compileNimble(icm_model)

waic_out <- calculateWAIC(
  mcmc = as.matrix(icm_samples_clean),
  model = C_icm_model,
  nburnin = 0,
  thin = 1
)

waic_out

################################################################################
################################################################################
################################################################################
