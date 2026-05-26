### Integrated Community Model (ICM)
### Results exploration
### Last updated: May 20, 2026
### xprockox@gmail.com

################################################################################
#########----------------- Results exploration -----------------################
################################################################################

library(tidyverse)
library(MCMCvis)
library(coda)
library(cowplot)
library(popbio)
library(stringr)

# load model results
load("data/outputs/ICM_parallel_output_2026-05-20.RData")

# import elk data
elk_dat_n <- read.csv("data/elk_abundanceEstimates_stages.csv")
elk_dat_n$n_female <- elk_dat_n$n_cow + (elk_dat_n$n_calf / 2)
elk_dat_n <- elk_dat_n[elk_dat_n$year %in% community_years,]

# select either "NR" or "full" for northern range vs. total YELL poopulation
wolf_range <- "full"

# import wolf data
wolf_pop <- switch(
  wolf_range,
  NR = read.csv("data/wolf_nr_pop.csv"),
  full = read.csv("data/wolf_full_park_pop.csv"),
  stop("wolf_range must be 'NR' or 'full'")
)
wolf_pop <- wolf_pop[wolf_pop$seasonal.year %in% community_years,]

################################################################################
#########---------------------- Settings -----------------------################
################################################################################

ci_prob <- c(0.025, 0.975)

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
  beta0_calfSurv = "Calf intercept",
  beta1_calfSurv_wolfN = "Calf wolf effect",
  beta2_calfSurv_wintPPT = "Calf winter precipitation effect",
  beta3_calfSurv_grizN = "Calf grizzly effect",
  beta4_calfSurv_elkN = "Calf elk density effect",
  beta0_yaSurv = "YA intercept",
  beta1_yaSurv_wolfN = "YA wolf effect",
  beta2_yaSurv_wintPPT = "YA winter precipitation effect",
  beta3_yaSurv_grizN = "YA grizzly effect",
  beta4_yaSurv_elkN = "YA elk density effect",
  beta0_oaSurv = "OA intercept",
  beta1_oaSurv_wolfN = "OA wolf effect",
  beta2_oaSurv_wintPPT = "OA winter precipitation effect",
  beta3_oaSurv_grizN = "OA grizzly effect",
  beta4_oaSurv_elkN = "OA elk density effect",
  beta0_wpupSurv = "Wolf pup intercept",
  beta1_wpupSurv_elkN = "Wolf pup elk effect",
  beta2_wpupSurv_bisonN = "Wolf pup bison effect",
  beta3_wpupSurv_wolfN = "Wolf pup wolf density effect",
  beta0_wadSurv = "Wolf adult intercept",
  beta1_wadSurv_elkN = "Wolf adult elk effect",
  beta2_wadSurv_bisonN = "Wolf adult bison effect",
  beta3_wadSurv_wolfN = "Wolf adult wolf density effect"
)

################################################################################
#########---------------------- Helpers -------------------------################
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

reg_coef_rhat <- get_rhat_issues(MCMCsummary(
  icm_clean,
  params = c(
    "beta0_calfSurv", "beta1_calfSurv_wolfN", "beta2_calfSurv_wintPPT", "beta3_calfSurv_grizN", "beta4_calfSurv_elkN",
    "beta0_yaSurv", "beta1_yaSurv_wolfN", "beta2_yaSurv_wintPPT", "beta3_yaSurv_grizN", "beta4_yaSurv_elkN",
    "beta0_oaSurv", "beta1_oaSurv_wolfN", "beta2_oaSurv_wintPPT", "beta3_oaSurv_grizN", "beta4_oaSurv_elkN",
    "beta0_wpupSurv", "beta1_wpupSurv_elkN", "beta2_wpupSurv_bisonN", "beta3_wpupSurv_wolfN",
    "beta0_wadSurv", "beta1_wadSurv_elkN", "beta2_wadSurv_bisonN", "beta3_wadSurv_wolfN"
  )
))

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

elk_vrate_plot <- plot_vital_rates(
  elk_vrates2,
  title = "Elk posterior time-varying vital rates (95% credible intervals)"
)

wolf_vrate_plot <- plot_vital_rates(
  wolf_vrates2,
  title = "Wolf posterior time-varying vital rates (95% credible intervals)"
)

################################################################################
##########----------------- Plotting data builders -----------------############
################################################################################

elk_surv_pts <- bind_rows(
  elk_vrates2 %>%
    filter(param_base == "Calf survival (s_c)") %>%
    transmute(year, elk_surv = mean, elk_low = low, elk_high = high, stage = "Calf survival"),
  elk_vrates2 %>%
    filter(param_base == "Young Adult survival (s_ya)") %>%
    transmute(year, elk_surv = mean, elk_low = low, elk_high = high, stage = "Young adult survival"),
  elk_vrates2 %>%
    filter(param_base == "Old Adult survival (s_oa)") %>%
    transmute(year, elk_surv = mean, elk_low = low, elk_high = high, stage = "Old adult survival")
)

wolf_surv_pts <- bind_rows(
  wolf_vrates2 %>%
    filter(param_base == "Pup survival (s_p)") %>%
    transmute(year, wolf_surv = mean, wolf_low = low, wolf_high = high, stage = "Pup survival"),
  wolf_vrates2 %>%
    filter(param_base == "Adult survival (s_a)") %>%
    transmute(year, wolf_surv = mean, wolf_low = low, wolf_high = high, stage = "Adult survival")
)

wolf_pts <- summarize_indexed_draws(post_mat, "wolf_N_tot", community_years, "wolf_N_tot") %>%
  rename(wolf_low = low, wolf_high = high)

elk_abund_pts <- summarize_indexed_draws(post_mat, "elk_N_female", community_years, "elk_N_female") %>%
  rename(elk_low_x = low, elk_high_x = high)

wolf_abund_pts <- summarize_indexed_draws(post_mat, "wolf_N_tot", community_years, "wolf_N_tot") %>%
  rename(wolfN_low = low, wolfN_high = high)

plot_df_elk <- elk_surv_pts %>%
  left_join(wolf_pts, by = "year") %>%
  left_join(elk_abund_pts, by = "year") %>%
  left_join(covars %>% select(year, winter_ppt_mm, griz_N), by = "year")

wolf_plot_df <- wolf_surv_pts %>%
  left_join(
    summarize_indexed_draws(post_mat, "elk_N_female", community_years, "elk_N_female") %>%
      rename(elk_low = low, elk_high = high),
    by = "year"
  ) %>%
  left_join(wolf_abund_pts, by = "year") %>%
  left_join(covars %>% select(year, NR_Bison), by = "year")

################################################################################
##########-------------- Dynamic regression plot system ---------------#########
################################################################################

elk_model_specs <- list(
  "Calf survival" = list(
    intercept = "beta0_calfSurv",
    terms = c(
      wolf_N_tot = "beta1_calfSurv_wolfN",
      winter_ppt_mm = "beta2_calfSurv_wintPPT",
      griz_N = "beta3_calfSurv_grizN",
      elk_N_female = "beta4_calfSurv_elkN"
    )
  ),
  "Young adult survival" = list(
    intercept = "beta0_yaSurv",
    terms = c(
      wolf_N_tot = "beta1_yaSurv_wolfN",
      winter_ppt_mm = "beta2_yaSurv_wintPPT",
      griz_N = "beta3_yaSurv_grizN",
      elk_N_female = "beta4_yaSurv_elkN"
    )
  ),
  "Old adult survival" = list(
    intercept = "beta0_oaSurv",
    terms = c(
      wolf_N_tot = "beta1_oaSurv_wolfN",
      winter_ppt_mm = "beta2_oaSurv_wintPPT",
      griz_N = "beta3_oaSurv_grizN",
      elk_N_female = "beta4_oaSurv_elkN"
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
  winter_ppt_mm = list(
    raw_col = "winter_ppt_mm",
    raw_df = covars,
    std_fun = function(x) (x - mean(covars$winter_ppt_mm, na.rm = TRUE)) / sd(covars$winter_ppt_mm, na.rm = TRUE),
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
    label = "Winter precipitation"
  ),
  griz_N = list(
    raw_col = "griz_N",
    raw_df = covars,
    std_fun = function(x) (x - mean(covars$griz_N, na.rm = TRUE)) / sd(covars$griz_N, na.rm = TRUE),
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
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
  NR_Bison = list(
    raw_col = "NR_Bison",
    raw_df = covars,
    std_fun = function(x) (x - mean(covars$NR_Bison, na.rm = TRUE)) / sd(covars$NR_Bison, na.rm = TRUE),
    held_value = 0,
    xmin = NULL,
    xmax = NULL,
    label = "Bison abundance"
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
    geom_errorbar(aes(ymin = .data[[ymin_var]], ymax = .data[[ymax_var]]), width = 0) +
    geom_point(shape = 21, fill = "gold", color = "black", size = 2) +
    facet_wrap(~stage, scales = "free_y") +
    theme_classic() +
    labs(
      x = x_label,
      y = y_label,
      title = title,
      subtitle = subtitle
    )
  
  if (!is.null(xlow_var) && !is.null(xhigh_var)) {
    p <- p + geom_errorbarh(aes(xmin = .data[[xlow_var]], xmax = .data[[xhigh_var]]), height = 0)
  }
  
  p
}

################################################################################
############------------ Reproducing regression plots --------------############
################################################################################

line_df_wolf_elk <- make_effect_curve(post_mat, elk_model_specs, "wolf_N_tot", predictor_specs)
line_df_ppt_elk <- make_effect_curve(post_mat, elk_model_specs, "winter_ppt_mm", predictor_specs)
line_df_griz_elk <- make_effect_curve(post_mat, elk_model_specs, "griz_N", predictor_specs)
line_df_elkN_elk <- make_effect_curve(post_mat, elk_model_specs, "elk_N_female", predictor_specs)

line_df_elk_wolf <- make_effect_curve(post_mat, wolf_model_specs, "elk_N_female", predictor_specs)
line_df_bison_wolf <- make_effect_curve(post_mat, wolf_model_specs, "NR_Bison", predictor_specs)
line_df_wolf_wolf <- make_effect_curve(post_mat, wolf_model_specs, "wolf_N_tot", predictor_specs)

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
  subtitle = "Winter precipitation, grizzly abundance, and elk abundance held constant",
  fill_color = "#6F263D",
  line_color = "#6F263D"
)

ppt_plot_elk <- plot_effect(
  plot_df = plot_df_elk,
  curve_df = line_df_ppt_elk,
  x_var = "winter_ppt_mm",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  x_label = "Winter precipitation",
  y_label = "Elk survival",
  title = "Estimated effect of winter precipitation on elk survival",
  subtitle = "Wolf abundance, grizzly abundance, and elk abundance held constant",
  fill_color = "#236192",
  line_color = "#236192"
)

griz_plot_elk <- plot_effect(
  plot_df = plot_df_elk,
  curve_df = line_df_griz_elk,
  x_var = "griz_N",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  x_label = "Grizzly abundance",
  y_label = "Elk survival",
  title = "Estimated effect of grizzly abundance on elk survival",
  subtitle = "Wolf abundance, winter precipitation, and elk abundance held constant",
  fill_color = "#4B7F52",
  line_color = "#4B7F52"
)

elkN_plot_elk <- plot_effect(
  plot_df = plot_df_elk,
  curve_df = line_df_elkN_elk,
  x_var = "elk_N_female",
  y_var = "elk_surv",
  ymin_var = "elk_low",
  ymax_var = "elk_high",
  xlow_var = "elk_low_x",
  xhigh_var = "elk_high_x",
  x_label = "Elk abundance",
  y_label = "Elk survival",
  title = "Estimated effect of elk abundance on elk survival",
  subtitle = "Wolf abundance, winter precipitation, and grizzly abundance held constant",
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
  subtitle = "Bison abundance and wolf abundance held constant",
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
  subtitle = "Elk abundance and wolf abundance held constant",
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
  title = "Estimated effect of wolf abundance on wolf survival",
  subtitle = "Elk abundance and bison abundance held constant",
  fill_color = "#4B7F52",
  line_color = "#4B7F52"
)

################################################################################
##########----------------- Coefficient density plots -----------------#########
################################################################################

make_coef_density_plot <- function(post_mat, coef_names, label_map, keep, fill_color, title, ncol) {
  df <- as_tibble(post_mat[, coef_names, drop = FALSE]) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(parameter = recode(parameter, !!!label_map))
  
  ggplot(df %>% filter(parameter %in% keep), aes(x = value)) +
    geom_density(fill = fill_color, alpha = 0.45) +
    geom_vline(xintercept = 0, linetype = 2) +
    facet_wrap(~parameter, scales = "free", ncol = ncol) +
    theme_classic() +
    labs(
      x = "Posterior value",
      y = "Density",
      title = title
    )
}

coef_plot <- make_coef_density_plot(
  post_mat = post_mat,
  coef_names = c(
    "beta0_calfSurv", "beta1_calfSurv_wolfN", "beta2_calfSurv_wintPPT", "beta3_calfSurv_grizN", "beta4_calfSurv_elkN",
    "beta0_yaSurv", "beta1_yaSurv_wolfN", "beta2_yaSurv_wintPPT", "beta3_yaSurv_grizN", "beta4_yaSurv_elkN",
    "beta0_oaSurv", "beta1_oaSurv_wolfN", "beta2_oaSurv_wintPPT", "beta3_oaSurv_grizN", "beta4_oaSurv_elkN"
  ),
  label_map = pretty_coef_labels,
  keep = c("Calf intercept", "Calf grizzly effect", "OA intercept", "OA grizzly effect", "YA intercept", "YA grizzly effect"),
  fill_color = "#236192",
  title = "Posterior distributions of elk regression coefficients",
  ncol = 6
)

wolf_coef_plot <- make_coef_density_plot(
  post_mat = post_mat,
  coef_names = c(
    "beta0_wpupSurv", "beta1_wpupSurv_elkN", "beta2_wpupSurv_bisonN", "beta3_wpupSurv_wolfN",
    "beta0_wadSurv", "beta1_wadSurv_elkN", "beta2_wadSurv_bisonN", "beta3_wadSurv_wolfN"
  ),
  label_map = pretty_coef_labels,
  keep = c("Wolf pup intercept", "Wolf pup bison effect", "Wolf adult intercept", "Wolf adult bison effect"),
  fill_color = "#6F263D",
  title = "Posterior distributions of wolf regression coefficients",
  ncol = 4
)

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
##########----------------- Wolf change summaries -----------------#############
################################################################################

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

cat("Wolf abundance summary\n")
cat("----------------------\n")
cat("Years:", wolf_first_year, "to", wolf_last_year, "\n")
cat("First-year total wolves:", round(wolf_first_total, 0), "\n")
cat("Last-year total wolves:", round(wolf_last_total, 0), "\n")
cat("Percent change:", round(wolf_pct_change, 1), "%\n")
cat("Geometric lambda:", round(wolf_lambda_geom, 3), "\n")
cat("Annual percent change:", round(wolf_annual_pct_change, 1), "%\n\n")

cat("Wolf stage structure summary\n")
cat("----------------------------\n")
cat("Pup proportion:", wolf_first_year, "=", round(100 * wolf_first_prop_pup, 1),
    "%;", wolf_last_year, "=", round(100 * wolf_last_prop_pup, 1), "%\n")
cat("Adult proportion:", wolf_first_year, "=", round(100 * wolf_first_prop_adult, 1),
    "%;", wolf_last_year, "=", round(100 * wolf_last_prop_adult, 1), "%\n")

################################################################################
################################ Final panels ##################################
################################################################################

elk_panel <- plot_grid(
  wolf_plot_elk,
  ppt_plot_elk,
  griz_plot_elk,
  elkN_plot_elk,
  coef_plot,
  ncol = 1,
  rel_heights = c(1.1, 1.1, 1.1, 1.1, 1),
  labels = c("A", "B", "C", "D", "E")
)

wolf_panel <- plot_grid(
  wolf_main_plot,
  bison_plot_wolf,
  wolfN_plot_wolf,
  wolf_coef_plot,
  ncol = 1,
  rel_heights = c(1.2, 1.2, 1.2, 1),
  labels = c("F", "G", "H", "I")
)

all_regression_plots <- plot_grid(
  elk_panel,
  wolf_panel,
  ncol = 2,
  rel_widths = c(1, 1)
)

all_regression_plots