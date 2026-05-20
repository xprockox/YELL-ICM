### Integrated Community Model (ICM)
### Results exploration
### Last updated: May 20, 2026
### xprockox@gmail.com

################################################################################
############---------- Packages and user settings ---------#####################
################################################################################
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(MCMCvis)
library(ggplot2)
library(nimble)
library(coda)
library(stringr)
library(cowplot)
library(popbio)

################################################################################
##########---------------- Load model results ---------------###################
################################################################################

load('data/outputs/ICM_parallel_output_2026-05-19.RData')

################################################################################
###########--------------------- Results ----------------------#################
################################################################################

# combine chains into one matrix so we can identify bad parameters
post_mat_raw <- as.matrix(icm_samples)

# identify parameters with any NA / NaN / Inf
bad_cols <- which(
  apply(post_mat_raw, 2, function(x) any(!is.finite(x)))
)

bad_params <- colnames(post_mat_raw)[bad_cols]

bad_summary <- data.frame(
  parameter = bad_params,
  n_bad = sapply(bad_cols, function(i) sum(!is.finite(post_mat_raw[, i]))),
  n_total = nrow(post_mat_raw)
)

bad_summary

# keep only parameters that are finite in all chains
good_cols <- colnames(post_mat_raw)[-bad_cols]

chain_samples_clean <- lapply(chain_samples, function(ch) {
  ch_mat <- as.matrix(ch)
  mcmc(ch_mat[, good_cols, drop = FALSE])
})

icm_clean <- mcmc.list(chain_samples_clean)

# summary table with Rhat and n.eff retained
icm_summary <- MCMCsummary(icm_clean)
post_sum <- as.data.frame(icm_summary)
post_sum$param <- rownames(post_sum)

# pooled posterior matrix for plotting / prediction
post_mat <- as.matrix(icm_clean)

# quick check
round(icm_summary, 2)

################################################################################
###########------------ Checking model convergence ------------#################
################################################################################

# elk abundances
elk_N_summ <- MCMCsummary(
  icm_clean,
  params = c("elk_N_1y", "elk_N_ya", "elk_N_oa", "elk_N_female")
)

bad_elk_N_rhats <- data.frame(
  param = rownames(elk_N_summ),
  Rhat = elk_N_summ[, "Rhat"]
)

bad_elk_N_rhats <- bad_elk_N_rhats[bad_elk_N_rhats$Rhat > 1.1, ]
bad_elk_N_rhats

# elk vital rates
elk_vrates <- MCMCsummary(
  icm_clean,
  params = c("elk_s_c", "elk_s_ya", "elk_s_oa", "elk_p_13", "elk_f_ya", "elk_f_oa")
) 

bad_elk_vrate_rhats <- data.frame(
  param = rownames(elk_vrates),
  Rhat = elk_vrates[, "Rhat"]
)

bad_elk_vrate_rhats <- bad_elk_vrate_rhats[bad_elk_vrate_rhats$Rhat > 1.1, ]
bad_elk_vrate_rhats

# wolf abundances
wolf_N_summ <- MCMCsummary(
  icm_clean,
  params = c("wolf_N_p", "wolf_N_a", "wolf_N_tot")
)

bad_wolf_N_rhats <- data.frame(
  param = rownames(wolf_N_summ),
  Rhat = wolf_N_summ[, "Rhat"]
)

bad_wolf_N_rhats <- bad_wolf_N_rhats[bad_wolf_N_rhats$Rhat > 1.1, ]
bad_wolf_N_rhats

# wolf vital rates
wolf_vrates <- MCMCsummary(
  icm_clean,
  params = c("wolf_s_p", "wolf_s_a", "wolf_f")
) 

bad_wolf_vrate_rhats <- data.frame(
  param = rownames(wolf_vrates),
  Rhat = wolf_vrates[, "Rhat"]
)

bad_wolf_vrate_rhats <- bad_wolf_vrate_rhats[bad_wolf_vrate_rhats$Rhat > 1.1, ]
bad_wolf_vrate_rhats

# regression coefficients
reg_coefs <- MCMCsummary(
  icm_clean,
  params = c(
    # elk reg coefs
    'beta0_calfSurv', 'beta1_calfSurv_wolfN', 'beta2_calfSurv_wintPPT', 'beta3_calfSurv_grizN', 'beta4_calfSurv_elkN',
    'beta0_yaSurv', 'beta1_yaSurv_wolfN', 'beta2_yaSurv_wintPPT', 'beta3_yaSurv_grizN', 'beta4_yaSurv_elkN',
    'beta0_oaSurv', 'beta1_oaSurv_wolfN', 'beta2_oaSurv_wintPPT', 'beta3_oaSurv_grizN', 'beta4_oaSurv_elkN',
    # wolf reg coefs
    'beta0_wpupSurv', 'beta1_wpupSurv_elkN', 'beta2_wpupSurv_bisonN', 'beta3_wpupSurv_wolfN',
    'beta0_wadSurv', 'beta1_wadSurv_elkN', 'beta2_wadSurv_bisonN', 'beta3_wadSurv_wolfN'
  )
)

bad_reg_coefs <- data.frame(
  param = rownames(reg_coefs),
  Rhat = reg_coefs[, "Rhat"]
)

bad_reg_coefs <- bad_reg_coefs[bad_reg_coefs$Rhat > 1.1, ]
bad_reg_coefs

################################################################################
###########--------- Plotting elk abundance posteriors --------#################
################################################################################

elk_N_summ <- elk_N_summ %>%
  rownames_to_column("param") %>%
  mutate(
    stage = str_extract(param, "^elk_N_[a-zA-Z0-9_]+"),
    t = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
  ) %>%
  select(stage, t, mean = mean, low = `2.5%`, high = `97.5%`) %>%
  arrange(stage, t)

elk_N_summ <- elk_N_summ %>%
  mutate(year = community_years[t])

elk_N_summ <- elk_N_summ %>%
  mutate(stage = recode(
    stage,
    "elk_N_1y" = "Yearling",
    "elk_N_ya" = "Young Adult",
    "elk_N_oa" = "Old Adult",
    "elk_N_female" = "Total Females"
  ))

elk_dat_long <- elk_dat_n %>%
  pivot_longer(
    cols = -c(year),
    names_to = "stage",
    values_to = "value"
  ) %>%
  mutate(stage = recode(
    stage,
    "n_calf" = "Yearling",
    "n_cow_youngadult" = "Young Adult",
    "n_cow_oldadult" = "Old Adult",
    "n_female" = "Total Females"
  ))

elk_dat_long$stage <- factor(elk_dat_long$stage, levels = unique(elk_N_summ$stage))

elk_dat_long <- elk_dat_long[
  elk_dat_long$stage %in% c("Yearling", "Young Adult", "Old Adult", "Total Females"),
]

elk_validation_plot <- ggplot(elk_N_summ, aes(x = year, y = mean, group = stage)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = stage), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(
    data = elk_dat_long[elk_dat_long$stage == "Total Females", ],
    aes(y = value),
    color = "red",
    size = 2
  ) +
  geom_line(
    data = elk_dat_long[elk_dat_long$stage == "Total Females", ],
    aes(y = value),
    color = "red",
    linetype = 2
  ) +
  facet_wrap(~stage, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Abundance",
    title = "Elk posterior population estimates with validation data",
    subtitle = "Ribbon = 95% credible interval, Line = posterior mean, Red = observed"
  ) +
  theme(legend.position = "none")

elk_validation_plot

################################################################################
###########--------- Plotting elk vital rate posteriors -------#################
################################################################################

elk_vrates <- elk_vrates %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  rename(mean = mean, low = `2.5%`, high = `97.5%`)

elk_vrates2 <- elk_vrates %>%
  mutate(
    year_index = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    rate = str_extract(param, "^[^\\[]+")
  )

elk_vrates2$rate <- factor(
  elk_vrates2$rate,
  levels = c("elk_s_c", "elk_s_ya", "elk_s_oa", "elk_p_13", "elk_f_ya", "elk_f_oa"),
  labels = c(
    "Calf survival (s_c)",
    "Young Adult survival (s_ya)",
    "Old Adult survival (s_oa)",
    "Young→Old transition (p_13)",
    "Fecundity (young) (f_ya)",
    "Fecundity (old) (f_oa)"
  )
)

elk_vrates2$year <- NA
elk_vrates2$year[elk_vrates2$rate %in% c(
  "Calf survival (s_c)",
  "Young Adult survival (s_ya)",
  "Old Adult survival (s_oa)",
  "Young→Old transition (p_13)"
)] <- community_years[elk_vrates2$year_index[
  elk_vrates2$rate %in% c(
    "Calf survival (s_c)",
    "Young Adult survival (s_ya)",
    "Old Adult survival (s_oa)",
    "Young→Old transition (p_13)"
  )
]]

elk_vrates2$year[elk_vrates2$rate %in% c(
  "Fecundity (young) (f_ya)",
  "Fecundity (old) (f_oa)"
)] <- community_years[-length(community_years)][elk_vrates2$year_index[
  elk_vrates2$rate %in% c(
    "Fecundity (young) (f_ya)",
    "Fecundity (old) (f_oa)"
  )
]]

elk_vrate_plot <- ggplot(elk_vrates2, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = rate), alpha = 0.2) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~rate, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Estimated value",
    title = "Elk posterior time-varying vital rates (95% credible intervals)"
  ) +
  theme(legend.position = "none")

elk_vrate_plot

################################################################################
###########--------- Plotting wolf abundance posteriors -------#################
################################################################################

wolf_N_summ <- wolf_N_summ %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  mutate(
    stage = str_extract(param, "^wolf_N_[a-zA-Z0-9_]+"),
    t = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
  ) %>%
  select(stage, t, mean = mean, low = `2.5%`, high = `97.5%`) %>%
  arrange(stage, t) %>%
  mutate(year = community_years[t]) %>%
  mutate(stage = recode(
    stage,
    "wolf_N_p" = "Pups",
    "wolf_N_a" = "Adults",
    "wolf_N_tot" = "Total Wolves"
  ))

wolf_dat_long <- wolf_pop %>%
  select(seasonal.year, dec_pups, dec_adults, total_abundance) %>%
  rename(year = seasonal.year) %>%
  pivot_longer(
    cols = -c(year),
    names_to = "stage",
    values_to = "value"
  ) %>%
  mutate(stage = recode(
    stage,
    "dec_pups" = "Pups",
    "dec_adults" = "Adults",
    "total_abundance" = "Total Wolves"
  ))

wolf_dat_long$stage <- factor(wolf_dat_long$stage, levels = unique(wolf_N_summ$stage))

wolf_validation_plot <- ggplot(wolf_N_summ, aes(x = year, y = mean, group = stage)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = stage), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(
    data = wolf_dat_long,
    aes(y = value),
    color = "red",
    size = 2
  ) +
  geom_line(
    data = wolf_dat_long,
    aes(y = value),
    color = "red",
    linetype = 2
  ) +
  facet_wrap(~stage, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Abundance",
    title = "Wolf posterior population estimates with validation data",
    subtitle = "Ribbon = 95% credible interval, Line = posterior mean, Red = observed"
  ) +
  theme(legend.position = "none")

wolf_validation_plot

################################################################################
##########--------- Plotting wolf vital rate posteriors -------#################
################################################################################

wolf_vrates <- wolf_vrates %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  rename(mean = mean, low = `2.5%`, high = `97.5%`)

wolf_vrates2 <- wolf_vrates %>%
  mutate(
    year_index = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    rate = str_extract(param, "^[^\\[]+")
  )

wolf_vrates2$rate <- factor(
  wolf_vrates2$rate,
  levels = c("wolf_s_p", "wolf_s_a", "wolf_f"),
  labels = c(
    "Pup survival (s_p)",
    "Adult survival (s_a)",
    "Fecundity (f)"
  )
)

wolf_vrates2$year <- NA
wolf_vrates2$year[wolf_vrates2$rate %in% c(
  "Pup survival (s_p)",
  "Adult survival (s_a)"
)] <- community_years[wolf_vrates2$year_index[
  wolf_vrates2$rate %in% c(
    "Pup survival (s_p)",
    "Adult survival (s_a)"
  )
]]

wolf_vrates2$year[wolf_vrates2$rate == "Fecundity (f)"] <- community_years[-length(community_years)][
  wolf_vrates2$year_index[wolf_vrates2$rate == "Fecundity (f)"]
]

wolf_vrate_plot <- ggplot(wolf_vrates2, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = rate), alpha = 0.2) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~rate, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Estimated value",
    title = "Wolf posterior time-varying vital rates (95% credible intervals)"
  ) +
  theme(legend.position = "none")

wolf_vrate_plot

################################################################################
#########--------- Rebuild plotting dataframe for regressions ------#########
################################################################################

# yearly wolf abundance summaries
wolf_pts <- wolf_N_summ %>%
  filter(stage == "Total Wolves") %>%
  transmute(
    year,
    wolf_N_tot = mean,
    wolf_low = low,
    wolf_high = high
  )

# yearly elk abundance summaries
elk_abund_pts <- elk_N_summ %>%
  filter(stage == "Total Females") %>%
  transmute(
    year,
    elk_N_female = mean,
    elk_low_x = low,
    elk_high_x = high
  )

# yearly elk survival summaries
elk_calf_pts <- elk_vrates2 %>%
  filter(rate == "Calf survival (s_c)") %>%
  transmute(
    year,
    elk_surv = mean,
    elk_low = low,
    elk_high = high,
    stage = "Calf survival"
  )

elk_ya_pts <- elk_vrates2 %>%
  filter(rate == "Young Adult survival (s_ya)") %>%
  transmute(
    year,
    elk_surv = mean,
    elk_low = low,
    elk_high = high,
    stage = "Young adult survival"
  )

elk_oa_pts <- elk_vrates2 %>%
  filter(rate == "Old Adult survival (s_oa)") %>%
  transmute(
    year,
    elk_surv = mean,
    elk_low = low,
    elk_high = high,
    stage = "Old adult survival"
  )

elk_surv_pts <- bind_rows(elk_calf_pts, elk_ya_pts, elk_oa_pts)

# final plotting dataframe
plot_df_elk <- elk_surv_pts %>%
  left_join(wolf_pts, by = "year") %>%
  left_join(elk_abund_pts, by = "year") %>%
  left_join(
    covars %>%
      select(year, winter_ppt_mm, griz_N),
    by = "year"
  )

# quick check
summary(plot_df_elk$elk_N_female)
sum(is.na(plot_df_elk$elk_N_female))

################################################################################
########---------------- Wolf abundance effect on elk ----------------##########
################################################################################

# x grid on raw wolf abundance scale
x_grid_wolf <- seq(
  min(plot_df_elk$wolf_N_tot, na.rm = TRUE),
  max(plot_df_elk$wolf_N_tot, na.rm = TRUE),
  length.out = 200
)

# standardize internally because model used standardized wolf abundance
x_grid_wolf_std <- (x_grid_wolf - icm_constants$wolf_tot_mean) / icm_constants$wolf_tot_sd

# hold other predictors constant at their mean standardized values
winterPPT_mean_std <- mean(covars_std$winter_ppt_mm, na.rm = TRUE)
grizN_mean_std <- mean(covars_std$griz_N, na.rm = TRUE)
elkN_mean_std <- 0

# fitted values for calf survival
pred_calf_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_calfSurv"] +
      post_mat[, "beta1_calfSurv_wolfN"] * x +
      post_mat[, "beta2_calfSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_calfSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_calfSurv_elkN"] * elkN_mean_std
  )
})

line_calf_wolf <- data.frame(
  x = x_grid_wolf,
  elk_surv = apply(pred_calf_wolf, 2, mean),
  elk_low = apply(pred_calf_wolf, 2, quantile, probs = 0.025, na.rm = TRUE),
  elk_high = apply(pred_calf_wolf, 2, quantile, probs = 0.975, na.rm = TRUE),
  stage = "Calf survival"
)

# fitted values for young adult survival
pred_ya_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_yaSurv"] +
      post_mat[, "beta1_yaSurv_wolfN"] * x +
      post_mat[, "beta2_yaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_yaSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_yaSurv_elkN"] * elkN_mean_std
  )
})

line_ya_wolf <- data.frame(
  x = x_grid_wolf,
  elk_surv = apply(pred_ya_wolf, 2, mean),
  elk_low = apply(pred_ya_wolf, 2, quantile, probs = 0.025, na.rm = TRUE),
  elk_high = apply(pred_ya_wolf, 2, quantile, probs = 0.975, na.rm = TRUE),
  stage = "Young adult survival"
)

# fitted values for old adult survival
pred_oa_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_oaSurv"] +
      post_mat[, "beta1_oaSurv_wolfN"] * x +
      post_mat[, "beta2_oaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_oaSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_oaSurv_elkN"] * elkN_mean_std
  )
})

line_oa_wolf <- data.frame(
  x = x_grid_wolf,
  elk_surv = apply(pred_oa_wolf, 2, mean),
  elk_low = apply(pred_oa_wolf, 2, quantile, probs = 0.025, na.rm = TRUE),
  elk_high = apply(pred_oa_wolf, 2, quantile, probs = 0.975, na.rm = TRUE),
  stage = "Old adult survival"
)

line_df_wolf_elk <- bind_rows(line_calf_wolf, line_ya_wolf, line_oa_wolf)

wolf_plot_elk <- ggplot(plot_df_elk, aes(x = wolf_N_tot, y = elk_surv)) +
  geom_ribbon(
    data = line_df_wolf_elk,
    aes(x = x, ymin = elk_low, ymax = elk_high),
    inherit.aes = FALSE,
    fill = "#6F263D",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df_wolf_elk,
    aes(x = x, y = elk_surv),
    inherit.aes = FALSE,
    color = "#6F263D",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = elk_low, ymax = elk_high), width = 0) +
  geom_errorbarh(aes(xmin = wolf_low, xmax = wolf_high), height = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Wolf abundance",
    y = "Elk survival",
    title = "Estimated effect of wolf abundance on elk survival",
    subtitle = "Winter precipitation, grizzly abundance, and elk abundance held constant"
  )

wolf_plot_elk

################################################################################
#######------------- Winter precipitation effect on elk ----------------########
################################################################################

# x grid on raw winter precipitation scale
x_grid_ppt_raw <- seq(
  min(covars$winter_ppt_mm, na.rm = TRUE),
  max(covars$winter_ppt_mm, na.rm = TRUE),
  length.out = 200
)

# convert raw ppt grid to standardized values used in the model
ppt_mean_raw <- mean(covars$winter_ppt_mm, na.rm = TRUE)
ppt_sd_raw <- sd(covars$winter_ppt_mm, na.rm = TRUE)
x_grid_ppt_std <- (x_grid_ppt_raw - ppt_mean_raw) / ppt_sd_raw

# hold other predictors constant at mean standardized values
wolf_mean_std <- 0
grizN_mean_std <- 0
elkN_mean_std <- 0

pred_calf_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_calfSurv"] +
      post_mat[, "beta1_calfSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_calfSurv_wintPPT"] * x +
      post_mat[, "beta3_calfSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_calfSurv_elkN"] * elkN_mean_std
  )
})

pred_ya_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_yaSurv"] +
      post_mat[, "beta1_yaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_yaSurv_wintPPT"] * x +
      post_mat[, "beta3_yaSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_yaSurv_elkN"] * elkN_mean_std
  )
})

pred_oa_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_oaSurv"] +
      post_mat[, "beta1_oaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_oaSurv_wintPPT"] * x +
      post_mat[, "beta3_oaSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_oaSurv_elkN"] * elkN_mean_std
  )
})

line_df_ppt_elk <- bind_rows(
  data.frame(
    x = x_grid_ppt_raw,
    elk_surv = apply(pred_calf_ppt, 2, mean),
    elk_low = apply(pred_calf_ppt, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_calf_ppt, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Calf survival"
  ),
  data.frame(
    x = x_grid_ppt_raw,
    elk_surv = apply(pred_ya_ppt, 2, mean),
    elk_low = apply(pred_ya_ppt, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_ya_ppt, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Young adult survival"
  ),
  data.frame(
    x = x_grid_ppt_raw,
    elk_surv = apply(pred_oa_ppt, 2, mean),
    elk_low = apply(pred_oa_ppt, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_oa_ppt, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Old adult survival"
  )
)

ppt_plot_elk <- ggplot(plot_df_elk, aes(x = winter_ppt_mm, y = elk_surv)) +
  geom_ribbon(
    data = line_df_ppt_elk,
    aes(x = x, ymin = elk_low, ymax = elk_high),
    inherit.aes = FALSE,
    fill = "#236192",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df_ppt_elk,
    aes(x = x, y = elk_surv),
    inherit.aes = FALSE,
    color = "#236192",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = elk_low, ymax = elk_high), width = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Winter precipitation",
    y = "Elk survival",
    title = "Estimated effect of winter precipitation on elk survival",
    subtitle = "Wolf abundance, grizzly abundance, and elk abundance held constant"
  )

ppt_plot_elk

################################################################################
########------------- Grizzly abundance effect on elk -----------------##########
################################################################################

# x grid on raw grizzly abundance scale
x_grid_griz_raw <- seq(
  min(plot_df_elk$griz_N, na.rm = TRUE),
  max(plot_df_elk$griz_N, na.rm = TRUE),
  length.out = 200
)

# convert raw grizzly grid to standardized values used in the model
griz_mean_raw <- mean(covars$griz_N, na.rm = TRUE)
griz_sd_raw <- sd(covars$griz_N, na.rm = TRUE)
x_grid_griz_std <- (x_grid_griz_raw - griz_mean_raw) / griz_sd_raw

# hold other predictors constant at mean standardized values
wolf_mean_std <- 0
winterPPT_mean_std <- 0
elkN_mean_std <- 0

pred_calf_griz <- sapply(x_grid_griz_std, function(x) {
  plogis(
    post_mat[, "beta0_calfSurv"] +
      post_mat[, "beta1_calfSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_calfSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_calfSurv_grizN"] * x +
      post_mat[, "beta4_calfSurv_elkN"] * elkN_mean_std
  )
})

pred_ya_griz <- sapply(x_grid_griz_std, function(x) {
  plogis(
    post_mat[, "beta0_yaSurv"] +
      post_mat[, "beta1_yaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_yaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_yaSurv_grizN"] * x +
      post_mat[, "beta4_yaSurv_elkN"] * elkN_mean_std
  )
})

pred_oa_griz <- sapply(x_grid_griz_std, function(x) {
  plogis(
    post_mat[, "beta0_oaSurv"] +
      post_mat[, "beta1_oaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_oaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_oaSurv_grizN"] * x +
      post_mat[, "beta4_oaSurv_elkN"] * elkN_mean_std
  )
})

line_df_griz_elk <- bind_rows(
  data.frame(
    x = x_grid_griz_raw,
    elk_surv = apply(pred_calf_griz, 2, mean),
    elk_low = apply(pred_calf_griz, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_calf_griz, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Calf survival"
  ),
  data.frame(
    x = x_grid_griz_raw,
    elk_surv = apply(pred_ya_griz, 2, mean),
    elk_low = apply(pred_ya_griz, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_ya_griz, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Young adult survival"
  ),
  data.frame(
    x = x_grid_griz_raw,
    elk_surv = apply(pred_oa_griz, 2, mean),
    elk_low = apply(pred_oa_griz, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_oa_griz, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Old adult survival"
  )
)

griz_plot_elk <- ggplot(plot_df_elk, aes(x = griz_N, y = elk_surv)) +
  geom_ribbon(
    data = line_df_griz_elk,
    aes(x = x, ymin = elk_low, ymax = elk_high),
    inherit.aes = FALSE,
    fill = "#4B7F52",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df_griz_elk,
    aes(x = x, y = elk_surv),
    inherit.aes = FALSE,
    color = "#4B7F52",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = elk_low, ymax = elk_high), width = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Grizzly abundance",
    y = "Elk survival",
    title = "Estimated effect of grizzly abundance on elk survival",
    subtitle = "Wolf abundance, winter precipitation, and elk abundance held constant"
  )

griz_plot_elk

################################################################################
########------------- Elk abundance effect on elk survival -------------#########
################################################################################

# x grid on raw elk abundance scale
x_grid_elk_raw <- seq(
  min(plot_df_elk$elk_N_female, na.rm = TRUE),
  max(plot_df_elk$elk_N_female, na.rm = TRUE),
  length.out = 200
)

# convert raw elk abundance grid to standardized values used in the model
elk_mean_raw <- icm_constants$elk_N_female_mean
elk_sd_raw <- icm_constants$elk_N_female_sd
x_grid_elk_std <- (x_grid_elk_raw - elk_mean_raw) / elk_sd_raw

# hold other predictors constant at mean standardized values
wolf_mean_std <- 0
winterPPT_mean_std <- 0
grizN_mean_std <- 0

# fitted values for calf survival
pred_calf_elkN <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_calfSurv"] +
      post_mat[, "beta1_calfSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_calfSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_calfSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_calfSurv_elkN"] * x
  )
})

# fitted values for young adult survival
pred_ya_elkN <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_yaSurv"] +
      post_mat[, "beta1_yaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_yaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_yaSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_yaSurv_elkN"] * x
  )
})

# fitted values for old adult survival
pred_oa_elkN <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_oaSurv"] +
      post_mat[, "beta1_oaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_oaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_oaSurv_grizN"] * grizN_mean_std +
      post_mat[, "beta4_oaSurv_elkN"] * x
  )
})

line_df_elkN_elk <- bind_rows(
  data.frame(
    x = x_grid_elk_raw,
    elk_surv = apply(pred_calf_elkN, 2, mean),
    elk_low = apply(pred_calf_elkN, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_calf_elkN, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Calf survival"
  ),
  data.frame(
    x = x_grid_elk_raw,
    elk_surv = apply(pred_ya_elkN, 2, mean),
    elk_low = apply(pred_ya_elkN, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_ya_elkN, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Young adult survival"
  ),
  data.frame(
    x = x_grid_elk_raw,
    elk_surv = apply(pred_oa_elkN, 2, mean),
    elk_low = apply(pred_oa_elkN, 2, quantile, probs = 0.025, na.rm = TRUE),
    elk_high = apply(pred_oa_elkN, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Old adult survival"
  )
)

elkN_plot_elk <- ggplot(plot_df_elk, aes(x = elk_N_female, y = elk_surv)) +
  geom_ribbon(
    data = line_df_elkN_elk,
    aes(x = x, ymin = elk_low, ymax = elk_high),
    inherit.aes = FALSE,
    fill = "#8C510A",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df_elkN_elk,
    aes(x = x, y = elk_surv),
    inherit.aes = FALSE,
    color = "#8C510A",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = elk_low, ymax = elk_high), width = 0) +
  geom_errorbarh(aes(xmin = elk_low_x, xmax = elk_high_x), height = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Elk abundance",
    y = "Elk survival",
    title = "Estimated effect of elk abundance on elk survival",
    subtitle = "Wolf abundance, winter precipitation, and grizzly abundance held constant"
  )

elkN_plot_elk

################################################################################
##########---------- Wolf regression plots + coefficient densities ----##########
################################################################################

# yearly elk abundance summaries for wolf regressions
elk_abund_pts <- data.frame(
  year = elk_dat_n$year,
  elk_N_female = apply(post_mat[, grep("^elk_N_female\\[", colnames(post_mat))], 2, mean),
  elk_low = apply(post_mat[, grep("^elk_N_female\\[", colnames(post_mat))], 2, quantile, probs = 0.025, na.rm = TRUE),
  elk_high = apply(post_mat[, grep("^elk_N_female\\[", colnames(post_mat))], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# yearly wolf abundance summaries
wolf_abund_pts <- data.frame(
  year = elk_dat_n$year,
  wolf_N_tot = apply(post_mat[, grep("^wolf_N_tot\\[", colnames(post_mat))], 2, mean),
  wolfN_low = apply(post_mat[, grep("^wolf_N_tot\\[", colnames(post_mat))], 2, quantile, probs = 0.025, na.rm = TRUE),
  wolfN_high = apply(post_mat[, grep("^wolf_N_tot\\[", colnames(post_mat))], 2, quantile, probs = 0.975, na.rm = TRUE)
)

# yearly wolf survival summaries
wolf_pup_pts <- wolf_vrates2 %>%
  filter(rate == "Pup survival (s_p)") %>%
  transmute(
    year,
    wolf_surv = mean,
    wolf_low = low,
    wolf_high = high,
    stage = "Pup survival"
  )

wolf_ad_pts <- wolf_vrates2 %>%
  filter(rate == "Adult survival (s_a)") %>%
  transmute(
    year,
    wolf_surv = mean,
    wolf_low = low,
    wolf_high = high,
    stage = "Adult survival"
  )

wolf_surv_pts <- bind_rows(wolf_pup_pts, wolf_ad_pts)

# plotting dataframe
wolf_plot_df <- wolf_surv_pts %>%
  left_join(elk_abund_pts, by = "year") %>%
  left_join(wolf_abund_pts, by = "year") %>%
  left_join(
    covars %>%
      select(year, NR_Bison),
    by = "year"
  )

################################################################################
########------------- Elk abundance effect on wolf survival -------------#######
################################################################################

# x grid on raw elk female abundance scale
x_grid_elk <- seq(
  min(wolf_plot_df$elk_N_female, na.rm = TRUE),
  max(wolf_plot_df$elk_N_female, na.rm = TRUE),
  length.out = 200
)

# standardize internally because model used standardized elk abundance
x_grid_elk_std <- (x_grid_elk - icm_constants$elk_N_female_mean) / icm_constants$elk_N_female_sd

# hold bison and wolf abundance constant at mean standardized values
bisonN_mean_std <- mean(covars_std$NR_Bison, na.rm = TRUE)
wolf_mean_std <- mean(
  (wolf_plot_df$wolf_N_tot - icm_constants$wolf_tot_mean) / icm_constants$wolf_tot_sd,
  na.rm = TRUE
)

pred_wpup_elk <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_wpupSurv"] +
      post_mat[, "beta1_wpupSurv_elkN"] * x +
      post_mat[, "beta2_wpupSurv_bisonN"] * bisonN_mean_std +
      post_mat[, "beta3_wpupSurv_wolfN"] * wolf_mean_std
  )
})

pred_wad_elk <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_wadSurv"] +
      post_mat[, "beta1_wadSurv_elkN"] * x +
      post_mat[, "beta2_wadSurv_bisonN"] * bisonN_mean_std +
      post_mat[, "beta3_wadSurv_wolfN"] * wolf_mean_std
  )
})

line_df_elk_wolf <- bind_rows(
  data.frame(
    x = x_grid_elk,
    wolf_surv = apply(pred_wpup_elk, 2, mean),
    wolf_low = apply(pred_wpup_elk, 2, quantile, probs = 0.025, na.rm = TRUE),
    wolf_high = apply(pred_wpup_elk, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Pup survival"
  ),
  data.frame(
    x = x_grid_elk,
    wolf_surv = apply(pred_wad_elk, 2, mean),
    wolf_low = apply(pred_wad_elk, 2, quantile, probs = 0.025, na.rm = TRUE),
    wolf_high = apply(pred_wad_elk, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Adult survival"
  )
)

wolf_main_plot <- ggplot(wolf_plot_df, aes(x = elk_N_female, y = wolf_surv)) +
  geom_ribbon(
    data = line_df_elk_wolf,
    aes(x = x, ymin = wolf_low, ymax = wolf_high),
    inherit.aes = FALSE,
    fill = "#236192",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df_elk_wolf,
    aes(x = x, y = wolf_surv),
    inherit.aes = FALSE,
    color = "#236192",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = wolf_low, ymax = wolf_high), width = 0) +
  geom_errorbarh(aes(xmin = elk_low, xmax = elk_high), height = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Elk female abundance",
    y = "Wolf survival",
    title = "Estimated effect of elk abundance on wolf survival",
    subtitle = "Bison abundance and wolf abundance held at their means"
  )

wolf_main_plot

################################################################################
##########--------------- Bison effect on wolf survival ----------------########
################################################################################

# x grid on raw bison abundance scale
x_grid_bison_raw <- seq(
  min(covars$NR_Bison, na.rm = TRUE),
  max(covars$NR_Bison, na.rm = TRUE),
  length.out = 200
)

# convert raw bison grid to standardized values used in the model
bison_mean_raw <- mean(covars$NR_Bison, na.rm = TRUE)
bison_sd_raw <- sd(covars$NR_Bison, na.rm = TRUE)
x_grid_bison_std <- (x_grid_bison_raw - bison_mean_raw) / bison_sd_raw

# hold elk abundance and wolf abundance constant at mean standardized values
elk_mean_std <- mean(
  (wolf_plot_df$elk_N_female - icm_constants$elk_N_female_mean) / icm_constants$elk_N_female_sd,
  na.rm = TRUE
)

pred_wpup_bison <- sapply(x_grid_bison_std, function(x) {
  plogis(
    post_mat[, "beta0_wpupSurv"] +
      post_mat[, "beta1_wpupSurv_elkN"] * elk_mean_std +
      post_mat[, "beta2_wpupSurv_bisonN"] * x +
      post_mat[, "beta3_wpupSurv_wolfN"] * wolf_mean_std
  )
})

pred_wad_bison <- sapply(x_grid_bison_std, function(x) {
  plogis(
    post_mat[, "beta0_wadSurv"] +
      post_mat[, "beta1_wadSurv_elkN"] * elk_mean_std +
      post_mat[, "beta2_wadSurv_bisonN"] * x +
      post_mat[, "beta3_wadSurv_wolfN"] * wolf_mean_std
  )
})

line_df_bison_wolf <- bind_rows(
  data.frame(
    x = x_grid_bison_raw,
    wolf_surv = apply(pred_wpup_bison, 2, mean),
    wolf_low = apply(pred_wpup_bison, 2, quantile, probs = 0.025, na.rm = TRUE),
    wolf_high = apply(pred_wpup_bison, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Pup survival"
  ),
  data.frame(
    x = x_grid_bison_raw,
    wolf_surv = apply(pred_wad_bison, 2, mean),
    wolf_low = apply(pred_wad_bison, 2, quantile, probs = 0.025, na.rm = TRUE),
    wolf_high = apply(pred_wad_bison, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Adult survival"
  )
)

bison_plot_wolf <- ggplot(wolf_plot_df, aes(x = NR_Bison, y = wolf_surv)) +
  geom_ribbon(
    data = line_df_bison_wolf,
    aes(x = x, ymin = wolf_low, ymax = wolf_high),
    inherit.aes = FALSE,
    fill = "#6F263D",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df_bison_wolf,
    aes(x = x, y = wolf_surv),
    inherit.aes = FALSE,
    color = "#6F263D",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = wolf_low, ymax = wolf_high), width = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Bison abundance",
    y = "Wolf survival",
    title = "Estimated effect of bison abundance on wolf survival",
    subtitle = "Elk abundance and wolf abundance held at their means"
  )

bison_plot_wolf

################################################################################
##########------------ Wolf abundance effect on wolf survival ----------------##
################################################################################

# x grid on raw wolf abundance scale
x_grid_wolf_raw <- seq(
  min(wolf_plot_df$wolf_N_tot, na.rm = TRUE),
  max(wolf_plot_df$wolf_N_tot, na.rm = TRUE),
  length.out = 200
)

# standardize internally because model used standardized wolf abundance
x_grid_wolf_std <- (x_grid_wolf_raw - icm_constants$wolf_tot_mean) / icm_constants$wolf_tot_sd

# hold elk abundance and bison abundance constant at mean standardized values
pred_wpup_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_wpupSurv"] +
      post_mat[, "beta1_wpupSurv_elkN"] * elk_mean_std +
      post_mat[, "beta2_wpupSurv_bisonN"] * bisonN_mean_std +
      post_mat[, "beta3_wpupSurv_wolfN"] * x
  )
})

pred_wad_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_wadSurv"] +
      post_mat[, "beta1_wadSurv_elkN"] * elk_mean_std +
      post_mat[, "beta2_wadSurv_bisonN"] * bisonN_mean_std +
      post_mat[, "beta3_wadSurv_wolfN"] * x
  )
})

line_df_wolf_wolf <- bind_rows(
  data.frame(
    x = x_grid_wolf_raw,
    wolf_surv = apply(pred_wpup_wolf, 2, mean),
    wolf_low = apply(pred_wpup_wolf, 2, quantile, probs = 0.025, na.rm = TRUE),
    wolf_high = apply(pred_wpup_wolf, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Pup survival"
  ),
  data.frame(
    x = x_grid_wolf_raw,
    wolf_surv = apply(pred_wad_wolf, 2, mean),
    wolf_low = apply(pred_wad_wolf, 2, quantile, probs = 0.025, na.rm = TRUE),
    wolf_high = apply(pred_wad_wolf, 2, quantile, probs = 0.975, na.rm = TRUE),
    stage = "Adult survival"
  )
)

wolfN_plot_wolf <- ggplot(wolf_plot_df, aes(x = wolf_N_tot, y = wolf_surv)) +
  geom_ribbon(
    data = line_df_wolf_wolf,
    aes(x = x, ymin = wolf_low, ymax = wolf_high),
    inherit.aes = FALSE,
    fill = "#4B7F52",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df_wolf_wolf,
    aes(x = x, y = wolf_surv),
    inherit.aes = FALSE,
    color = "#4B7F52",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = wolf_low, ymax = wolf_high), width = 0) +
  geom_errorbarh(aes(xmin = wolfN_low, xmax = wolfN_high), height = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Wolf abundance",
    y = "Wolf survival",
    title = "Estimated effect of wolf abundance on wolf survival",
    subtitle = "Elk abundance and bison abundance held at their means"
  )

wolfN_plot_wolf

################################################################################
########--------------- Posterior coefficient densities ----------------########
################################################################################
################################################################################
#
##
###
####
#####
######
# specify coefficients of interest for elk density plots
coefs <- c("Calf intercept", "Calf grizzly effect",
           "OA intercept", "OA grizzly effect",
           "YA intercept", "YA grizzly effect")
######
#####
####
###
##
#

################################################################################

# elk coefficient densities
beta_df <- data.frame(
  beta0_calfSurv = post_mat[, "beta0_calfSurv"],
  beta1_calfSurv_wolfN = post_mat[, "beta1_calfSurv_wolfN"],
  beta2_calfSurv_wintPPT = post_mat[, "beta2_calfSurv_wintPPT"],
  beta3_calfSurv_grizN = post_mat[, "beta3_calfSurv_grizN"],
  beta0_yaSurv = post_mat[, "beta0_yaSurv"],
  beta1_yaSurv_wolfN = post_mat[, "beta1_yaSurv_wolfN"],
  beta2_yaSurv_wintPPT = post_mat[, "beta2_yaSurv_wintPPT"],
  beta3_yaSurv_grizN = post_mat[, "beta3_yaSurv_grizN"],
  beta0_oaSurv = post_mat[, "beta0_oaSurv"],
  beta1_oaSurv_wolfN = post_mat[, "beta1_oaSurv_wolfN"],
  beta2_oaSurv_wintPPT = post_mat[, "beta2_oaSurv_wintPPT"],
  beta3_oaSurv_grizN = post_mat[, "beta3_oaSurv_grizN"]
)

beta_long <- beta_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

beta_long$parameter <- factor(
  beta_long$parameter,
  levels = c(
    "beta0_calfSurv", "beta1_calfSurv_wolfN", "beta2_calfSurv_wintPPT", "beta3_calfSurv_grizN",
    "beta0_yaSurv", "beta1_yaSurv_wolfN", "beta2_yaSurv_wintPPT", "beta3_yaSurv_grizN",
    "beta0_oaSurv", "beta1_oaSurv_wolfN", "beta2_oaSurv_wintPPT", "beta3_oaSurv_grizN"
  ),
  labels = c(
    "Calf intercept", "Calf wolf effect", "Calf winterPPT effect", "Calf grizzly effect",
    "YA intercept", "YA wolf effect", "YA winterPPT effect", "YA grizzly effect",
    "OA intercept", "OA wolf effect", "OA winterPPT effect", "OA grizzly effect"
  )
)

coef_plot <- ggplot(beta_long[beta_long$parameter %in% coefs,], aes(x = value)) +
  geom_density(fill = "#236192", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~parameter, scales = "free", ncol = 6) +
  theme_classic() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distributions of elk regression coefficients"
  )

coef_plot

################################################################################

#
##
###
####
#####
######
# specify coefficients of interest for elk density plots
coefs_wolf <- c("Wolf pup intercept", "Wolf pup bison effect",
                "Wolf adult intercept", "Wolf adult bison effect")
######
#####
####
###
##
#

################################################################################

# wolf coefficient densities
wolf_beta_df <- data.frame(
  beta0_wpupSurv = post_mat[, "beta0_wpupSurv"],
  beta1_wpupSurv_elkN = post_mat[, "beta1_wpupSurv_elkN"],
  beta2_wpupSurv_bisonN = post_mat[, "beta2_wpupSurv_bisonN"],
  beta0_wadSurv = post_mat[, "beta0_wadSurv"],
  beta1_wadSurv_elkN = post_mat[, "beta1_wadSurv_elkN"],
  beta2_wadSurv_bisonN = post_mat[, "beta2_wadSurv_bisonN"]
)

wolf_beta_long <- wolf_beta_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

wolf_beta_long$parameter <- factor(
  wolf_beta_long$parameter,
  levels = c(
    "beta0_wpupSurv", "beta1_wpupSurv_elkN", "beta2_wpupSurv_bisonN",
    "beta0_wadSurv", "beta1_wadSurv_elkN", "beta2_wadSurv_bisonN"
  ),
  labels = c(
    "Wolf pup intercept", "Wolf pup elk effect", "Wolf pup bison effect",
    "Wolf adult intercept", "Wolf adult elk effect", "Wolf adult bison effect"
  )
)

wolf_coef_plot <- ggplot(wolf_beta_long[wolf_beta_long$parameter %in% coefs_wolf,], aes(x = value)) +
  geom_density(fill = "#6F263D", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~parameter, scales = "free", ncol = 4) +
  theme_classic() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distributions of wolf regression coefficients"
  )

wolf_coef_plot

################################################################################
#########------------------------ Combine plots -------------------------#######
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

wolf_final_plot <- plot_grid(
  wolf_main_plot,
  bison_plot_wolf,
  wolfN_plot_wolf,
  wolf_coef_plot,
  ncol = 1,
  rel_heights = c(1.2, 1.2, 1.2, 1),
  labels = c("E", "F", "G", "H")
)

all_regression_plots <- plot_grid(
  elk_panel,
  wolf_final_plot,
  ncol = 2,
  rel_widths = c(1, 1)
)

all_regression_plots

################################################################################
##########------ Summarize coefficient evidence using nested CIs ------##########
################################################################################

# This chunk classifies each coefficient as:
# - strong evidence: 95% CI excludes 0
# - moderate evidence: 95% CI includes 0, but 80% CI excludes 0
# - weak evidence: 80% CI includes 0, but 50% CI excludes 0
# - little/no evidence: 50% CI includes 0

# assumes you already have:
# post_mat = posterior draws matrix
# OR another matrix/data.frame of posterior samples with parameters in columns

# ------------------------------------------------------------------------------
# 1) choose coefficients to summarize
# ------------------------------------------------------------------------------

coef_names <- c(
  # elk coefficients
  "beta0_calfSurv", "beta1_calfSurv_wolfN", "beta2_calfSurv_wintPPT", "beta3_calfSurv_grizN", "beta4_calfSurv_elkN",
  "beta0_yaSurv", "beta1_yaSurv_wolfN", "beta2_yaSurv_wintPPT", "beta3_yaSurv_grizN", "beta4_yaSurv_elkN",
  "beta0_oaSurv", "beta1_oaSurv_wolfN", "beta2_oaSurv_wintPPT", "beta3_oaSurv_grizN", "beta4_oaSurv_elkN",
  
  # wolf coefficients
  "beta0_wpupSurv", "beta1_wpupSurv_elkN", "beta2_wpupSurv_bisonN", "beta3_wpupSurv_wolfN",
  "beta0_wadSurv", "beta1_wadSurv_elkN", "beta2_wadSurv_bisonN", "beta3_wadSurv_wolfN"
)

coef_names <- coef_names[coef_names %in% colnames(post_mat)]

# ------------------------------------------------------------------------------
# 2) optional pretty labels
# ------------------------------------------------------------------------------

pretty_labels <- c(
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

# ------------------------------------------------------------------------------
# 3) helper function
# ------------------------------------------------------------------------------

classify_evidence <- function(draws) {
  ci95 <- quantile(draws, probs = c(0.025, 0.975), na.rm = TRUE)
  ci80 <- quantile(draws, probs = c(0.10, 0.90), na.rm = TRUE)
  ci50 <- quantile(draws, probs = c(0.25, 0.75), na.rm = TRUE)
  
  excludes_zero_95 <- ci95[1] > 0 | ci95[2] < 0
  excludes_zero_80 <- ci80[1] > 0 | ci80[2] < 0
  excludes_zero_50 <- ci50[1] > 0 | ci50[2] < 0
  
  evidence <- dplyr::case_when(
    excludes_zero_95 ~ "Strong",
    !excludes_zero_95 & excludes_zero_80 ~ "Moderate",
    !excludes_zero_80 & excludes_zero_50 ~ "Weak",
    TRUE ~ "Little/none"
  )
  
  direction <- dplyr::case_when(
    median(draws, na.rm = TRUE) > 0 ~ "Positive",
    median(draws, na.rm = TRUE) < 0 ~ "Negative",
    TRUE ~ "Neutral"
  )
  
  data.frame(
    mean = mean(draws, na.rm = TRUE),
    median = median(draws, na.rm = TRUE),
    ci95_low = ci95[1],
    ci95_high = ci95[2],
    ci80_low = ci80[1],
    ci80_high = ci80[2],
    ci50_low = ci50[1],
    ci50_high = ci50[2],
    excludes_zero_95 = excludes_zero_95,
    excludes_zero_80 = excludes_zero_80,
    excludes_zero_50 = excludes_zero_50,
    direction = direction,
    evidence = evidence
  )
}

# ------------------------------------------------------------------------------
# 4) apply to all coefficients
# ------------------------------------------------------------------------------

coef_evidence <- lapply(coef_names, function(param) {
  out <- classify_evidence(post_mat[, param])
  out$parameter <- param
  out$label <- ifelse(param %in% names(pretty_labels), pretty_labels[param], param)
  out
})

coef_evidence <- bind_rows(coef_evidence) %>%
  select(
    parameter,
    label,
    mean,
    median,
    ci95_low,
    ci95_high,
    ci80_low,
    ci80_high,
    ci50_low,
    ci50_high,
    direction,
    evidence,
    excludes_zero_95,
    excludes_zero_80,
    excludes_zero_50
  ) %>%
  arrange(
    factor(evidence, levels = c("Strong", "Moderate", "Weak", "Little/none")),
    desc(abs(median))
  )

# ------------------------------------------------------------------------------
# 5) print concise summary table
# ------------------------------------------------------------------------------

coef_evidence %>%
  mutate(
    across(c(mean, median, ci95_low, ci95_high, ci80_low, ci80_high, ci50_low, ci50_high), round, 3)
  )

# ------------------------------------------------------------------------------
# 6) split into groups if useful
# ------------------------------------------------------------------------------

strong_effects <- coef_evidence %>% filter(evidence == "Strong")
moderate_effects <- coef_evidence %>% filter(evidence == "Moderate")
weak_effects <- coef_evidence %>% filter(evidence == "Weak")
little_effects <- coef_evidence %>% filter(evidence == "Little/none")

strong_effects
moderate_effects
weak_effects
little_effects

# ------------------------------------------------------------------------------
# 7) optional compact text summary
# ------------------------------------------------------------------------------

cat("\nSTRONG EVIDENCE\n")
print(strong_effects %>% select(label, median, ci95_low, ci95_high, direction))

cat("\nMODERATE EVIDENCE\n")
print(moderate_effects %>% select(label, median, ci80_low, ci80_high, direction))

cat("\nWEAK EVIDENCE\n")
print(weak_effects %>% select(label, median, ci50_low, ci50_high, direction))

cat("\nLITTLE / NO EVIDENCE\n")
print(little_effects %>% select(label, median, ci50_low, ci50_high, direction))

################################################################################
##########------------------ Elasticity analysis ---------------------##########
################################################################################

# ------------------------------------------------------------------------------
# Build annual elk matrices, lambda, and matrix-transition elasticities
# ------------------------------------------------------------------------------

elk_rates_wide <- elk_vrates2 %>%
  filter(rate %in% c(
    "Calf survival (s_c)",
    "Young Adult survival (s_ya)",
    "Old Adult survival (s_oa)",
    "Young→Old transition (p_13)",
    "Fecundity (young) (f_ya)",
    "Fecundity (old) (f_oa)"
  )) %>%
  select(year, rate, mean) %>%
  pivot_wider(names_from = rate, values_from = mean) %>%
  drop_na()

elk_elasticity_list <- lapply(1:nrow(elk_rates_wide), function(i) {
  
  s_c  <- elk_rates_wide$`Calf survival (s_c)`[i]
  s_ya <- elk_rates_wide$`Young Adult survival (s_ya)`[i]
  s_oa <- elk_rates_wide$`Old Adult survival (s_oa)`[i]
  p_13 <- elk_rates_wide$`Young→Old transition (p_13)`[i]
  f_ya <- elk_rates_wide$`Fecundity (young) (f_ya)`[i]
  f_oa <- elk_rates_wide$`Fecundity (old) (f_oa)`[i]
  
  A <- matrix(
    c(
      0,      f_ya * s_c,          f_oa * s_c,
      s_ya,   s_ya * (1 - p_13),   0,
      0,      s_ya * p_13,         s_oa
    ),
    nrow = 3,
    byrow = TRUE
  )
  
  E <- popbio::elasticity(A)
  lambda <- popbio::lambda(A)
  
  data.frame(
    year = elk_rates_wide$year[i],
    lambda = lambda,
    a12 = E[1, 2],
    a13 = E[1, 3],
    a21 = E[2, 1],
    a22 = E[2, 2],
    a32 = E[3, 2],
    a33 = E[3, 3]
  )
})

elk_elasticity_df <- bind_rows(elk_elasticity_list)

# ------------------------------------------------------------------------------
# Clean long-format elasticity dataframe with readable transition labels
# ------------------------------------------------------------------------------

elk_elasticity_long <- elk_elasticity_df %>%
  pivot_longer(
    cols = c(a12, a13, a21, a22, a32, a33),
    names_to = "transition",
    values_to = "elasticity"
  ) %>%
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

# ------------------------------------------------------------------------------
# 1) Lambda through time
# ------------------------------------------------------------------------------

lambda_plot <- ggplot(elk_elasticity_df, aes(x = year, y = lambda)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    x = "Year",
    y = expression(lambda),
    title = "Annual elk population growth rate from projection matrix",
    subtitle = "Dashed line at lambda = 1 indicates stable population size"
  )

lambda_plot

# ------------------------------------------------------------------------------
# 2) Elasticity of lambda to matrix transitions through time
# ------------------------------------------------------------------------------

elasticity_time_plot <- ggplot(elk_elasticity_long, aes(x = year, y = elasticity)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ transition, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Year",
    y = "Elasticity",
    title = "Elasticity of elk population growth to matrix transitions"
  ) +
  theme(legend.position = "none")

elasticity_time_plot

# ------------------------------------------------------------------------------
# 3) Mean elasticity across years
# ------------------------------------------------------------------------------

elk_elasticity_summary <- elk_elasticity_long %>%
  group_by(transition) %>%
  summarise(
    mean_elasticity = mean(elasticity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_elasticity))

elasticity_bar_plot <- ggplot(elk_elasticity_summary,
                              aes(x = reorder(transition, mean_elasticity),
                                  y = mean_elasticity)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean elasticity",
    title = "Mean elasticity of elk population growth to matrix transitions"
  )

elasticity_bar_plot

# ------------------------------------------------------------------------------
# 4) Lambda through time and mean elasticity across years
# ------------------------------------------------------------------------------

elasticity_combo <- plot_grid(
  lambda_plot,
  elasticity_bar_plot,
  ncol = 1,
  rel_heights = c(1, 1.1),
  labels = c("A", "B")
)

elasticity_combo

################################################################################
##########------------- Wolf population change summaries -------------###########
################################################################################

# pull posterior abundance summaries by stage
wolf_total_df <- wolf_N_summ %>%
  filter(stage == "Total Wolves") %>%
  arrange(year)

wolf_pup_df <- wolf_N_summ %>%
  filter(stage == "Pups") %>%
  arrange(year)

wolf_adult_df <- wolf_N_summ %>%
  filter(stage == "Adults") %>%
  arrange(year)

# first and last year values
wolf_first_year <- min(wolf_total_df$year)
wolf_last_year <- max(wolf_total_df$year)

wolf_first_total <- wolf_total_df$mean[wolf_total_df$year == wolf_first_year]
wolf_last_total <- wolf_total_df$mean[wolf_total_df$year == wolf_last_year]

# percent change in total abundance
wolf_pct_change <- 100 * (wolf_last_total - wolf_first_total) / wolf_first_total

# geometric annual growth rate lambda
wolf_n_intervals <- wolf_last_year - wolf_first_year
wolf_lambda_geom <- (wolf_last_total / wolf_first_total)^(1 / wolf_n_intervals)

# annual percent change
wolf_annual_pct_change <- (wolf_lambda_geom - 1) * 100

################################################################################
##########------------- Wolf stage structure summaries ---------------###########
################################################################################

# merge pup and adult summaries
wolf_stage_df <- wolf_pup_df %>%
  select(year, pup_mean = mean) %>%
  left_join(
    wolf_adult_df %>%
      select(year, adult_mean = mean),
    by = "year"
  ) %>%
  mutate(
    total_mean = pup_mean + adult_mean,
    prop_pup = pup_mean / total_mean,
    prop_adult = adult_mean / total_mean
  )

# first and last year stage structure
wolf_first_prop_pup <- wolf_stage_df$prop_pup[wolf_stage_df$year == wolf_first_year]
wolf_last_prop_pup <- wolf_stage_df$prop_pup[wolf_stage_df$year == wolf_last_year]

wolf_first_prop_adult <- wolf_stage_df$prop_adult[wolf_stage_df$year == wolf_first_year]
wolf_last_prop_adult <- wolf_stage_df$prop_adult[wolf_stage_df$year == wolf_last_year]

################################################################################
##########---------------------- Print results ------------------------###########
################################################################################

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