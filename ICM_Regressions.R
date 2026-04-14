### Integrated Community Model (ICM) - Stage 2: Regressions
### Uses estimated abundances, demographic rates from stage 1 in regressions
### Last updated: Apr. 13, 2026
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
############--------------- Load covariates ---------------#####################
################################################################################

# load estimated demographic data from ICM stage 1
dat <- read.csv('data/outputs/wolf_elk_ICM_estimates_2026Apr13.csv')

# then covariates
annual_prism <- read.csv('data/covariates/prism_annual_precip_tmean.csv')
bison <- read.csv('data/covariates/NR_Bison_Abundance.csv') %>%
  rename(year = Year)

# bind data into one dataframe 'covars"
covars <- left_join(annual_prism, bison)

# trim covariate data to shared years 
covars <- covars[covars$year %in% community_years,]

# standardize all covariates except year
covars_std <- covars %>%
  mutate(
    across(-year, ~ as.numeric(scale(.)))
  )

# then bind to the wolf and elk abundances/demographic rates
stage2_dat <- dat %>%
  left_join(covars_std, by = "year")

stage2_dat

################################################################################
###########---------------- Regression models -----------------#################
################################################################################

icm_code <- nimbleCode({
  
  ##########---------------- REGRESSIONS ----------------#############
  
  # Priors for elk regression coefficients
  beta0_calfSurv ~ dnorm(qlogis(0.22), 1 / 0.3^2) # mean = 0.22
  beta1_calfSurv ~ dnorm(0, 1 / 0.3^2)
  
  beta0_yaSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2) # mean = 0.90
  beta1_yaSurv ~ dnorm(0, 1 / 0.3^2)
  
  beta0_oaSurv ~ dnorm(qlogis(0.80), 1 / 0.3^2) # mean = 0.80
  beta1_oaSurv ~ dnorm(0, 1 / 0.3^2)
  
  # Priors for elk random year-effect SDs
  sigma_calf ~ dunif(0, 0.5)
  tau_calf <- 1 / (sigma_calf^2)
  
  sigma_ya ~ dunif(0, 0.3)
  tau_ya <- 1 / (sigma_ya^2)
  
  sigma_oa ~ dunif(0, 0.3)
  tau_oa <- 1 / (sigma_oa^2)
  
  # Priors for regression coefficients
  beta0_wpupSurv ~ dnorm(qlogis(0.5), 1 / 0.3^2) # mean = 0.22
  beta1_wpupSurv ~ dnorm(0, 1 / 0.3^2)
  
  beta0_wadSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2) # mean = 0.90
  beta1_wadSurv ~ dnorm(0, 1 / 0.3^2)
  
  # Priors for wolf random year-effect SDs
  sigma_wpup ~ dunif(0, 0.5)
  tau_wpup <- 1 / (sigma_wpup^2)
  
  sigma_wad ~ dunif(0, 0.3)
  tau_wad <- 1 / (sigma_wad^2)
  
  # Year-specific regressions
  for (t in 1:n_years) {
    
    # elk random year effects
    eps_elk_s_c[t] ~ dnorm(0, tau_calf)
    eps_elk_s_ya[t] ~ dnorm(0, tau_ya)
    eps_elk_s_oa[t] ~ dnorm(0, tau_oa)
    
    # elk random year effects
    eps_wolf_s_p[t] ~ dnorm(0, tau_wpup)
    eps_wolf_s_a[t] ~ dnorm(0, tau_wad)
    
    # standardize wolf abundance
    wolf_N_tot_std[t] <- (wolf_N_tot[t] - wolf_tot_mean) / wolf_tot_sd
    
    # standardize elk abundance
    elk_N_female_std[t] <- (elk_N_female[t] - elk_N_female_mean) / elk_N_female_sd
    
    # elk regression models
    logit(elk_s_c[t])  <- beta0_calfSurv + beta1_calfSurv * wolf_N_tot_std[t] + eps_elk_s_c[t]
    logit(elk_s_ya[t]) <- beta0_yaSurv   + beta1_yaSurv   * wolf_N_tot_std[t] + eps_elk_s_ya[t]
    logit(elk_s_oa[t]) <- beta0_oaSurv   + beta1_oaSurv   * wolf_N_tot_std[t] + eps_elk_s_oa[t]
    
    # wolf regression models
    logit(wolf_s_p[t])  <- beta0_wpupSurv + beta1_wpupSurv * elk_N_female_std[t] + eps_wolf_s_p[t]
    logit(wolf_s_a[t]) <- beta0_wadSurv   + beta1_wadSurv   * elk_N_female_std[t] + eps_wolf_s_a[t]
  }
})

icm_constants <- list(
  n_years = n_years,
  # these are added to allow the abundance metrics to be standardized
  wolf_tot_mean = mean(wolf_pop$total_abundance, na.rm = TRUE),
  wolf_tot_sd = sd(wolf_pop$total_abundance, na.rm = TRUE),
  elk_N_female_mean = mean(elk_dat_n$n_female, na.rm=TRUE),
  elk_N_female_sd =sd(elk_dat_n$n_female, na.rm=TRUE)
  )

params <- c(
  # elk regression coefficients
  'beta0_calfSurv', 'beta1_calfSurv',
  'beta0_yaSurv', 'beta1_yaSurv',
  'beta0_oaSurv', 'beta1_oaSurv',
  
  "sigma_calf", "sigma_ya", "sigma_oa",
  "eps_elk_s_c", "eps_elk_s_ya", "eps_elk_s_oa",
  
  # wolf regression coefficients
  'beta0_wpupSurv', 'beta1_wpupSurv',
  'beta0_wadSurv', 'beta1_wadSurv',
  
  "sigma_wpup", "sigma_wad",
  "eps_wolf_s_p", "eps_wolf_s_a"
)


################################################################################
###########--------------------- Results ----------------------#################
################################################################################

# simple summary
round(icm_mod$summary$all.chains, 2)

# pull nimble summary table
post_sum <- as.data.frame(icm_mod$summary$all.chains)
post_sum$param <- rownames(post_sum)

# =============================================================================
# Clean icm_mod$samples so MCMCvis works
# =============================================================================

# 1) extract each chain as a plain matrix
mats <- lapply(icm_mod$samples, function(ch) as.matrix(ch))

# 2) trim all chains to the same length
lens <- sapply(mats, nrow)
L <- min(lens)

mats_trim <- lapply(mats, function(M) {
  tail(M, L)
})

# 3) rebuild a consistent mcmc.list
icm_mlist <- mcmc.list(lapply(mats_trim, function(M) {
  mcmc(M, start = 1, end = L, thin = 1)
}))

# 4) remove any parameters with NA/NaN in any chain
mats2 <- lapply(icm_mlist, as.matrix)

keep_cols <- Reduce(intersect, lapply(mats2, function(M) {
  colnames(M)[colSums(is.na(M) | is.nan(M)) == 0]
}))

icm_clean <- mcmc.list(lapply(mats2, function(M) {
  mcmc(M[, keep_cols, drop = FALSE], start = 1, end = nrow(M), thin = 1)
}))

# quick check
# round(MCMCsummary(icm_clean, params = "all"), 2)

################################################################################
###########------------ Checking model convergence ------------#################
################################################################################

# regression coefficients
reg_coefs <- MCMCsummary(
  icm_clean,
  params = c(
    "beta0_calfSurv", "beta1_calfSurv",
    "beta0_yaSurv", "beta1_yaSurv",
    "beta0_oaSurv", "beta1_oaSurv",
    "beta0_wpupSurv", "beta1_wpupSurv",
    "beta0_wadSurv", "beta1_wadSurv"
  )
)

bad_reg_coefs <- data.frame(
  param = rownames(reg_coefs),
  Rhat = reg_coefs[, "Rhat"]
)

bad_reg_coefs <- bad_reg_coefs[bad_reg_coefs$Rhat > 1.1, ]
bad_reg_coefs

################################################################################
##########---------- Regression plots + coefficient densities ---------##########
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

# combine elk survival summaries
elk_surv_pts <- bind_rows(elk_calf_pts, elk_ya_pts, elk_oa_pts)

# join elk and wolf estimates
plot_df <- elk_surv_pts %>%
  left_join(wolf_pts, by = "year")

# pull posterior draws from cleaned chains
post_mat <- do.call(rbind, lapply(icm_clean, as.matrix))

# x grid on raw wolf abundance scale
x_grid <- seq(
  min(plot_df$wolf_N_tot, na.rm = TRUE),
  max(plot_df$wolf_N_tot, na.rm = TRUE),
  length.out = 200
)

# standardize internally because model used standardized wolf abundance
x_grid_std <- (x_grid - icm_constants$wolf_tot_mean) / icm_constants$wolf_tot_sd

# fitted values for calf survival
pred_calf <- sapply(x_grid_std, function(x) {
  plogis(post_mat[, "beta0_calfSurv"] + post_mat[, "beta1_calfSurv"] * x)
})

line_calf <- data.frame(
  wolf_N_tot = x_grid,
  elk_surv = apply(pred_calf, 2, mean),
  elk_low = apply(pred_calf, 2, quantile, probs = 0.025),
  elk_high = apply(pred_calf, 2, quantile, probs = 0.975),
  stage = "Calf survival"
)

# fitted values for young adult survival
pred_ya <- sapply(x_grid_std, function(x) {
  plogis(post_mat[, "beta0_yaSurv"] + post_mat[, "beta1_yaSurv"] * x)
})

line_ya <- data.frame(
  wolf_N_tot = x_grid,
  elk_surv = apply(pred_ya, 2, mean),
  elk_low = apply(pred_ya, 2, quantile, probs = 0.025),
  elk_high = apply(pred_ya, 2, quantile, probs = 0.975),
  stage = "Young adult survival"
)

# fitted values for old adult survival
pred_oa <- sapply(x_grid_std, function(x) {
  plogis(post_mat[, "beta0_oaSurv"] + post_mat[, "beta1_oaSurv"] * x)
})

line_oa <- data.frame(
  wolf_N_tot = x_grid,
  elk_surv = apply(pred_oa, 2, mean),
  elk_low = apply(pred_oa, 2, quantile, probs = 0.025),
  elk_high = apply(pred_oa, 2, quantile, probs = 0.975),
  stage = "Old adult survival"
)

# combine survival values from all stages
line_df <- bind_rows(line_calf, line_ya, line_oa)

# first plot: wolf abundance vs. all three stages' survival 
main_plot <- ggplot(plot_df, aes(x = wolf_N_tot, y = elk_surv)) +
  geom_ribbon(
    data = line_df,
    aes(x = wolf_N_tot, ymin = elk_low, ymax = elk_high),
    inherit.aes = FALSE,
    fill = "#6F263D",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df,
    aes(x = wolf_N_tot, y = elk_surv),
    inherit.aes = FALSE,
    color = "#6F263D",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = elk_low, ymax = elk_high), width = 0) +
  geom_errorbarh(aes(xmin = wolf_low, xmax = wolf_high), height = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Wolf abundance",
    y = "Elk survival",
    title = "Estimated effect of wolf abundance on elk survival"
  )

# view main plot
main_plot

# then extract posterior densities of coefficients
beta_df <- data.frame(
  beta0_calfSurv = post_mat[, "beta0_calfSurv"],
  beta1_calfSurv = post_mat[, "beta1_calfSurv"],
  beta0_yaSurv = post_mat[, "beta0_yaSurv"],
  beta1_yaSurv = post_mat[, "beta1_yaSurv"],
  beta0_oaSurv = post_mat[, "beta0_oaSurv"],
  beta1_oaSurv = post_mat[, "beta1_oaSurv"]
)

# reshape dataframe
beta_long <- beta_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

# relabel
beta_long$parameter <- factor(
  beta_long$parameter,
  levels = c(
    "beta0_calfSurv", "beta1_calfSurv",
    "beta0_yaSurv", "beta1_yaSurv",
    "beta0_oaSurv", "beta1_oaSurv"
  ),
  labels = c(
    "Calf intercept", "Calf wolf effect",
    "YA intercept", "YA wolf effect",
    "OA intercept", "OA wolf effect"
  )
)

# generate coefficient estimate density plots
coef_plot <- ggplot(beta_long, aes(x = value)) +
  geom_density(fill = "#236192", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~parameter, scales = "free", ncol = 2) +
  theme_classic() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distributions of regression coefficients"
  )

# view coef plot
coef_plot

# piece the main_plot and coef_plot together 
final_plot <- plot_grid(
  main_plot,
  coef_plot,
  ncol = 1,
  rel_heights = c(2, 1.5),
  labels = c("A", "C")
)

# view final plot
final_plot

################################################################################
##########---------- Wolf regression plots + coefficient densities ----##########
################################################################################

# yearly elk abundance summaries for wolf regressions
elk_abund_pts <- elk_N_summ %>%
  filter(stage == "Total Females") %>%
  transmute(
    year,
    elk_N_female = mean,
    elk_low = low,
    elk_high = high
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

# combine wolf survival summaries
wolf_surv_pts <- bind_rows(wolf_pup_pts, wolf_ad_pts)

# join wolf survival and elk abundance
wolf_plot_df <- wolf_surv_pts %>%
  left_join(elk_abund_pts, by = "year")

# x grid on raw elk female abundance scale
x_grid_wolf <- seq(
  min(wolf_plot_df$elk_N_female, na.rm = TRUE),
  max(wolf_plot_df$elk_N_female, na.rm = TRUE),
  length.out = 200
)

# standardize internally because model used standardized elk abundance
x_grid_elk_std <- (x_grid_wolf - icm_constants$elk_N_female_mean) / icm_constants$elk_N_female_sd

# fitted values for wolf pup survival
pred_wpup <- sapply(x_grid_elk_std, function(x) {
  plogis(post_mat[, "beta0_wpupSurv"] + post_mat[, "beta1_wpupSurv"] * x)
})

line_wpup <- data.frame(
  elk_N_female = x_grid_wolf,
  wolf_surv = apply(pred_wpup, 2, mean),
  wolf_low = apply(pred_wpup, 2, quantile, probs = 0.025),
  wolf_high = apply(pred_wpup, 2, quantile, probs = 0.975),
  stage = "Pup survival"
)

# fitted values for wolf adult survival
pred_wad <- sapply(x_grid_elk_std, function(x) {
  plogis(post_mat[, "beta0_wadSurv"] + post_mat[, "beta1_wadSurv"] * x)
})

line_wad <- data.frame(
  elk_N_female = x_grid_wolf,
  wolf_surv = apply(pred_wad, 2, mean),
  wolf_low = apply(pred_wad, 2, quantile, probs = 0.025),
  wolf_high = apply(pred_wad, 2, quantile, probs = 0.975),
  stage = "Adult survival"
)

# combine wolf fitted values
wolf_line_df <- bind_rows(line_wpup, line_wad)

# main wolf regression plot
wolf_main_plot <- ggplot(wolf_plot_df, aes(x = elk_N_female, y = wolf_surv)) +
  geom_ribbon(
    data = wolf_line_df,
    aes(x = elk_N_female, ymin = wolf_low, ymax = wolf_high),
    inherit.aes = FALSE,
    fill = "#236192",
    alpha = 0.25
  ) +
  geom_line(
    data = wolf_line_df,
    aes(x = elk_N_female, y = wolf_surv),
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
    title = "Estimated effect of elk abundance on wolf survival"
  )

wolf_main_plot

# posterior densities of wolf regression coefficients
wolf_beta_df <- data.frame(
  beta0_wpupSurv = post_mat[, "beta0_wpupSurv"],
  beta1_wpupSurv = post_mat[, "beta1_wpupSurv"],
  beta0_wadSurv = post_mat[, "beta0_wadSurv"],
  beta1_wadSurv = post_mat[, "beta1_wadSurv"]
)

wolf_beta_long <- wolf_beta_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

wolf_beta_long$parameter <- factor(
  wolf_beta_long$parameter,
  levels = c(
    "beta0_wpupSurv", "beta1_wpupSurv",
    "beta0_wadSurv", "beta1_wadSurv"
  ),
  labels = c(
    "Wolf pup intercept", "Wolf pup elk effect",
    "Wolf adult intercept", "Wolf adult elk effect"
  )
)

wolf_coef_plot <- ggplot(wolf_beta_long, aes(x = value)) +
  geom_density(fill = "#6F263D", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~parameter, scales = "free", ncol = 2) +
  theme_classic() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distributions of wolf regression coefficients"
  )

wolf_coef_plot

# combine main wolf regression plot + coefficient densities
wolf_final_plot <- plot_grid(
  wolf_main_plot,
  wolf_coef_plot,
  ncol = 1,
  rel_heights = c(2, 1.4),
  labels = c("B", "D")
)

wolf_final_plot


################################################################################
##########------- Combining wolf and elk regression results ----------##########
################################################################################

all_regression_plots <- plot_grid(
  final_plot,
  wolf_final_plot,
  ncol = 2,
  rel_heights = c(1, 1)
)

all_regression_plots

