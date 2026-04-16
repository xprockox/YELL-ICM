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
##########--------------- Load data and clean ---------------###################
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
covars <- covars[covars$year %in% dat$year,]

# standardize all covariates except year
covars_std <- covars %>%
  mutate(
    across(-year, ~ as.numeric(scale(.)))
  )

# then bind to the wolf and elk abundances/demographic rates to the covariates
dat2 <- dat %>%
  select(-X) %>%
  left_join(covars_std, by = "year")

# function to force all rate values to be bounded by 0.000001 and 0.999999 
# (logis functions e.g. qlogis() break when using 0 or 1)
clamp01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

# filter all probability data through clamp01 function and rename for cleaner code
reg_df <- dat2 %>%
  transmute(
    year,
    
    elk_s_c = clamp01(mean_Calf.survival..s_c.),
    elk_s_ya = clamp01(mean_Young.Adult.survival..s_ya.),
    elk_s_oa = clamp01(mean_Old.Adult.survival..s_oa.),
    
    wolf_s_p = clamp01(mean_Pup.survival..s_p.),
    wolf_s_a = clamp01(mean_Adult.survival..s_a.),
    
    wolf_N_tot = wolf_N_tot,
    elk_N_female = elk_N_female,
    
    annual_ppt_mm_raw = covars$annual_ppt_mm[match(year, covars$year)],
    summer_ppt_mm_raw = covars$summer_ppt_mm[match(year, covars$year)],
    winter_ppt_mm_raw = covars$winter_ppt_mm[match(year, covars$year)],
    summer_tmean_c_raw = covars$summer_tmean_c[match(year, covars$year)],
    winter_tmean_c_raw = covars$winter_tmean_c[match(year, covars$year)],
    NR_Bison_raw = covars$NR_Bison[match(year, covars$year)],
    
    annual_ppt_mm_std = annual_ppt_mm,
    summer_ppt_mm_std = summer_ppt_mm,
    winter_ppt_mm_std = winter_ppt_mm,
    summer_tmean_c_std = summer_tmean_c,
    winter_tmean_c_std = winter_tmean_c,
    NR_Bison_std = NR_Bison
  ) %>%
  mutate(
    wolf_N_tot_std = as.numeric(scale(wolf_N_tot)),
    elk_N_female_std = as.numeric(scale(elk_N_female))
  )

################################################################################
###########---------------- Regression models -----------------#################
################################################################################

stage2_code <- nimbleCode({
  
  ################ ELK PRIORS ################
  
  # calf survival priors
  beta0_calf ~ dnorm(qlogis(0.22), 1 / 0.3^2)
  beta1_calf_wolf ~ dnorm(0, 1 / 0.3^2)
  beta2_calf_winterppt ~ dnorm(0, 1 / 0.3^2)
  
  sigma_calf ~ dunif(0, 0.5)
  tau_calf <- 1 / (sigma_calf^2)
  
  # young adult survival priors
  beta0_ya ~ dnorm(qlogis(0.90), 1 / 0.3^2)
  beta1_ya_wolf ~ dnorm(0, 1 / 0.3^2)
  beta2_ya_winterppt ~ dnorm(0, 1 / 0.3^2)
  
  sigma_ya ~ dunif(0, 0.3)
  tau_ya <- 1 / (sigma_ya^2)
  
  # old adult survival priors
  beta0_oa ~ dnorm(qlogis(0.80), 1 / 0.3^2)
  beta1_oa_wolf ~ dnorm(0, 1 / 0.3^2)
  beta2_oa_winterppt ~ dnorm(0, 1 / 0.3^2)
  
  sigma_oa ~ dunif(0, 0.3)
  tau_oa <- 1 / (sigma_oa^2)
  
  ################ WOLF PRIORS ################
  
  # pup survival priors
  beta0_wp ~ dnorm(qlogis(0.5), 1 / 0.3^2)
  beta1_wp_elk ~ dnorm(0, 1 / 0.3^2)
  beta2_wp_bison ~ dnorm(0, 1 / 0.3^2)
  
  sigma_wp ~ dunif(0, 0.5)
  tau_wp <- 1 / (sigma_wp^2)
  
  # adult survival priors
  beta0_wa ~ dnorm(qlogis(0.9), 1 / 0.3^2)
  beta1_wa_elk ~ dnorm(0, 1 / 0.3^2)
  beta2_wa_bison ~ dnorm(0, 1 / 0.3^2)
  
  sigma_wa ~ dunif(0, 0.3)
  tau_wa <- 1 / (sigma_wa^2)
  
  ################ REGRESSIONS ################
  
  for (t in 1:n_years) {
    
    # elk regression means
    mu_calf[t] <- beta0_calf +
      beta1_calf_wolf * wolf_N_tot_std[t] +
      beta2_calf_winterppt * winter_ppt_mm[t]
    
    mu_ya[t] <- beta0_ya +
      beta1_ya_wolf * wolf_N_tot_std[t] +
      beta2_ya_winterppt * winter_ppt_mm[t]
    
    mu_oa[t] <- beta0_oa +
      beta1_oa_wolf * wolf_N_tot_std[t] +
      beta2_oa_winterppt * winter_ppt_mm[t]
    
    # wolf regression means
    mu_wp[t] <- beta0_wp +
      beta1_wp_elk * elk_N_female_std[t] +
      beta2_wp_bison * NR_Bison[t]
    
    mu_wa[t] <- beta0_wa +
      beta1_wa_elk * elk_N_female_std[t] +
      beta2_wa_bison * NR_Bison[t]
    
    # observed Stage 1 estimates on logit scale
    logit_elk_s_c[t] ~ dnorm(mu_calf[t], tau_calf)
    logit_elk_s_ya[t] ~ dnorm(mu_ya[t], tau_ya)
    logit_elk_s_oa[t] ~ dnorm(mu_oa[t], tau_oa)
    
    logit_wolf_s_p[t] ~ dnorm(mu_wp[t], tau_wp)
    logit_wolf_s_a[t] ~ dnorm(mu_wa[t], tau_wa)
  }
})

################################################################################
###########---------------- Constants and data ----------------#################
################################################################################

stage2_constants <- list(
  n_years = nrow(reg_df)
)

stage2_data <- list(
  # elk responses
  logit_elk_s_c = qlogis(reg_df$elk_s_c),
  logit_elk_s_ya = qlogis(reg_df$elk_s_ya),
  logit_elk_s_oa = qlogis(reg_df$elk_s_oa),
  
  # wolf responses
  logit_wolf_s_p = qlogis(reg_df$wolf_s_p),
  logit_wolf_s_a = qlogis(reg_df$wolf_s_a),
  
  # elk predictors
  wolf_N_tot_std = reg_df$wolf_N_tot_std,
  winter_ppt_mm = reg_df$winter_ppt_mm_std,
  
  # wolf predictors
  elk_N_female_std = reg_df$elk_N_female_std,
  NR_Bison = reg_df$NR_Bison_std
)

################################################################################
###########-------------------- Parameters --------------------#################
################################################################################

stage2_params <- c(
  # elk
  "beta0_calf", "beta1_calf_wolf", "beta2_calf_winterppt",
  "beta0_ya", "beta1_ya_wolf", "beta2_ya_winterppt",
  "beta0_oa", "beta1_oa_wolf", "beta2_oa_winterppt",
  "sigma_calf", "sigma_ya", "sigma_oa",
  
  # wolf
  "beta0_wp", "beta1_wp_elk", "beta2_wp_bison",
  "beta0_wa", "beta1_wa_elk", "beta2_wa_bison",
  "sigma_wp", "sigma_wa"
)

################################################################################
###########------------------- Run model ----------------------#################
################################################################################

set.seed(17)
iters = 100000
burn = 40000
nchains = 3
thin = 4

start_time <- Sys.time()

stage2_mod <- nimbleMCMC(
  code = stage2_code,
  data = stage2_data,
  constants = stage2_constants,
  monitors = stage2_params,
  nchains = nchains,
  niter = iters,
  nburnin = burn,
  thin = thin,
  summary = TRUE
)

end_time <- Sys.time()
run_time <- end_time - start_time

print(paste0('Model runtime: ',
             round(run_time, 2),
             ' ',
             units(run_time)))

# SAVE OUTPUT
# stop('The following line will overwrite data. Are you sure you would like to proceed?')
save.image('data/outputs/ICM_Regression_environment_2026-04-14.RData')

# load('data/outputs/ICM_Regression_environment_2026-04-14.RData')

################################################################################
###########--------------------- Results ----------------------#################
################################################################################

# 1) extract each chain as a plain matrix
mats <- lapply(stage2_mod$samples, function(ch) as.matrix(ch))

# 2) trim all chains to same length
lens <- sapply(mats, nrow)
L <- min(lens)

mats_trim <- lapply(mats, function(M) {
  tail(M, L)
})

# 3) rebuild clean mcmc.list
stage2_clean <- mcmc.list(lapply(mats_trim, function(M) {
  mcmc(M, start = 1, end = L, thin = 1)
}))

# 4) coefficient summaries
reg_coefs <- MCMCsummary(
  stage2_clean,
  params = stage2_params
)

round(reg_coefs, 2)

################################################################################
###########------------ Checking model convergence ------------#################
################################################################################

bad_reg_coefs <- data.frame(
  param = rownames(reg_coefs),
  Rhat = reg_coefs[, "Rhat"]
)

bad_reg_coefs <- bad_reg_coefs[bad_reg_coefs$Rhat > 1.1, ]
bad_reg_coefs

################################################################################
########---------- Elk regression plots + coefficient densities -------#########
################################################################################

plot_df_elk <- bind_rows(
  reg_df %>%
    transmute(
      year,
      elk_surv = elk_s_c,
      wolf_N_tot = wolf_N_tot,
      winter_ppt_mm = winter_ppt_mm_raw,
      stage = "Calf survival"
    ),
  reg_df %>%
    transmute(
      year,
      elk_surv = elk_s_ya,
      wolf_N_tot = wolf_N_tot,
      winter_ppt_mm = winter_ppt_mm_raw,
      stage = "Young adult survival"
    ),
  reg_df %>%
    transmute(
      year,
      elk_surv = elk_s_oa,
      wolf_N_tot = wolf_N_tot,
      winter_ppt_mm = winter_ppt_mm_raw,
      stage = "Old adult survival"
    )
)

post_mat <- do.call(rbind, lapply(stage2_clean, as.matrix))

################################################################################
########---------------- Wolf abundance effect on elk ----------------##########
################################################################################

x_grid_wolf <- seq(
  min(reg_df$wolf_N_tot, na.rm = TRUE),
  max(reg_df$wolf_N_tot, na.rm = TRUE),
  length.out = 200
)

wolf_mean <- mean(reg_df$wolf_N_tot, na.rm = TRUE)
wolf_sd <- sd(reg_df$wolf_N_tot, na.rm = TRUE)
x_grid_wolf_std <- (x_grid_wolf - wolf_mean) / wolf_sd

winter_ppt_mean <- mean(reg_df$winter_ppt_mm_std, na.rm = TRUE)

pred_calf_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_calf"] +
      post_mat[, "beta1_calf_wolf"] * x +
      post_mat[, "beta2_calf_winterppt"] * winter_ppt_mean
  )
})

pred_ya_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_ya"] +
      post_mat[, "beta1_ya_wolf"] * x +
      post_mat[, "beta2_ya_winterppt"] * winter_ppt_mean
  )
})

pred_oa_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_oa"] +
      post_mat[, "beta1_oa_wolf"] * x +
      post_mat[, "beta2_oa_winterppt"] * winter_ppt_mean
  )
})

line_df_wolf_elk <- bind_rows(
  data.frame(
    x = x_grid_wolf,
    elk_surv = apply(pred_calf_wolf, 2, mean),
    elk_low = apply(pred_calf_wolf, 2, quantile, probs = 0.025),
    elk_high = apply(pred_calf_wolf, 2, quantile, probs = 0.975),
    stage = "Calf survival"
  ),
  data.frame(
    x = x_grid_wolf,
    elk_surv = apply(pred_ya_wolf, 2, mean),
    elk_low = apply(pred_ya_wolf, 2, quantile, probs = 0.025),
    elk_high = apply(pred_ya_wolf, 2, quantile, probs = 0.975),
    stage = "Young adult survival"
  ),
  data.frame(
    x = x_grid_wolf,
    elk_surv = apply(pred_oa_wolf, 2, mean),
    elk_low = apply(pred_oa_wolf, 2, quantile, probs = 0.025),
    elk_high = apply(pred_oa_wolf, 2, quantile, probs = 0.975),
    stage = "Old adult survival"
  )
)

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
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Wolf abundance",
    y = "Elk survival",
    title = "Estimated effect of wolf abundance on elk survival",
    subtitle = "Winter precipitation held at its mean"
  )

wolf_plot_elk

################################################################################
#######------------- Winter precipitation effect on elk ----------------########
################################################################################

x_grid_ppt_raw <- seq(
  min(reg_df$winter_ppt_mm_raw, na.rm = TRUE),
  max(reg_df$winter_ppt_mm_raw, na.rm = TRUE),
  length.out = 200
)

ppt_mean_raw <- mean(reg_df$winter_ppt_mm_raw, na.rm = TRUE)
ppt_sd_raw <- sd(reg_df$winter_ppt_mm_raw, na.rm = TRUE)

x_grid_ppt_std <- (x_grid_ppt_raw - ppt_mean_raw) / ppt_sd_raw

wolf_mean_std <- mean(reg_df$wolf_N_tot_std, na.rm = TRUE)

pred_calf_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_calf"] +
      post_mat[, "beta1_calf_wolf"] * wolf_mean_std +
      post_mat[, "beta2_calf_winterppt"] * x
  )
})

pred_ya_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_ya"] +
      post_mat[, "beta1_ya_wolf"] * wolf_mean_std +
      post_mat[, "beta2_ya_winterppt"] * x
  )
})

pred_oa_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_oa"] +
      post_mat[, "beta1_oa_wolf"] * wolf_mean_std +
      post_mat[, "beta2_oa_winterppt"] * x
  )
})

line_df_ppt_elk <- bind_rows(
  data.frame(
    x = x_grid_ppt_raw,
    elk_surv = apply(pred_calf_ppt, 2, mean),
    elk_low = apply(pred_calf_ppt, 2, quantile, probs = 0.025),
    elk_high = apply(pred_calf_ppt, 2, quantile, probs = 0.975),
    stage = "Calf survival"
  ),
  data.frame(
    x = x_grid_ppt_raw,
    elk_surv = apply(pred_ya_ppt, 2, mean),
    elk_low = apply(pred_ya_ppt, 2, quantile, probs = 0.025),
    elk_high = apply(pred_ya_ppt, 2, quantile, probs = 0.975),
    stage = "Young adult survival"
  ),
  data.frame(
    x = x_grid_ppt_raw,
    elk_surv = apply(pred_oa_ppt, 2, mean),
    elk_low = apply(pred_oa_ppt, 2, quantile, probs = 0.025),
    elk_high = apply(pred_oa_ppt, 2, quantile, probs = 0.975),
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
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Winter precipitation",
    y = "Elk survival",
    title = "Estimated effect of winter precipitation on elk survival",
    subtitle = "Wolf abundance held at its mean"
  )

ppt_plot_elk

################################################################################
########---------- Wolf regression plots + coefficient densities -------########
################################################################################

plot_df_wolf <- bind_rows(
  reg_df %>%
    transmute(
      year,
      wolf_surv = wolf_s_p,
      elk_N_female = elk_N_female,
      NR_Bison = NR_Bison_raw,
      stage = "Pup survival"
    ),
  reg_df %>%
    transmute(
      year,
      wolf_surv = wolf_s_a,
      elk_N_female = elk_N_female,
      NR_Bison = NR_Bison_raw,
      stage = "Adult survival"
    )
)

################################################################################
########------------- Elk abundance effect on wolf survival -------------#######
################################################################################

x_grid_elk <- seq(
  min(reg_df$elk_N_female, na.rm = TRUE),
  max(reg_df$elk_N_female, na.rm = TRUE),
  length.out = 200
)

elk_mean <- mean(reg_df$elk_N_female, na.rm = TRUE)
elk_sd <- sd(reg_df$elk_N_female, na.rm = TRUE)
x_grid_elk_std <- (x_grid_elk - elk_mean) / elk_sd

bison_mean <- mean(reg_df$NR_Bison_std, na.rm = TRUE)

pred_wp_elk <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_wp"] +
      post_mat[, "beta1_wp_elk"] * x +
      post_mat[, "beta2_wp_bison"] * bison_mean
  )
})

pred_wa_elk <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_wa"] +
      post_mat[, "beta1_wa_elk"] * x +
      post_mat[, "beta2_wa_bison"] * bison_mean
  )
})

line_df_elk_wolf <- bind_rows(
  data.frame(
    x = x_grid_elk,
    wolf_surv = apply(pred_wp_elk, 2, mean),
    wolf_low = apply(pred_wp_elk, 2, quantile, probs = 0.025),
    wolf_high = apply(pred_wp_elk, 2, quantile, probs = 0.975),
    stage = "Pup survival"
  ),
  data.frame(
    x = x_grid_elk,
    wolf_surv = apply(pred_wa_elk, 2, mean),
    wolf_low = apply(pred_wa_elk, 2, quantile, probs = 0.025),
    wolf_high = apply(pred_wa_elk, 2, quantile, probs = 0.975),
    stage = "Adult survival"
  )
)

elk_plot_wolf <- ggplot(plot_df_wolf, aes(x = elk_N_female, y = wolf_surv)) +
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
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Elk female abundance",
    y = "Wolf survival",
    title = "Estimated effect of elk abundance on wolf survival",
    subtitle = "Bison abundance held at its mean"
  )

elk_plot_wolf

################################################################################
###########-------------- Bison effect on wolf survival ----------------########
################################################################################

x_grid_bison_raw <- seq(
  min(reg_df$NR_Bison_raw, na.rm = TRUE),
  max(reg_df$NR_Bison_raw, na.rm = TRUE),
  length.out = 200
)

bison_mean_raw <- mean(reg_df$NR_Bison_raw, na.rm = TRUE)
bison_sd_raw <- sd(reg_df$NR_Bison_raw, na.rm = TRUE)

x_grid_bison_std <- (x_grid_bison_raw - bison_mean_raw) / bison_sd_raw

# hold elk abundance constant at its mean standardized value
elk_mean_std <- mean(reg_df$elk_N_female_std, na.rm = TRUE)

pred_wp_bison <- sapply(x_grid_bison_std, function(x) {
  plogis(
    post_mat[, "beta0_wp"] +
      post_mat[, "beta1_wp_elk"] * elk_mean_std +
      post_mat[, "beta2_wp_bison"] * x
  )
})

pred_wa_bison <- sapply(x_grid_bison_std, function(x) {
  plogis(
    post_mat[, "beta0_wa"] +
      post_mat[, "beta1_wa_elk"] * elk_mean_std +
      post_mat[, "beta2_wa_bison"] * x
  )
})

line_df_bison_wolf <- bind_rows(
  data.frame(
    x = x_grid_bison_raw,
    wolf_surv = apply(pred_wp_bison, 2, mean),
    wolf_low = apply(pred_wp_bison, 2, quantile, probs = 0.025),
    wolf_high = apply(pred_wp_bison, 2, quantile, probs = 0.975),
    stage = "Pup survival"
  ),
  data.frame(
    x = x_grid_bison_raw,
    wolf_surv = apply(pred_wa_bison, 2, mean),
    wolf_low = apply(pred_wa_bison, 2, quantile, probs = 0.025),
    wolf_high = apply(pred_wa_bison, 2, quantile, probs = 0.975),
    stage = "Adult survival"
  )
)

bison_plot_wolf <- ggplot(plot_df_wolf, aes(x = NR_Bison, y = wolf_surv)) +
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
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Bison abundance",
    y = "Wolf survival",
    title = "Estimated effect of bison abundance on wolf survival",
    subtitle = "Elk abundance held at its mean"
  )

bison_plot_wolf

################################################################################
##########--------------- Posterior coefficient densities ----------------######
################################################################################

beta_df <- data.frame(
  # elk
  beta0_calf = post_mat[, "beta0_calf"],
  beta1_calf_wolf = post_mat[, "beta1_calf_wolf"],
  beta2_calf_winterppt = post_mat[, "beta2_calf_winterppt"],
  
  beta0_ya = post_mat[, "beta0_ya"],
  beta1_ya_wolf = post_mat[, "beta1_ya_wolf"],
  beta2_ya_winterppt = post_mat[, "beta2_ya_winterppt"],
  
  beta0_oa = post_mat[, "beta0_oa"],
  beta1_oa_wolf = post_mat[, "beta1_oa_wolf"],
  beta2_oa_winterppt = post_mat[, "beta2_oa_winterppt"],
  
  # wolf
  beta0_wp = post_mat[, "beta0_wp"],
  beta1_wp_elk = post_mat[, "beta1_wp_elk"],
  beta2_wp_bison = post_mat[, "beta2_wp_bison"],
  
  beta0_wa = post_mat[, "beta0_wa"],
  beta1_wa_elk = post_mat[, "beta1_wa_elk"],
  beta2_wa_bison = post_mat[, "beta2_wa_bison"]
)

beta_long <- beta_df %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

beta_long$parameter <- factor(
  beta_long$parameter,
  levels = c(
    "beta0_calf", "beta1_calf_wolf", "beta2_calf_winterppt",
    "beta0_ya", "beta1_ya_wolf", "beta2_ya_winterppt",
    "beta0_oa", "beta1_oa_wolf", "beta2_oa_winterppt",
    "beta0_wp", "beta1_wp_elk", "beta2_wp_bison",
    "beta0_wa", "beta1_wa_elk", "beta2_wa_bison"
  ),
  labels = c(
    "Calf intercept", "Calf wolf effect", "Calf winter ppt effect",
    "YA intercept", "YA wolf effect", "YA winter ppt effect",
    "OA intercept", "OA wolf effect", "OA winter ppt effect",
    "Wolf pup intercept", "Wolf pup elk effect", "Wolf pup bison effect",
    "Wolf adult intercept", "Wolf adult elk effect", "Wolf adult bison effect"
  )
)

coef_plot <- ggplot(beta_long, aes(x = value)) +
  geom_density(fill = "#236192", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  theme_classic() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distributions of regression coefficients"
  )

coef_plot

################################################################################
#########------------------------ Combine plots -------------------------#######
################################################################################

elk_panel <- plot_grid(
  wolf_plot_elk,
  ppt_plot_elk,
  ncol = 2,
  labels = c("A", "B")
)

wolf_panel <- plot_grid(
  elk_plot_wolf,
  bison_plot_wolf,
  ncol = 2,
  labels = c("C", "D")
)

top_row <- plot_grid(
  elk_panel,
  wolf_panel,
  ncol = 1
)

final_plot <- plot_grid(
  top_row,
  coef_plot,
  ncol = 1,
  rel_heights = c(1,1),
  labels = c("", "E")
)

top_row

coef_plot

final_plot

