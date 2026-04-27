### Integrated Community Model (ICM)
### Combines elk IPM + wolf IPM in one NIMBLE model
### Last updated: Apr. 16, 2026
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

#
##
###
#### # select whether you want to use the wolf population from the Northern Range:
#### # "NR"
#### # or whether you would like to include interior wolves:
#### # "full"
wolf_range <- "NR"
###
##
#

################################################################################
############---------------- Load elk data ----------------#####################
################################################################################

load("data/elk_adultSurvival_cjsMatrices.rData")

# immediately rename to avoid collisions
elk_z <- z
elk_y <- y
elk_is_class1 <- is_class1
elk_is_class2 <- is_class2
elk_first_seen <- first_seen

rm(z, y, is_class1, is_class2, first_seen)

elk_dat_n <- read.csv("data/elk_abundanceEstimates_stages.csv")
elk_dat_fec <- read.csv("data/elk_fecundity.csv")

# fix latent histories
elk_z_fixed <- elk_z
for (i in 1:nrow(elk_z)) {
  first_det <- which(elk_z[i, ] == 1)[1]
  last_det <- max(which(elk_z[i, ] == 1))
  if (!is.na(first_det) && !is.na(last_det) && first_det < last_det) {
    elk_z_fixed[i, first_det:last_det] <- 1
  }
}
elk_z <- elk_z_fixed
rm(elk_z_fixed)

# align years
elk_shared_years <- intersect(as.numeric(elk_dat_n$year), as.numeric(colnames(elk_z)))
elk_dat_n <- elk_dat_n[elk_dat_n$year %in% elk_shared_years, ]
elk_dat_fec <- elk_dat_fec[elk_dat_fec$year %in% elk_shared_years, ]

colnames(elk_is_class1) <- colnames(elk_z)
colnames(elk_is_class2) <- colnames(elk_z)

elk_z <- elk_z[, colnames(elk_z) %in% elk_shared_years, drop = FALSE]
elk_y <- elk_y[, colnames(elk_y) %in% elk_shared_years, drop = FALSE]
elk_is_class1 <- elk_is_class1[, colnames(elk_is_class1) %in% elk_shared_years, drop = FALSE]
elk_is_class2 <- elk_is_class2[, colnames(elk_is_class2) %in% elk_shared_years, drop = FALSE]

# female abundance
elk_dat_n$n_female <- elk_dat_n$n_cow + (elk_dat_n$n_calf / 2)

elk_n_years <- nrow(elk_dat_n)
elk_N_indiv <- nrow(elk_y)

################################################################################
############--------------- Load wolf data ----------------#####################
################################################################################

load("data/wolf_adultSurvival_cjsMatrices.rData")

# immediately rename to avoid collisions
wolf_z <- z
wolf_y <- y
wolf_is_class1 <- is_class1
wolf_is_class2 <- is_class2
wolf_first_seen <- first_seen

rm(z, y, is_class1, is_class2, first_seen)

wolf_pop <- if (wolf_range == "NR") {
  read.csv("data/wolf_nr_pop.csv")
} else if (wolf_range == "full") {
  read.csv("data/wolf_full_park_pop.csv")
} else {
  stop("wolf_range must be 'NR' or 'full'")
}

# align years
wolf_shared_years <- intersect(as.numeric(wolf_pop$seasonal.year), as.numeric(colnames(wolf_z)))
wolf_pop <- wolf_pop[wolf_pop$seasonal.year %in% wolf_shared_years, ]

colnames(wolf_is_class1) <- colnames(wolf_z)
colnames(wolf_is_class2) <- colnames(wolf_z)

wolf_z <- wolf_z[, colnames(wolf_z) %in% wolf_shared_years, drop = FALSE]
wolf_y <- wolf_y[, colnames(wolf_y) %in% wolf_shared_years, drop = FALSE]
wolf_is_class1 <- wolf_is_class1[, colnames(wolf_is_class1) %in% wolf_shared_years, drop = FALSE]
wolf_is_class2 <- wolf_is_class2[, colnames(wolf_is_class2) %in% wolf_shared_years, drop = FALSE]

# rebuild first_seen after trimming years
wolf_first_seen <- apply(wolf_y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) ncol(wolf_y) else first
})
wolf_first_seen <- as.integer(wolf_first_seen)

wolf_n_years <- nrow(wolf_pop)
wolf_N_indiv <- nrow(wolf_y)

################################################################################
########---------- Align elk and wolf data to same years ---------##############
################################################################################

# choose years shared by BOTH elk and wolf submodels
community_years <- intersect(elk_shared_years, wolf_shared_years)

# trim elk to community years
elk_dat_n <- elk_dat_n[elk_dat_n$year %in% community_years, ]
elk_dat_fec <- elk_dat_fec[elk_dat_fec$year %in% community_years, ]

elk_z <- elk_z[, colnames(elk_z) %in% community_years, drop = FALSE]
elk_y <- elk_y[, colnames(elk_y) %in% community_years, drop = FALSE]
elk_is_class1 <- elk_is_class1[, colnames(elk_is_class1) %in% community_years, drop = FALSE]
elk_is_class2 <- elk_is_class2[, colnames(elk_is_class2) %in% community_years, drop = FALSE]

# trim wolf to community years
wolf_pop <- wolf_pop[wolf_pop$seasonal.year %in% community_years, ]

wolf_z <- wolf_z[, colnames(wolf_z) %in% community_years, drop = FALSE]
wolf_y <- wolf_y[, colnames(wolf_y) %in% community_years, drop = FALSE]
wolf_is_class1 <- wolf_is_class1[, colnames(wolf_is_class1) %in% community_years, drop = FALSE]
wolf_is_class2 <- wolf_is_class2[, colnames(wolf_is_class2) %in% community_years, drop = FALSE]

wolf_first_seen <- apply(wolf_y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) ncol(wolf_y) else first
})
wolf_first_seen <- as.integer(wolf_first_seen)

# update year counts after trimming
elk_n_years <- nrow(elk_dat_n)
wolf_n_years <- nrow(wolf_pop)

if (elk_n_years != wolf_n_years) {
  stop("Elk and wolf year counts still do not match after alignment.")
}

n_years <- elk_n_years


################################################################################
############-------------- Load in covariates -------------#####################
################################################################################

# then covariates
annual_prism <- read.csv('data/covariates/prism_annual_precip_tmean.csv')
bison <- read.csv('data/covariates/NR_Bison_Abundance.csv') %>%
  rename(year = Year)
grizzly <- read.csv('data/covariates/grizzly_abundances_cleaned.csv') %>%
  rename(griz_N = N) %>%
  select(year, griz_N)

# bind data into one dataframe 'covars"
covars <- left_join(annual_prism, bison)
covars <- left_join(covars, grizzly)

# trim covariate data to shared years 
covars <- covars[covars$year %in% elk_dat_n$year,]

# standardize all covariates except year
covars_std <- covars %>%
  mutate(
    across(-year, ~ as.numeric(scale(.)))
  )

################################################################################
############---------------- ICM NIMBLE Code --------------#####################
################################################################################

icm_code <- nimbleCode({
  
  ##########---------------- ELK SUBMODEL ----------------#############
  
  # elk priors
  for (t in 1:n_years) {
    
    # logit(elk_s_c[t]) ~ dnorm(qlogis(0.22), 1 / 0.5^2)
    # logit(elk_s_ya[t]) ~ dnorm(qlogis(0.90), 1 / 0.5^2) # mean survival = 0.9
    # logit(elk_s_oa[t]) ~ dnorm(qlogis(0.80), 1 / 0.5^2) # mean survival = 0.8
    
    logit(elk_p_13[t]) ~ dnorm(qlogis(0.15), 1 / 0.5^2)
  }
  
  for (t in 1:(n_years - 1)) {
    elk_f_ya[t] ~ dbeta(1, 1) # assume elk cows can't have more than 1 calf each (f is bounded 0-1)
    elk_f_oa[t] ~ dbeta(1, 1)
  }
  
  elk_sigma_obs_female ~ dunif(0.05, 2)
  elk_tau_obs_female <- 1 / (elk_sigma_obs_female^2)
  
  # elk initial values
  elk_N_1y[1] ~ dunif(0, 5000) 
  elk_N_ya[1] ~ dunif(0, 15000)
  elk_N_oa[1] ~ dunif(0, 5000) 
  
  elk_N_female[1] <- elk_N_1y[1] + elk_N_ya[1] + elk_N_oa[1]
  elk_obs_female[1] ~ dlnorm(log(elk_N_female[1] + 1e-6), elk_tau_obs_female)
  
  # elk state-space model
  for (t in 1:(n_years - 1)) {
    
    elk_mu_1y[t + 1] <- elk_f_ya[t] * elk_s_c[t] * elk_N_ya[t] +
      elk_f_oa[t] * elk_s_c[t] * elk_N_oa[t]
    
    elk_mu_ya[t + 1] <- elk_s_ya[t] * elk_N_1y[t] +
      elk_s_ya[t] * (1 - elk_p_13[t]) * elk_N_ya[t]
    
    elk_mu_oa[t + 1] <- elk_s_ya[t] * elk_p_13[t] * elk_N_ya[t] +
      elk_s_oa[t] * elk_N_oa[t]
    
    elk_N_1y[t + 1] ~ dpois(max(1e-6, elk_mu_1y[t + 1]))
    elk_N_ya[t + 1] ~ dpois(max(1e-6, elk_mu_ya[t + 1]))
    elk_N_oa[t + 1] ~ dpois(max(1e-6, elk_mu_oa[t + 1]))
    
    elk_N_female[t + 1] <- elk_N_1y[t + 1] + elk_N_ya[t + 1] + elk_N_oa[t + 1]
    elk_obs_female[t + 1] ~ dlnorm(log(elk_N_female[t + 1] + 1e-6), elk_tau_obs_female)
  }
  
  # elk CJS
  for (t in 1:n_years) {
    elk_p_det[t] ~ dunif(0, 1)
  }
  
  for (i in 1:elk_N_indiv) {
    elk_z[i, 1] ~ dbern(equals(1, elk_first_seen[i]))
    
    for (t in 2:n_years) {
      elk_phi[i, t] <- elk_is_class1[i, t - 1] * elk_s_ya[t - 1] +
        elk_is_class2[i, t - 1] * elk_s_oa[t - 1] +
        1e-10
      
      elk_is_first[i, t] <- equals(t, elk_first_seen[i])
      elk_after_first[i, t] <- step(t - elk_first_seen[i] - 0.5)
      elk_survive_to_t[i, t] <- elk_z[i, t - 1] * elk_phi[i, t]
      
      elk_gamma[i, t] <- elk_is_first[i, t] +
        elk_after_first[i, t] * (1 - elk_is_first[i, t]) * elk_survive_to_t[i, t]
      
      elk_z[i, t] ~ dbern(elk_gamma[i, t])
    }
    
    for (t in 1:n_years) {
      elk_y[i, t] ~ dbern(elk_p_det[t] * elk_z[i, t])
    }
  }
  
  # elk calf survival
  for (t in 1:(n_years - 1)) {
    elk_CCR_prob_young[t] <- elk_f_ya[t] * elk_s_c[t] # elk calf survival is the same between
    elk_CCR_prob_old[t] <- elk_f_oa[t] * elk_s_c[t]   # young and old cows
    
    elk_CCR_c_fromYoungCows[t] ~ dbin(elk_CCR_prob_young[t], elk_CCR_cow_youngadult[t])
    elk_CCR_c_fromOldCows[t] ~ dbin(elk_CCR_prob_old[t], elk_CCR_cow_oldadult[t])
    
    elk_CCR_c[t] <- elk_CCR_c_fromYoungCows[t] + elk_CCR_c_fromOldCows[t]
  }
  
  # elk pregnancy
  for (t in 1:(n_years - 1)) {
    elk_young_num_preg[t] ~ dbin(elk_f_ya[t], elk_young_num_capt[t])
    elk_old_num_preg[t] ~ dbin(elk_f_oa[t], elk_old_num_capt[t])
  }
  
  # elk growth
  for (t in 1:n_years) {
    elk_harvested_13yo[t] ~ dbin(elk_p_13[t], elk_harvested_ya[t])
  }
  
  ##########---------------- WOLF SUBMODEL ----------------#############
  
  # wolf priors
  for (t in 1:n_years) {
    # logit(wolf_s_p[t]) ~ dnorm(qlogis(0.5), 1 / 0.5^2) # mean = 0.5
    # logit(wolf_s_a[t]) ~ dnorm(qlogis(0.9), 1 / 0.5^2) # mean = 0.9
    wolf_f[t] ~ dgamma(2, 2)   # mean = 1, but tighter than gamma(1,1)
  }
  
  wolf_sigma_obs ~ dunif(0.05, 2)
  wolf_tau_obs <- 1 / (wolf_sigma_obs^2)
  
  # because this is a reintroduced population, the initial values are known.
  wolf_N_a[1] <- 8 # 14 wolves originally introduced the first year:
  wolf_N_p[1] <- 6 # 8 adults, 6 pups
  wolf_N_p_sum[1] <- 0 
  # there were no summer pups the first year because wolves were introduced in Jan 1995, 
  # and then the pups counted in summer 1995 would belong to the next year's cohort
  
  # wolf abundances and observation errors
  for (t in 1:(n_years-1)) {
    
    # expected adults in Dec of year t+1 are surviving adults from prev. year + pups
    wolf_mu_a[t + 1] <- wolf_s_a[t] * wolf_N_a[t] + wolf_N_p[t] + equals(t, 1) * 8 # 8 adults introduced in 1996
    wolf_N_a[t + 1] ~ dpois(max(1e-6, wolf_mu_a[t + 1]))
    
    # pups counted in Dec of year t produce pups in summer of t+1
    wolf_N_p_sum[t + 1] ~ dpois(max(1e-6, wolf_f[t] * wolf_N_a[t]))
  }
  
  # observation processes 
  for (t in 1:n_years){
    
    # total wolves in Dec of year t+1
    wolf_N_tot[t] <- wolf_N_a[t] + wolf_N_p[t]
    
    # observed summer pups in year t+1
    wolf_obs_p_sum[t] ~ dpois(max(1e-6, wolf_N_p_sum[t]))
    
    # december count observations
    wolf_obs_tot[t] ~ dlnorm(log(wolf_N_tot[t] + 1e-6), wolf_tau_obs)
    wolf_obs_p[t] ~ dlnorm(log(wolf_N_p[t] + 1e-6), wolf_tau_obs)
    wolf_obs_a[t] ~ dlnorm(log(wolf_N_a[t] + 1e-6), wolf_tau_obs)
  }
  
  # pup survival 
  for (t in 2:n_years) {
    wolf_N_p_bio[t] ~ dbin(wolf_s_p[t - 1], wolf_N_p_sum[t]) # pups from the biological process (reproduction + survival)
    wolf_N_p[t] <- wolf_N_p_bio[t] + equals(t, 2) * 9 # add 9 pups to the year 1996 for the introduced pups
  }
  
  # wolf CJS
  for (t in 1:n_years) {
    wolf_p_det[t] ~ dunif(0, 1)
  }
  
  for (i in 1:wolf_N_indiv) {
    wolf_z[i, 1] ~ dbern(equals(1, wolf_first_seen[i]))
    
    for (t in 2:n_years) {
      wolf_phi[i, t] <- wolf_is_class1[i, t - 1] * wolf_s_p[t - 1] +
        wolf_is_class2[i, t - 1] * wolf_s_a[t - 1] +
        1e-10
      
      wolf_is_first[i, t] <- equals(t, wolf_first_seen[i])
      wolf_after_first[i, t] <- step(t - wolf_first_seen[i] - 0.5)
      wolf_survive_to_t[i, t] <- wolf_z[i, t - 1] * wolf_phi[i, t]
      
      wolf_gamma[i, t] <- wolf_is_first[i, t] +
        wolf_after_first[i, t] * (1 - wolf_is_first[i, t]) * wolf_survive_to_t[i, t]
      
      wolf_z[i, t] ~ dbern(wolf_gamma[i, t])
    }
    
    for (t in 1:n_years) {
      wolf_y[i, t] ~ dbern(wolf_p_det[t] * wolf_z[i, t])
    }
  }
  
  ##########---------------- REGRESSIONS ----------------#############
  
  # Priors for elk regression coefficients
  beta0_calfSurv ~ dnorm(qlogis(0.22), 1 / 0.3^2) # mean = 0.22
  beta1_calfSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_calfSurv_wintPPT ~ dnorm(0, 1 / 0.3^2)
  beta3_calfSurv_grizN ~ dnorm(0, 1 / 0.3^2)
  
  beta0_yaSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2) # mean = 0.90
  beta1_yaSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_yaSurv_wintPPT ~ dnorm(0, 1 / 0.3^2)
  beta3_yaSurv_grizN ~ dnorm(0, 1 / 0.3^2)
  
  beta0_oaSurv ~ dnorm(qlogis(0.80), 1 / 0.3^2) # mean = 0.80
  beta1_oaSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_oaSurv_wintPPT ~ dnorm(0, 1 / 0.3^2)
  beta3_oaSurv_grizN ~ dnorm(0, 1 / 0.3^2)
  
  # Priors for elk random year-effect SDs
  sigma_calf ~ dunif(0, 0.5)
  tau_calf <- 1 / (sigma_calf^2)
  
  sigma_ya ~ dunif(0, 0.3)
  tau_ya <- 1 / (sigma_ya^2)
  
  sigma_oa ~ dunif(0, 0.3)
  tau_oa <- 1 / (sigma_oa^2)
  
  # Priors for wolf regression coefficients
  beta0_wpupSurv ~ dnorm(qlogis(0.5), 1 / 0.3^2) # mean = 0.22
  beta1_wpupSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  beta2_wpupSurv_bisonN ~ dnorm(0, 1 / 0.3^2)
  
  beta0_wadSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2) # mean = 0.90
  beta1_wadSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  beta2_wadSurv_bisonN ~ dnorm(0, 1 / 0.3^2)
  
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
    
    # wolf random year effects
    eps_wolf_s_p[t] ~ dnorm(0, tau_wpup)
    eps_wolf_s_a[t] ~ dnorm(0, tau_wad)
    
    # standardize wolf abundance
    wolf_N_tot_std[t] <- (wolf_N_tot[t] - wolf_tot_mean) / wolf_tot_sd
    
    # standardize elk abundance
    elk_N_female_std[t] <- (elk_N_female[t] - elk_N_female_mean) / elk_N_female_sd
    
    # elk regression models
    logit(elk_s_c[t])  <- 
      beta0_calfSurv + 
      beta1_calfSurv_wolfN * wolf_N_tot_std[t] + 
      beta2_calfSurv_wintPPT * wintPPT[t] + 
      beta3_calfSurv_grizN * grizN_std[t] +
      eps_elk_s_c[t]
    
    logit(elk_s_ya[t]) <- 
      beta0_yaSurv + 
      beta1_yaSurv_wolfN * wolf_N_tot_std[t] + 
      beta2_yaSurv_wintPPT * wintPPT[t] + 
      beta3_yaSurv_grizN * grizN_std[t] +
      eps_elk_s_ya[t]
    
    logit(elk_s_oa[t]) <- 
      beta0_oaSurv + 
      beta1_oaSurv_wolfN * wolf_N_tot_std[t] + 
      beta2_oaSurv_wintPPT * wintPPT[t] + 
      beta3_oaSurv_grizN * grizN_std[t] +
      eps_elk_s_oa[t]
    
    # wolf regression models
    logit(wolf_s_p[t])  <- 
      beta0_wpupSurv + 
      beta1_wpupSurv_elkN * elk_N_female_std[t] + 
      beta2_wpupSurv_bisonN * bisonN_std[t] +
      eps_wolf_s_p[t]
    
    logit(wolf_s_a[t]) <- 
      beta0_wadSurv + 
      beta1_wadSurv_elkN * elk_N_female_std[t] + 
      beta2_wadSurv_bisonN * bisonN_std[t] +
      eps_wolf_s_a[t]
  }
})

################################################################################
###########---------------- Constants and data ----------------#################
################################################################################

icm_constants <- list(
  n_years = n_years,
  elk_N_indiv = elk_N_indiv,
  wolf_N_indiv = wolf_N_indiv,
  # these are added to allow the abundance metrics to be standardized
  wolf_tot_mean = mean(wolf_pop$total_abundance, na.rm = TRUE),
  wolf_tot_sd = sd(wolf_pop$total_abundance, na.rm = TRUE),
  elk_N_female_mean = mean(elk_dat_n$n_female, na.rm=TRUE),
  elk_N_female_sd =sd(elk_dat_n$n_female, na.rm=TRUE)
)

icm_data <- list(
  # elk state-space
  elk_obs_female = elk_dat_n$n_female,
  
  # elk CJS
  elk_y = elk_y,
  elk_is_class1 = elk_is_class1,
  elk_is_class2 = elk_is_class2,
  elk_first_seen = elk_first_seen,
  
  # elk fecundity
  elk_young_num_preg = elk_dat_fec$young_num_preg,
  elk_young_num_capt = elk_dat_fec$young_num_capt,
  elk_old_num_preg = elk_dat_fec$old_num_preg,
  elk_old_num_capt = elk_dat_fec$old_num_capt,
  
  # elk calf survival
  elk_CCR_cow_youngadult = elk_dat_fec$n_cows_young,
  elk_CCR_cow_oldadult = elk_dat_fec$n_cows_old,
  
  # elk growth
  elk_harvested_13yo = elk_dat_fec$harvested_age13,
  elk_harvested_ya = elk_dat_fec$harvested_total,
  
  # wolf state-space
  wolf_obs_tot = wolf_pop$total_abundance,
  wolf_obs_p_sum = wolf_pop$summer_pups,
  wolf_obs_p = wolf_pop$dec_pups,
  wolf_obs_a = wolf_pop$dec_adults,
  
  # wolf CJS
  wolf_y = wolf_y,
  wolf_is_class1 = wolf_is_class1,
  wolf_is_class2 = wolf_is_class2,
  wolf_first_seen = wolf_first_seen,
  
  # covariates
  wintPPT = covars_std$winter_ppt_mm,
  bisonN_std = covars_std$NR_Bison,
  grizN_std = covars_std$griz_N
)

################################################################################
###########----------------- Initial values -------------------#################
################################################################################

make_icm_inits <- function() {
  
  ## elk inits
  elk_init_N1y <- ifelse(
    is.na(elk_dat_n$n_calf),
    pmax(1, round(mean(elk_dat_n$n_calf, na.rm = TRUE))),
    pmax(1, round(elk_dat_n$n_calf))
  )
  
  elk_init_Nya <- ifelse(
    is.na(elk_dat_n$n_cow_youngadult),
    pmax(1, round(mean(elk_dat_n$n_cow_youngadult, na.rm = TRUE))),
    pmax(1, round(elk_dat_n$n_cow_youngadult))
  )
  
  elk_init_Noa <- ifelse(
    is.na(elk_dat_n$n_cow_oldadult),
    pmax(1, round(mean(elk_dat_n$n_cow_oldadult, na.rm = TRUE))),
    pmax(1, round(elk_dat_n$n_cow_oldadult))
  )
  
  elk_z_init <- matrix(NA, nrow = nrow(elk_y), ncol = ncol(elk_y))
  for (i in 1:nrow(elk_y)) {
    detections <- which(elk_y[i, ] == 1)
    if (length(detections) > 0) {
      first_det <- min(detections)
      last_det <- max(detections)
      elk_z_init[i, first_det:last_det] <- 1L
      if (first_det > 1) {
        elk_z_init[i, 1:(first_det - 1)] <- 0L
      }
    }
  }
  
  ## wolf inits
  wolf_init_Ntot <- ifelse(
    is.na(wolf_pop$total_abundance),
    max(1, round(mean(wolf_pop$total_abundance, na.rm = TRUE))),
    pmax(1, round(wolf_pop$total_abundance))
  )
  
  wolf_init_Np <- ifelse(
    is.na(wolf_pop$summer_pups),
    max(1, round(mean(wolf_pop$summer_pups, na.rm = TRUE) * 0.5)),
    pmax(1, round(wolf_pop$summer_pups * 0.5))
  )
  
  wolf_init_Na <- pmax(1, round(wolf_init_Ntot - wolf_init_Np))
  
  wolf_z_init <- matrix(NA, nrow = nrow(wolf_y), ncol = ncol(wolf_y))
  for (i in 1:nrow(wolf_y)) {
    detections <- which(wolf_y[i, ] == 1)
    if (length(detections) > 0) {
      first_det <- min(detections)
      last_det <- max(detections)
      
      wolf_z_init[i, first_det:last_det] <- 1L
      
      if (first_det > 1) {
        wolf_z_init[i, 1:(first_det - 1)] <- 0L
      }
      
      if (last_det < ncol(wolf_y)) {
        wolf_z_init[i, (last_det + 1):ncol(wolf_y)] <- 0L
      }
    }
  }
  
  list(
    # elk
    # elk_s_c = rep(0.22, n_years), # these get removed when the regressions are there
    # elk_s_ya = rep(0.90, n_years), # because they are now deterministic
    # elk_s_oa = rep(0.80, n_years),
    elk_p_13 = rep(0.15, n_years),
    elk_f_ya = rep(0.76, n_years - 1),
    elk_f_oa = rep(0.64, n_years - 1),
    elk_sigma_obs_female = 0.30,
    elk_N_1y = pmax(1, elk_init_N1y),
    elk_N_ya = pmax(1, elk_init_Nya),
    elk_N_oa = pmax(1, elk_init_Noa),
    elk_p_det = runif(n_years, 0.6, 0.95),
    elk_z = elk_z_init,
    
    # wolf
    # wolf_s_p = rep(0.5, n_years),
    # wolf_s_a = rep(0.9, n_years),
    wolf_f = rep(1.0, n_years),
    wolf_sigma_obs = 0.2,
    wolf_N_p_sum = pmax(1, round(wolf_pop$summer_pups)),
    wolf_N_p = wolf_init_Np,
    wolf_N_a = wolf_init_Na,
    wolf_p_det = runif(n_years, 0.6, 0.95),
    wolf_z = wolf_z_init
  )
}

################################################################################
###########-------------------- Parameters --------------------#################
################################################################################

icm_params <- c(
  # elk
  "elk_s_c", "elk_s_ya", "elk_s_oa", "elk_p_13",
  "elk_f_ya", "elk_f_oa",
  "elk_p_det",
  "elk_N_1y", "elk_N_ya", "elk_N_oa", "elk_N_female",
  
  # wolf
  "wolf_s_p", "wolf_s_a", "wolf_f",
  "wolf_p_det",
  "wolf_N_p_sum", "wolf_N_p", "wolf_N_a", "wolf_N_tot",
  
  # elk regression coefficients
  'beta0_calfSurv', 'beta1_calfSurv_wolfN', 'beta2_calfSurv_wintPPT', 'beta3_calfSurv_grizN',
  'beta0_yaSurv', 'beta1_yaSurv_wolfN', 'beta2_yaSurv_wintPPT', 'beta3_yaSurv_grizN',
  'beta0_oaSurv', 'beta1_oaSurv_wolfN', 'beta2_oaSurv_wintPPT', 'beta3_oaSurv_grizN',
  
  "sigma_calf", "sigma_ya", "sigma_oa",
  "eps_elk_s_c", "eps_elk_s_ya", "eps_elk_s_oa",
  
  # wolf regression coefficients
  'beta0_wpupSurv', 'beta1_wpupSurv_elkN', 'beta2_wpupSurv_bisonN',
  'beta0_wadSurv', 'beta1_wadSurv_elkN', 'beta2_wadSurv_bisonN',
  
  "sigma_wpup", "sigma_wad",
  "eps_wolf_s_p", "eps_wolf_s_a"
)

################################################################################
###########------------------- Run model ----------------------#################
################################################################################

set.seed(17)
nc <- 3
ni <- 200000
nb <- 40000
th <- 4

start_time <- Sys.time()

icm_mod <- nimbleMCMC(
  code = icm_code,
  data = icm_data,
  constants = icm_constants,
  inits = make_icm_inits,
  monitors = icm_params,
  nchains = nc,
  niter = ni,
  nburnin = nb,
  thin = th,
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
save.image('data/outputs/ICM_environment_2026-04-27.RData')

# load('data/outputs/ICM_environment_2026-04-27.RData')

################################################################################
############------- Reload packages if data loaded in --------##################
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
    'beta0_calfSurv', 'beta1_calfSurv_wolfN', 'beta2_calfSurv_wintPPT', 'beta3_calfSurv_grizN',
    'beta0_yaSurv', 'beta1_yaSurv_wolfN', 'beta2_yaSurv_wintPPT', 'beta3_yaSurv_grizN',
    'beta0_oaSurv', 'beta1_oaSurv_wolfN', 'beta2_oaSurv_wintPPT', 'beta3_oaSurv_grizN',
    # wolf reg coefs
    'beta0_wpupSurv', 'beta1_wpupSurv_elkN', 'beta2_wpupSurv_bisonN',
    'beta0_wadSurv', 'beta1_wadSurv_elkN', 'beta2_wadSurv_bisonN'
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
#########---------- Regression plots + coefficient densities ---------##########
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

# join elk, wolf, winter ppt, and grizzly estimates
plot_df_elk <- elk_surv_pts %>%
  left_join(wolf_pts, by = "year") %>%
  left_join(
    covars %>%
      select(year, winter_ppt_mm, griz_N),
    by = "year"
  )

# pull posterior draws from cleaned chains
post_mat <- do.call(rbind, lapply(icm_clean, as.matrix))

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

# hold winter precipitation and grizzly abundance constant at mean standardized values
winterPPT_mean_std <- mean(covars_std$winter_ppt_mm, na.rm = TRUE)
grizN_mean_std <- mean(covars_std$griz_N, na.rm = TRUE)

# fitted values for calf survival
pred_calf_wolf <- sapply(x_grid_wolf_std, function(x) {
  plogis(
    post_mat[, "beta0_calfSurv"] +
      post_mat[, "beta1_calfSurv_wolfN"] * x +
      post_mat[, "beta2_calfSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_calfSurv_grizN"] * grizN_mean_std
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
      post_mat[, "beta3_yaSurv_grizN"] * grizN_mean_std
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
      post_mat[, "beta3_oaSurv_grizN"] * grizN_mean_std
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
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Wolf abundance",
    y = "Elk survival",
    title = "Estimated effect of wolf abundance on elk survival",
    subtitle = "Winter precipitation and grizzly abundance held at their means"
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

# hold wolf abundance and grizzly abundance constant at mean standardized values
wolf_mean_std <- mean((wolf_pts$wolf_N_tot - icm_constants$wolf_tot_mean) / icm_constants$wolf_tot_sd,
                      na.rm = TRUE)

pred_calf_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_calfSurv"] +
      post_mat[, "beta1_calfSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_calfSurv_wintPPT"] * x +
      post_mat[, "beta3_calfSurv_grizN"] * grizN_mean_std
  )
})

pred_ya_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_yaSurv"] +
      post_mat[, "beta1_yaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_yaSurv_wintPPT"] * x +
      post_mat[, "beta3_yaSurv_grizN"] * grizN_mean_std
  )
})

pred_oa_ppt <- sapply(x_grid_ppt_std, function(x) {
  plogis(
    post_mat[, "beta0_oaSurv"] +
      post_mat[, "beta1_oaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_oaSurv_wintPPT"] * x +
      post_mat[, "beta3_oaSurv_grizN"] * grizN_mean_std
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
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Winter precipitation",
    y = "Elk survival",
    title = "Estimated effect of winter precipitation on elk survival",
    subtitle = "Wolf and grizzly abundance held at their means"
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

# hold wolf abundance and winter precipitation constant at mean standardized values
pred_calf_griz <- sapply(x_grid_griz_std, function(x) {
  plogis(
    post_mat[, "beta0_calfSurv"] +
      post_mat[, "beta1_calfSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_calfSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_calfSurv_grizN"] * x
  )
})

pred_ya_griz <- sapply(x_grid_griz_std, function(x) {
  plogis(
    post_mat[, "beta0_yaSurv"] +
      post_mat[, "beta1_yaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_yaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_yaSurv_grizN"] * x
  )
})

pred_oa_griz <- sapply(x_grid_griz_std, function(x) {
  plogis(
    post_mat[, "beta0_oaSurv"] +
      post_mat[, "beta1_oaSurv_wolfN"] * wolf_mean_std +
      post_mat[, "beta2_oaSurv_wintPPT"] * winterPPT_mean_std +
      post_mat[, "beta3_oaSurv_grizN"] * x
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
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  facet_wrap(~stage, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Grizzly abundance",
    y = "Elk survival",
    title = "Estimated effect of grizzly abundance on elk survival",
    subtitle = "Wolf abundance and winter precipitation held at their means"
  )

griz_plot_elk

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

wolf_surv_pts <- bind_rows(wolf_pup_pts, wolf_ad_pts)

wolf_plot_df <- wolf_surv_pts %>%
  left_join(elk_abund_pts, by = "year") %>%
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

# hold bison abundance constant at mean standardized value
bisonN_mean_std <- mean(covars_std$NR_Bison, na.rm = TRUE)

pred_wpup_elk <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_wpupSurv"] +
      post_mat[, "beta1_wpupSurv_elkN"] * x +
      post_mat[, "beta2_wpupSurv_bisonN"] * bisonN_mean_std
  )
})

pred_wad_elk <- sapply(x_grid_elk_std, function(x) {
  plogis(
    post_mat[, "beta0_wadSurv"] +
      post_mat[, "beta1_wadSurv_elkN"] * x +
      post_mat[, "beta2_wadSurv_bisonN"] * bisonN_mean_std
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
    subtitle = "Bison abundance held at its mean"
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

# hold elk abundance constant at mean standardized value
elk_mean_std <- mean((elk_abund_pts$elk_N_female - icm_constants$elk_N_female_mean) /
                       icm_constants$elk_N_female_sd,
                     na.rm = TRUE)

pred_wpup_bison <- sapply(x_grid_bison_std, function(x) {
  plogis(
    post_mat[, "beta0_wpupSurv"] +
      post_mat[, "beta1_wpupSurv_elkN"] * elk_mean_std +
      post_mat[, "beta2_wpupSurv_bisonN"] * x
  )
})

pred_wad_bison <- sapply(x_grid_bison_std, function(x) {
  plogis(
    post_mat[, "beta0_wadSurv"] +
      post_mat[, "beta1_wadSurv_elkN"] * elk_mean_std +
      post_mat[, "beta2_wadSurv_bisonN"] * x
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
    subtitle = "Elk abundance held at its mean"
  )

bison_plot_wolf

################################################################################
########--------------- Posterior coefficient densities ----------------########
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

coef_plot <- ggplot(beta_long, aes(x = value)) +
  geom_density(fill = "#236192", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~parameter, scales = "free", ncol = 4) +
  theme_classic() +
  labs(
    x = "Posterior value",
    y = "Density",
    title = "Posterior distributions of elk regression coefficients"
  )

coef_plot

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

wolf_coef_plot <- ggplot(wolf_beta_long, aes(x = value)) +
  geom_density(fill = "#6F263D", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
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
  coef_plot,
  ncol = 1,
  rel_heights = c(1.1, 1.1, 1.1, 1),
  labels = c("A", "B", "C", "D")
)

wolf_final_plot <- plot_grid(
  wolf_main_plot,
  bison_plot_wolf,
  wolf_coef_plot,
  ncol = 1,
  rel_heights = c(1.2, 1.2, 1),
  labels = c("E", "F", "G")
)

all_regression_plots <- plot_grid(
  elk_panel,
  wolf_final_plot,
  ncol = 2,
  rel_widths = c(1, 1)
)

all_regression_plots

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