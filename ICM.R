### Integrated Community Model (ICM)
### Combines elk IPM + wolf IPM in one NIMBLE model
### Last updated: Mar. 26, 2026
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
############---------------- ICM NIMBLE Code --------------#####################
################################################################################

icm_code <- nimbleCode({
  
  ##########---------------- ELK SUBMODEL ----------------#############
  
  # elk priors
  for (t in 1:n_years) {
    
    # logit(elk_s_c[t]) ~ dnorm(qlogis(0.22), 1 / 0.5^2) 
    # logit(elk_s_ya[t]) ~ dnorm(qlogis(0.90), 1 / 0.5^2) # mean survival = 0.9
    # logit(elk_s_oa[t]) ~ dnorm(qlogis(0.80), 1 / 0.5^2) # mean survival = 0.8
    
    # ^ These are commented out so the regressions at the bottom can run.
    
    logit(elk_p_13[t]) ~ dnorm(qlogis(0.15), 1 / 0.5^2)
  }
  
  for (t in 1:(n_years - 1)) {
    elk_f_ya[t] ~ dbeta(1, 1) # assume elk cows can't have more than 1 calf each (f is bounded 0-1)
    elk_f_oa[t] ~ dbeta(1, 1)
  }
  
  elk_sigma_obs_female ~ dunif(0.05, 2)
  elk_tau_obs_female <- 1 / (elk_sigma_obs_female^2)
  
  # elk initial values
  elk_N_1y[1] ~ dpois(2500) # mean values from rough estimations
  elk_N_ya[1] ~ dpois(10000) # using total abundance and classification (cow:calf) data
  elk_N_oa[1] ~ dpois(1000) # and hunter harvest (age) data for prop young vs. old adults
  
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
    logit(wolf_s_p[t]) ~ dnorm(qlogis(0.5), 1 / 0.5^2) # mean = 0.5
    logit(wolf_s_a[t]) ~ dnorm(qlogis(0.9), 1 / 0.5^2) # mean = 0.9
    wolf_f[t] ~ dgamma(2, 2)   # mean = 1, but tighter than gamma(1,1)
  }
  
  wolf_sigma_obs ~ dunif(0.05, 2)
  wolf_tau_obs <- 1 / (wolf_sigma_obs^2)
  
  wolf_N_a[1] ~ dpois(14) # 14 wolves originally introduced the first year
  
  for (t in 1:n_years) {
    
    # expected summer pups
    wolf_N_p_sum[t] ~ dpois(max(1e-6, wolf_f[t] * wolf_N_a[t]))

    # observed summer pup counts
    wolf_obs_p_sum[t] ~ dpois(max(1e-6, wolf_N_p_sum[t]))
    
    # summer pups surviving to December
    wolf_N_p[t] ~ dbin(wolf_s_p[t], wolf_N_p_sum[t])
    
    # December total
    wolf_N_tot[t] <- wolf_N_p[t] + wolf_N_a[t]
    wolf_obs_tot[t] ~ dlnorm(log(wolf_N_tot[t] + 1e-6), wolf_tau_obs)
  }
  
  for (t in 1:(n_years - 1)) {
    wolf_mu_a[t + 1] <- wolf_N_p[t] + wolf_s_a[t] * wolf_N_a[t]
    wolf_N_a[t + 1] ~ dpois(max(1e-6, wolf_mu_a[t + 1]))
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
  
  # Priors
  beta0_calfSurv ~ dnorm(qlogis(0.22), 1 / 0.5^2) # centered on 0.22
  beta0_yaSurv ~ dnorm(qlogis(0.90), 1 / 0.5^2) # centered on 0.90
  beta0_oaSurv ~ dnorm(qlogis(0.80), 1 / 0.5^2) # centered on 0.80
  
  beta1_calfSurv ~ dnorm(0, sd=sqrt(1/0.33))
  beta1_yaSurv ~ dnorm(0, sd=sqrt(1/0.33))
  beta1_oaSurv ~ dnorm(0, sd=sqrt(1/0.33))
  
  # Likelihood
  for (t in 1:n_years){
    
    # first standardize wolf abundance
    wolf_N_tot_std[t] <- (wolf_N_tot[t] - wolf_tot_mean) / wolf_tot_sd
    
    # then build regression models
    logit(elk_s_c[t]) <- beta0_calfSurv + beta1_calfSurv*wolf_N_tot_std[t] 
    logit(elk_s_ya[t]) <- beta0_yaSurv + beta1_yaSurv*wolf_N_tot_std[t] 
    logit(elk_s_oa[t]) <- beta0_oaSurv + beta1_oaSurv*wolf_N_tot_std[t] 
  }
})

################################################################################
###########---------------- Constants and data ----------------#################
################################################################################

icm_constants <- list(
  n_years = n_years,
  elk_N_indiv = elk_N_indiv,
  wolf_N_indiv = wolf_N_indiv,
  # these are added to allow the wolf abundance metric to be standardized
  wolf_tot_mean = mean(wolf_pop$total_abundance, na.rm = TRUE),
  wolf_tot_sd = sd(wolf_pop$total_abundance, na.rm = TRUE)
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
  
  # wolf CJS
  wolf_y = wolf_y,
  wolf_is_class1 = wolf_is_class1,
  wolf_is_class2 = wolf_is_class2,
  wolf_first_seen = wolf_first_seen
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
    elk_s_c = rep(0.22, n_years),
    elk_s_ya = rep(0.90, n_years),
    elk_s_oa = rep(0.80, n_years),
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
    wolf_s_p = rep(0.5, n_years),
    wolf_s_a = rep(0.9, n_years),
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
  
  # regression coefficients
  'beta0_calfSurv', 'beta1_calfSurv',
  'beta0_yaSurv', 'beta1_yaSurv',
  'beta0_oaSurv', 'beta1_oaSurv'
)

################################################################################
###########------------------- Run model ----------------------#################
################################################################################

set.seed(17)
nc <- 3
ni <- 100000
nb <- 20000
th <- 4

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

# SAVE OUTPUT
#stop('The following line will overwrite data. Are you sure you would like to proceed?')
#save.image('data/outputs/ICM_environment_2026-03-30.RData')

load('data/outputs/ICM_environment_2026-03-31.RData')

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
round(MCMCsummary(icm_clean, params = "all"), 2)

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
  params = c("beta0", 'beta1')
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
##########---------- Regression plot + coefficient densities ---------##########
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

# yearly elk calf survival summaries
elk_pts <- elk_vrates2 %>%
  filter(rate == "Calf survival (s_c)") %>%
  transmute(
    year,
    elk_s_c = mean,
    elk_low = low,
    elk_high = high
  )

# combine into one plotting dataframe
plot_df <- left_join(wolf_pts, elk_pts, by = "year")

# posterior draws from cleaned chains
post_mat <- do.call(rbind, lapply(icm_clean, as.matrix))

# x grid on raw wolf abundance scale
x_grid <- seq(
  min(plot_df$wolf_N_tot, na.rm = TRUE),
  max(plot_df$wolf_N_tot, na.rm = TRUE),
  length.out = 200
)

# standardize internally because model used standardized wolf abundance
x_grid_std <- (x_grid - icm_constants$wolf_tot_mean) / icm_constants$wolf_tot_sd

# fitted values for every posterior draw across x grid
pred_mat <- sapply(x_grid_std, function(x) {
  plogis(post_mat[, "beta0"] + post_mat[, "beta1"] * x)
})

# summarize into mean and 95% credible ribbon
line_df <- data.frame(
  wolf_N_tot = x_grid,
  elk_s_c = apply(pred_mat, 2, mean),
  elk_low = apply(pred_mat, 2, quantile, probs = 0.025),
  elk_high = apply(pred_mat, 2, quantile, probs = 0.975)
)

# main plot with ribbon
main_plot <- ggplot(plot_df, aes(x = wolf_N_tot, y = elk_s_c)) +
  geom_ribbon(
    data = line_df,
    aes(x = wolf_N_tot, ymin = elk_low, ymax = elk_high),
    inherit.aes = FALSE,
    fill = "#6F263D",
    alpha = 0.25
  ) +
  geom_line(
    data = line_df,
    aes(x = wolf_N_tot, y = elk_s_c),
    inherit.aes = FALSE,
    color = "#6F263D",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = elk_low, ymax = elk_high), width = 0) +
  geom_errorbarh(aes(xmin = wolf_low, xmax = wolf_high), height = 0) +
  geom_point(shape = 21, fill = "gold", color = "black", size = 2.5) +
  theme_classic() +
  labs(
    x = "Wolf abundance",
    y = "Elk calf survival",
    title = "Elk calf survival vs. wolf abundance"
  )

# posterior densities of coefficients
beta_df <- data.frame(
  beta0 = post_mat[, "beta0"],
  beta1 = post_mat[, "beta1"]
)

beta0_plot <- ggplot(beta_df, aes(x = beta0)) +
  geom_density(fill = "#236192", alpha = 0.45) +
  theme_classic() +
  labs(
    x = "beta0",
    y = "Density",
    title = "Posterior of intercept"
  )

beta1_plot <- ggplot(beta_df, aes(x = beta1)) +
  geom_density(fill = "#236192", alpha = 0.45) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic() +
  labs(
    x = "beta1",
    y = "Density",
    title = "Posterior of wolf effect"
  )

# combine bottom row
bottom_row <- plot_grid(
  beta0_plot,
  beta1_plot,
  ncol = 2,
  labels = c("B", "C")
)

# final combined figure
final_plot <- plot_grid(
  main_plot,
  bottom_row,
  ncol = 1,
  rel_heights = c(2, 1),
  labels = c("A", "")
)

final_plot