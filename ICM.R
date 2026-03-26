### Integrated Community Model (ICM)
### Combines elk IPM + wolf IPM in one NIMBLE model
### Last updated: Mar. 26, 2026
### xprockox@gmail.com

################################################################################
############################## PACKAGES #########################################
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
############################## ELK DATA #########################################
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
############################## WOLF DATA ########################################
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
###################### COMMUNITY ALIGNMENT ######################################
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
############################## ICM CODE #########################################
################################################################################

icm_code <- nimbleCode({
  
  ############################ ELK SUBMODEL ################################
  
  # elk priors
  for (t in 1:n_years) {
    logit(elk_s_c[t]) ~ dnorm(qlogis(0.22), 1 / 0.5^2)
    logit(elk_s_ya[t]) ~ dnorm(qlogis(0.90), 1 / 0.5^2)
    logit(elk_s_oa[t]) ~ dnorm(qlogis(0.80), 1 / 0.5^2)
    logit(elk_p_13[t]) ~ dnorm(qlogis(0.15), 1 / 0.5^2)
  }
  
  for (t in 1:(n_years - 1)) {
    elk_f_ya[t] ~ dbeta(1, 1)
    elk_f_oa[t] ~ dbeta(1, 1)
  }
  
  elk_sigma_obs_female ~ dunif(0.05, 2)
  elk_tau_obs_female <- 1 / (elk_sigma_obs_female^2)
  
  elk_lambda_init_1y ~ dgamma(11.1, 0.00454)
  elk_lambda_init_ya ~ dgamma(11.1, 0.00118)
  elk_lambda_init_oa ~ dgamma(11.1, 0.0111)
  elk_lambda_init_female ~ dgamma(11.1, 0.000863)
  
  # elk initial values
  elk_N_1y[1] ~ dpois(elk_lambda_init_1y)
  elk_N_ya[1] ~ dpois(elk_lambda_init_ya)
  elk_N_oa[1] ~ dpois(elk_lambda_init_oa)
  
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
      
      elk_z[i, t] ~ dbern(
        equals(t, elk_first_seen[i]) +
          step(t - elk_first_seen[i] - 0.5) *
          (1 - equals(t, elk_first_seen[i])) *
          elk_z[i, t - 1] * elk_phi[i, t]
      )
    }
    
    for (t in 1:n_years) {
      elk_y[i, t] ~ dbern(elk_p_det[t] * elk_z[i, t])
    }
  }
  
  # elk calf survival 
  for (t in 1:(n_years - 1)) {
    elk_CCR_c_fromYoungCows[t] ~ dbin(elk_f_ya[t] * elk_s_c[t], elk_CCR_cow_youngadult[t])
    elk_CCR_c_fromOldCows[t] ~ dbin(elk_f_oa[t] * elk_s_c[t], elk_CCR_cow_oldadult[t])
    elk_CCR_c[t] <- elk_CCR_c_fromYoungCows[t] + elk_CCR_c_fromOldCows[t]
  }
  
  # elk pregnancy 
  for (t in 1:(n_years - 1)) {
    elk_young_num_preg[t] ~ dbin(elk_f_ya[t], elk_young_num_capt[t])
    elk_old_num_preg[t] ~ dbin(elk_f_oa[t], elk_old_num_capt[t])
  }
  
  # elk growth (p_13)
  for (t in 1:n_years) {
    elk_harvested_13yo[t] ~ dbinom(elk_p_13[t], elk_harvested_ya[t])
  }
  
  ############################ WOLF SUBMODEL ###############################

  # wolf priors
  for (t in 1:n_years) {
    logit(wolf_s_p[t]) ~ dnorm(qlogis(0.5), 1 / 0.5^2)
    logit(wolf_s_a[t]) ~ dnorm(qlogis(0.9), 1 / 0.5^2)
  }
  
  for (t in 1:(n_years - 1)) {
    wolf_f[t] ~ dbeta(1, 1)
  }
  
  wolf_sigma_obs ~ dunif(0.05, 2)
  wolf_tau_obs <- 1 / (wolf_sigma_obs^2)
  
  wolf_lambda_init_a ~ dgamma(10, 1)
  wolf_N_a[1] ~ dpois(wolf_lambda_init_a)
  
  for (t in 1:n_years) {
    wolf_N_p[t] ~ dbin(wolf_s_p[t], wolf_obs_p_sum[t])
    wolf_N_tot[t] <- wolf_N_p[t] + wolf_N_a[t]
    wolf_obs_tot[t] ~ dlnorm(log(wolf_N_tot[t] + 1e-6), wolf_tau_obs)
  }
  
  # wolf state-space 
  for (t in 1:(n_years - 1)) {
    
    wolf_mu_a[t + 1] <- wolf_N_p[t] + wolf_s_a[t] * wolf_N_a[t]
    wolf_N_a[t + 1] ~ dpois(max(1e-6, wolf_mu_a[t + 1]))
    
    wolf_lambda_p_sum[t + 1] <- wolf_f[t] * wolf_N_a[t]
    wolf_obs_p_sum[t + 1] ~ dpois(max(1e-6, wolf_lambda_p_sum[t + 1]))
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
      
      wolf_z[i, t] ~ dbern(
        equals(t, wolf_first_seen[i]) +
          step(t - wolf_first_seen[i] - 0.5) *
          (1 - equals(t, wolf_first_seen[i])) *
          wolf_z[i, t - 1] * wolf_phi[i, t]
      )
    }
    
    for (t in 1:n_years) {
      wolf_y[i, t] ~ dbern(wolf_p_det[t] * wolf_z[i, t])
    }
  }
  
  ##########################################################################
  ######################## COMMUNITY HOOKS #################################
  ##########################################################################
  
  # This section is just for derived quantities right now.
  # Once you're ready for true coupling, this is where you can connect the systems.
  
  for (t in 1:n_years) {
    wolf_to_elk_ratio[t] <- wolf_N_tot[t] / (elk_N_female[t] + 1e-6)
  }
})

################################################################################
########################### CONSTANTS AND DATA ##################################
################################################################################

icm_constants <- list(
  n_years = n_years,
  elk_N_indiv = elk_N_indiv,
  wolf_N_indiv = wolf_N_indiv
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
############################## INITIAL VALUES ###################################
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
    elk_lambda_init_1y = max(1, round(elk_init_N1y[1])),
    elk_lambda_init_ya = max(1, round(elk_init_Nya[1])),
    elk_lambda_init_oa = max(1, round(elk_init_Noa[1])),
    elk_lambda_init_female = max(1, round(elk_init_N1y[1] + elk_init_Nya[1] + elk_init_Noa[1])),
    elk_N_1y = pmax(1, elk_init_N1y),
    elk_N_ya = pmax(1, elk_init_Nya),
    elk_N_oa = pmax(1, elk_init_Noa),
    elk_p_det = runif(n_years, 0.6, 0.95),
    elk_z = elk_z_init,
    
    # wolf
    wolf_s_p = rep(0.5, n_years),
    wolf_s_a = rep(0.9, n_years),
    wolf_f = rep(0.5, n_years - 1),
    wolf_sigma_obs = 0.2,
    wolf_lambda_init_a = max(1, wolf_init_Na[1]),
    wolf_N_p = wolf_init_Np,
    wolf_N_a = wolf_init_Na,
    wolf_p_det = runif(n_years, 0.6, 0.95),
    wolf_z = wolf_z_init
  )
}

################################################################################
############################ PARAMETERS #########################################
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
  "wolf_N_p", "wolf_N_a", "wolf_N_tot"
)

################################################################################
############################### RUN MODEL #######################################
################################################################################

set.seed(17)
nc <- 3
ni <- 1000
nb <- 20
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

################################################################################
############################## QUICK CHECK ######################################
################################################################################

round(icm_mod$summary$all.chains, 2)