### Integrated Community Model (ICM)
### Combines elk IPM + wolf IPM in one NIMBLE model
### Last updated: May 19, 2026
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
library(parallel)

#
##
###
#### # select whether you want to use the wolf population from the Northern Range:
#### # "NR"
#### # or whether you would like to include interior wolves:
#### # "full"
wolf_range <- "full"
###
##
#

################################################################################
##########------------------- Load elk data -----------------###################
################################################################################

load("data/elk_adultSurvival_cjsMatrices.rData")

# immediately rename to avoid collisions with wolf data
elk_z <- z
elk_y <- y
elk_is_class1 <- is_class1
elk_is_class2 <- is_class2
elk_first_seen <- first_seen

rm(z, y, is_class1, is_class2, first_seen)

colnames(elk_is_class1) <- colnames(elk_z)
colnames(elk_is_class2) <- colnames(elk_z)

elk_dat_n <- read.csv("data/elk_abundanceEstimates_stages.csv")
elk_dat_fec <- read.csv("data/elk_fecundity.csv")

elk_dat_n$n_female <- elk_dat_n$n_cow + (elk_dat_n$n_calf / 2)

elk_shared_years <- intersect(as.numeric(elk_dat_n$year), as.numeric(colnames(elk_z)))

################################################################################
############------------------ Load wolf data -----------------#################
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

# drop the first three years (optional)
wolf_pop <- wolf_pop[wolf_pop$seasonal.year %in% c(1998:2025),]

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
elk_dat_n <- elk_dat_n %>%
  filter(year %in% community_years) %>%
  arrange(match(year, community_years))

elk_dat_fec <- elk_dat_fec %>%
  filter(year %in% community_years) %>%
  arrange(match(year, community_years))

elk_z <- elk_z[, match(community_years, colnames(elk_z)), drop = FALSE]
elk_y <- elk_y[, match(community_years, colnames(elk_y)), drop = FALSE]
elk_is_class1 <- elk_is_class1[, match(community_years, colnames(elk_is_class1)), drop = FALSE]
elk_is_class2 <- elk_is_class2[, match(community_years, colnames(elk_is_class2)), drop = FALSE]

elk_first_seen <- apply(elk_y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) ncol(elk_y) else first
})
elk_first_seen <- as.integer(elk_first_seen)

# trim wolf to community years
wolf_pop <- wolf_pop %>%
  filter(seasonal.year %in% community_years) %>%
  arrange(match(seasonal.year, community_years))

wolf_z <- wolf_z[, match(community_years, colnames(wolf_z)), drop = FALSE]
wolf_y <- wolf_y[, match(community_years, colnames(wolf_y)), drop = FALSE]
wolf_is_class1 <- wolf_is_class1[, match(community_years, colnames(wolf_is_class1)), drop = FALSE]
wolf_is_class2 <- wolf_is_class2[, match(community_years, colnames(wolf_is_class2)), drop = FALSE]

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
##########---------------- M-array helper function -----------------############
################################################################################

# build an M-array for animals first released in a given class
# class_mat should be a 0/1 matrix same dimension as y
# it assigns each individual to a class at each occasion
make_class_marray <- function(y_mat, class_mat, class_name = "class") {
  y_mat <- as.matrix(y_mat)
  class_mat <- as.matrix(class_mat)
  
  storage.mode(y_mat) <- "numeric"
  storage.mode(class_mat) <- "numeric"
  
  if (!all(dim(y_mat) == dim(class_mat))) {
    stop("y_mat and class_mat must have the same dimensions.")
  }
  
  n_ind <- nrow(y_mat)
  n_occasions <- ncol(y_mat)
  
  first_seen_vec <- apply(y_mat, 1, function(x) {
    dets <- which(x == 1)
    if (length(dets) == 0) NA_integer_ else min(dets)
  })
  
  # keep only individuals first detected while in this class
  keep <- rep(FALSE, n_ind)
  
  for (i in seq_len(n_ind)) {
    f <- first_seen_vec[i]
    if (is.na(f) || f >= n_occasions) next
    keep[i] <- class_mat[i, f] == 1
  }
  
  y_sub <- y_mat[keep, , drop = FALSE]
  first_seen_sub <- first_seen_vec[keep]
  
  m_array <- matrix(
    0,
    nrow = n_occasions - 1,
    ncol = n_occasions,
    dimnames = list(
      paste0(class_name, "_rel_", 1:(n_occasions - 1)),
      c(paste0("t", 2:n_occasions), "never")
    )
  )
  
  for (i in seq_len(nrow(y_sub))) {
    f <- first_seen_sub[i]
    
    if (is.na(f) || f >= n_occasions) next
    
    later <- which(y_sub[i, (f + 1):n_occasions] == 1)
    
    if (length(later) == 0) {
      m_array[f, n_occasions] <- m_array[f, n_occasions] + 1
    } else {
      first_recap <- f + later[1]
      m_array[f, first_recap - 1] <- m_array[f, first_recap - 1] + 1
    }
  }
  
  releases <- rowSums(m_array)
  
  list(
    m_array = m_array,
    releases = releases,
    first_seen = first_seen_sub,
    keep = keep,
    n_individuals = sum(keep)
  )
}

################################################################################
##########------------- Build class-specific elk M-arrays ----------------######
################################################################################

# young adults = elk_is_class1
elk_m_ya <- make_class_marray(
  y_mat = elk_y,
  class_mat = elk_is_class1,
  class_name = "young_adult"
)

# old adults = elk_is_class2
elk_m_oa <- make_class_marray(
  y_mat = elk_y,
  class_mat = elk_is_class2,
  class_name = "old_adult"
)

# then write parts of the M-array dfs to specific objects
elk_marray_ya <- elk_m_ya$m_array
elk_rel_ya <- elk_m_ya$releases

elk_marray_oa <- elk_m_oa$m_array
elk_rel_oa <- elk_m_oa$releases

################################################################################
##########------------ Build class-specific wolf M-arrays -----------###########
################################################################################

# pups = wolf_is_class1
wolf_m_p <- make_class_marray(
  y_mat = wolf_y,
  class_mat = wolf_is_class1,
  class_name = "pup"
)

# adults = wolf_is_class2
wolf_m_a <- make_class_marray(
  y_mat = wolf_y,
  class_mat = wolf_is_class2,
  class_name = "adult"
)

# write parts of the M-array objects to specific inputs
wolf_marray_p <- wolf_m_p$m_array
wolf_rel_p <- wolf_m_p$releases
wolf_marray_a <- wolf_m_a$m_array
wolf_rel_a <- wolf_m_a$releases

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
covars <- covars %>%
  filter(year %in% community_years) %>%
  arrange(match(year, community_years))

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
  
  # elk M-array detection model
  for (t in 1:n_years) {
    elk_p_det[t] ~ dunif(0, 1)
  }
  
  # young-adult release cohorts
  for (r in 1:(n_years - 1)) {
    
    for (j in 1:(n_years - 1)) {
      elk_marr_prob_ya[r, j] <-
        equals(j, r) * elk_s_ya[r] * elk_p_det[r + 1] +
        step(j - r - 0.5) *
        prod(elk_s_ya[r:j]) *
        prod(1 - elk_p_det[(r + 1):j]) *
        elk_p_det[j + 1]
    }
    
    elk_marr_prob_ya[r, n_years] <- 1 - sum(elk_marr_prob_ya[r, 1:(n_years - 1)])
    elk_marray_ya[r, 1:n_years] ~ dmulti(elk_marr_prob_ya[r, 1:n_years], elk_rel_ya[r])
  }
  
  # old-adult release cohorts
  for (r in 1:(n_years - 1)) {
    
    for (j in 1:(n_years - 1)) {
      elk_marr_prob_oa[r, j] <-
        equals(j, r) * elk_s_oa[r] * elk_p_det[r + 1] +
        step(j - r - 0.5) *
        prod(elk_s_oa[r:j]) *
        prod(1 - elk_p_det[(r + 1):j]) *
        elk_p_det[j + 1]
    }
    
    elk_marr_prob_oa[r, n_years] <- 1 - sum(elk_marr_prob_oa[r, 1:(n_years - 1)])
    elk_marray_oa[r, 1:n_years] ~ dmulti(elk_marr_prob_oa[r, 1:n_years], elk_rel_oa[r])
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
  
  # wolf M-array detection model
  for (t in 1:n_years) {
    wolf_p_det[t] ~ dunif(0, 1)
  }
  
  # pup release cohorts
  for (r in 1:(n_years - 1)) {
    
    for (j in 1:(n_years - 1)) {
      wolf_marr_prob_p[r, j] <-
        equals(j, r) * wolf_s_p[r] * wolf_p_det[r + 1] +
        step(j - r - 0.5) *
        wolf_s_p[r] *
        prod(wolf_s_a[(r + 1):j]) *
        prod(1 - wolf_p_det[(r + 1):j]) *
        wolf_p_det[j + 1]
    }
    
    wolf_marr_prob_p[r, n_years] <- 1 - sum(wolf_marr_prob_p[r, 1:(n_years - 1)])
    wolf_marray_p[r, 1:n_years] ~ dmulti(wolf_marr_prob_p[r, 1:n_years], wolf_rel_p[r])
  }
  
  # adult release cohorts
  for (r in 1:(n_years - 1)) {
    
    for (j in 1:(n_years - 1)) {
      wolf_marr_prob_a[r, j] <-
        equals(j, r) * wolf_s_a[r] * wolf_p_det[r + 1] +
        step(j - r - 0.5) *
        prod(wolf_s_a[r:j]) *
        prod(1 - wolf_p_det[(r + 1):j]) *
        wolf_p_det[j + 1]
    }
    
    wolf_marr_prob_a[r, n_years] <- 1 - sum(wolf_marr_prob_a[r, 1:(n_years - 1)])
    wolf_marray_a[r, 1:n_years] ~ dmulti(wolf_marr_prob_a[r, 1:n_years], wolf_rel_a[r])
  }
  
  ##########---------------- REGRESSIONS ----------------#############
  
  # Priors for elk regression coefficients
  beta0_calfSurv ~ dnorm(qlogis(0.22), 1 / 0.3^2) # mean = 0.22
  beta1_calfSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_calfSurv_wintPPT ~ dnorm(0, 1 / 0.3^2)
  beta3_calfSurv_grizN ~ dnorm(0, 1 / 0.3^2)
  beta4_calfSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  
  beta0_yaSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2) # mean = 0.90
  beta1_yaSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_yaSurv_wintPPT ~ dnorm(0, 1 / 0.3^2)
  beta3_yaSurv_grizN ~ dnorm(0, 1 / 0.3^2)
  beta4_yaSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  
  beta0_oaSurv ~ dnorm(qlogis(0.80), 1 / 0.3^2) # mean = 0.80
  beta1_oaSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_oaSurv_wintPPT ~ dnorm(0, 1 / 0.3^2)
  beta3_oaSurv_grizN ~ dnorm(0, 1 / 0.3^2)
  beta4_oaSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  
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
  beta3_wpupSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  
  beta0_wadSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2) # mean = 0.90
  beta1_wadSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  beta2_wadSurv_bisonN ~ dnorm(0, 1 / 0.3^2)
  beta3_wadSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  
  # Priors for wolf random year-effect SDs
  sigma_wpup ~ dunif(0, 0.5)
  tau_wpup <- 1 / (sigma_wpup^2)
  
  sigma_wad ~ dunif(0, 0.3)
  tau_wad <- 1 / (sigma_wad^2)
  
  # Initial year vital rates must be explicitly defined because
  # lagged predictors (t - 1) are used in the regression models.
  logit(elk_s_c[1]) ~ dnorm(qlogis(0.22), 1 / 0.5^2)
  logit(elk_s_ya[1]) ~ dnorm(qlogis(0.90), 1 / 0.5^2)
  logit(elk_s_oa[1]) ~ dnorm(qlogis(0.80), 1 / 0.5^2)
  
  logit(wolf_s_p[1]) ~ dnorm(qlogis(0.50), 1 / 0.5^2)
  logit(wolf_s_a[1]) ~ dnorm(qlogis(0.90), 1 / 0.5^2)
  
  # standardize abundance data for use in regressions
  for (t in 1:n_years) {
    wolf_N_tot_std[t] <- (wolf_N_tot[t] - wolf_tot_mean) / wolf_tot_sd
    elk_N_female_std[t] <- (elk_N_female[t] - elk_N_female_mean) / elk_N_female_sd
  }
  
  # first-year random effects fixed to zero (otherwise be NA and lead to issues extracting param ests.)
  eps_elk_s_c[1] <- 0
  eps_elk_s_ya[1] <- 0
  eps_elk_s_oa[1] <- 0
  eps_wolf_s_p[1] <- 0
  eps_wolf_s_a[1] <- 0
  
  # Year-specific regressions
  for (t in 2:n_years) {
    
    # elk random year effects
    eps_elk_s_c[t] ~ dnorm(0, tau_calf)
    eps_elk_s_ya[t] ~ dnorm(0, tau_ya)
    eps_elk_s_oa[t] ~ dnorm(0, tau_oa)
    
    # wolf random year effects
    eps_wolf_s_p[t] ~ dnorm(0, tau_wpup)
    eps_wolf_s_a[t] ~ dnorm(0, tau_wad)
    
    # elk regression models
    logit(elk_s_c[t])  <- 
      beta0_calfSurv + 
      beta1_calfSurv_wolfN * wolf_N_tot_std[t - 1] + 
      beta2_calfSurv_wintPPT * wintPPT[t] + 
      beta3_calfSurv_grizN * grizN_std[t - 1] +
      beta4_calfSurv_elkN * elk_N_female_std[t - 1] +
      eps_elk_s_c[t]
    
    logit(elk_s_ya[t]) <- 
      beta0_yaSurv + 
      beta1_yaSurv_wolfN * wolf_N_tot_std[t - 1] + 
      beta2_yaSurv_wintPPT * wintPPT[t] + 
      beta3_yaSurv_grizN * grizN_std[t - 1] +
      beta4_yaSurv_elkN * elk_N_female_std[t - 1] +
      eps_elk_s_ya[t]
    
    logit(elk_s_oa[t]) <- 
      beta0_oaSurv + 
      beta1_oaSurv_wolfN * wolf_N_tot_std[t - 1] + 
      beta2_oaSurv_wintPPT * wintPPT[t] + 
      beta3_oaSurv_grizN * grizN_std[t - 1] +
      beta4_oaSurv_elkN * elk_N_female_std[t - 1] +
      eps_elk_s_oa[t]
    
    # wolf regression models
    logit(wolf_s_p[t])  <- 
      beta0_wpupSurv + 
      beta1_wpupSurv_elkN * elk_N_female_std[t - 1] + 
      beta2_wpupSurv_bisonN * bisonN_std[t - 1] +
      beta3_wpupSurv_wolfN * wolf_N_tot_std[t - 1] +
      eps_wolf_s_p[t]
    
    logit(wolf_s_a[t]) <- 
      beta0_wadSurv + 
      beta1_wadSurv_elkN * elk_N_female_std[t - 1] + 
      beta2_wadSurv_bisonN * bisonN_std[t - 1] +
      beta3_wadSurv_wolfN * wolf_N_tot_std[t - 1] +
      eps_wolf_s_a[t]
  }
})

################################################################################
###########---------------- Constants and data ----------------#################
################################################################################

icm_constants <- list(
  n_years = n_years,
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
  elk_marray_ya = elk_marray_ya,
  elk_rel_ya = elk_rel_ya,
  elk_marray_oa = elk_marray_oa,
  elk_rel_oa = elk_rel_oa,
  
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
  wolf_marray_p = wolf_marray_p,
  wolf_rel_p = wolf_rel_p,
  wolf_marray_a = wolf_marray_a,
  wolf_rel_a = wolf_rel_a,
  
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
  
  wolf_init_Np_bio <- rep(0, n_years)
  
  if (n_years >= 2) {
    wolf_init_Np_bio[2] <- max(0, wolf_init_Np[2] - 9)
  }
  
  if (n_years >= 3) {
    wolf_init_Np_bio[3:n_years] <- wolf_init_Np[3:n_years]
  }
  
  wolf_init_Np_bio <- pmin(wolf_init_Np_bio, pmax(0, round(wolf_pop$summer_pups)))
  wolf_init_Np_bio[1] <- 0
  
  wolf_init_Na <- pmax(1, round(wolf_init_Ntot - wolf_init_Np))
  
  list(
    # elk
    elk_p_13 = rep(0.15, n_years),
    elk_f_ya = rep(0.76, n_years - 1),
    elk_f_oa = rep(0.64, n_years - 1),
    elk_sigma_obs_female = 0.30,
    elk_N_1y = pmax(1, elk_init_N1y),
    elk_N_ya = pmax(1, elk_init_Nya),
    elk_N_oa = pmax(1, elk_init_Noa),
    elk_p_det = runif(n_years, 0.6, 0.95),

    # wolf
    wolf_f = rep(1.0, n_years),
    wolf_sigma_obs = 0.2,
    wolf_N_p_sum = pmax(1, round(wolf_pop$summer_pups)),
    wolf_N_p = wolf_init_Np,
    wolf_N_p_bio = wolf_init_Np_bio,
    wolf_N_a = wolf_init_Na,
    wolf_p_det = runif(n_years, 0.6, 0.95),

    # elk regression coefficients
    beta0_calfSurv = qlogis(0.22),
    beta1_calfSurv_wolfN = 0,
    beta2_calfSurv_wintPPT = 0,
    beta3_calfSurv_grizN = 0,
    beta4_calfSurv_elkN = 0,
    
    beta0_yaSurv = qlogis(0.90),
    beta1_yaSurv_wolfN = 0,
    beta2_yaSurv_wintPPT = 0,
    beta3_yaSurv_grizN = 0,
    beta4_yaSurv_elkN = 0,
    
    beta0_oaSurv = qlogis(0.80),
    beta1_oaSurv_wolfN = 0,
    beta2_oaSurv_wintPPT = 0,
    beta3_oaSurv_grizN = 0,
    beta4_oaSurv_elkN = 0,
    
    sigma_calf = 0.1,
    sigma_ya = 0.1,
    sigma_oa = 0.1,
    
    eps_elk_s_c = c(0, rep(0, n_years - 1)),
    eps_elk_s_ya = c(0, rep(0, n_years - 1)),
    eps_elk_s_oa = c(0, rep(0, n_years - 1)),
    
    # wolf regression coefficients
    beta0_wpupSurv = qlogis(0.50),
    beta1_wpupSurv_elkN = 0,
    beta2_wpupSurv_bisonN = 0,
    beta3_wpupSurv_wolfN = 0,
    
    beta0_wadSurv = qlogis(0.90),
    beta1_wadSurv_elkN = 0,
    beta2_wadSurv_bisonN = 0,
    beta3_wadSurv_wolfN = 0,
    
    sigma_wpup = 0.1,
    sigma_wad = 0.1,
    
    eps_wolf_s_p = c(0, rep(0, n_years - 1)),
    eps_wolf_s_a = c(0, rep(0, n_years - 1)),
    
    # first-year survival values
    logit_wolf_s_p = c(qlogis(0.5), rep(NA, n_years - 1)),
    logit_wolf_s_a = c(qlogis(0.9), rep(NA, n_years - 1)),
    logit_elk_s_c = c(qlogis(0.22), rep(NA, n_years - 1)),
    logit_elk_s_ya = c(qlogis(0.90), rep(NA, n_years - 1)),
    logit_elk_s_oa = c(qlogis(0.80), rep(NA, n_years - 1))
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
  'beta0_calfSurv', 'beta1_calfSurv_wolfN', 'beta2_calfSurv_wintPPT', 'beta3_calfSurv_grizN', 'beta4_calfSurv_elkN',
  'beta0_yaSurv', 'beta1_yaSurv_wolfN', 'beta2_yaSurv_wintPPT', 'beta3_yaSurv_grizN', 'beta4_yaSurv_elkN',
  'beta0_oaSurv', 'beta1_oaSurv_wolfN', 'beta2_oaSurv_wintPPT', 'beta3_oaSurv_grizN', 'beta4_oaSurv_elkN',
  
  "sigma_calf", "sigma_ya", "sigma_oa",
  "eps_elk_s_c", "eps_elk_s_ya", "eps_elk_s_oa",
  
  # wolf regression coefficients
  'beta0_wpupSurv', 'beta1_wpupSurv_elkN', 'beta2_wpupSurv_bisonN', 'beta3_wpupSurv_wolfN',
  'beta0_wadSurv', 'beta1_wadSurv_elkN', 'beta2_wadSurv_bisonN', 'beta3_wadSurv_wolfN',
  
  "sigma_wpup", "sigma_wad",
  "eps_wolf_s_p", "eps_wolf_s_a"
)

################################################################################
###########------------- Run model in parallel by chain -------------###########
################################################################################

library(parallel)

set.seed(17)
nc <- 3
ni <- 200000
nb <- 40000
th <- 4

# precompute objects used by make_icm_inits() so workers do not need
# to look back into your global environment
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

wolf_summer_pups <- wolf_pop$summer_pups

wolf_init_Ntot <- ifelse(
  is.na(wolf_pop$total_abundance),
  max(1, round(mean(wolf_pop$total_abundance, na.rm = TRUE))),
  pmax(1, round(wolf_pop$total_abundance))
)

wolf_init_Np <- ifelse(
  is.na(wolf_summer_pups),
  max(1, round(mean(wolf_summer_pups, na.rm = TRUE) * 0.5)),
  pmax(1, round(wolf_summer_pups * 0.5))
)

wolf_init_Np_bio <- rep(0, n_years)

if (n_years >= 2) {
  wolf_init_Np_bio[2] <- max(0, wolf_init_Np[2] - 9)
}

if (n_years >= 3) {
  wolf_init_Np_bio[3:n_years] <- wolf_init_Np[3:n_years]
}

wolf_init_Np_bio <- pmin(wolf_init_Np_bio, pmax(0, round(wolf_summer_pups)))
wolf_init_Np_bio[1] <- 0

wolf_init_Na <- pmax(1, round(wolf_init_Ntot - wolf_init_Np))

make_icm_inits <- function() {
  list(
    # elk
    elk_p_13 = rep(0.15, n_years),
    elk_f_ya = rep(0.76, n_years - 1),
    elk_f_oa = rep(0.64, n_years - 1),
    elk_sigma_obs_female = 0.30,
    elk_N_1y = pmax(1, elk_init_N1y),
    elk_N_ya = pmax(1, elk_init_Nya),
    elk_N_oa = pmax(1, elk_init_Noa),
    elk_p_det = runif(n_years, 0.6, 0.95),
    
    # wolf
    wolf_f = rep(1.0, n_years),
    wolf_sigma_obs = 0.2,
    wolf_N_p_sum = pmax(1, round(wolf_summer_pups)),
    wolf_N_p = wolf_init_Np,
    wolf_N_p_bio = wolf_init_Np_bio,
    wolf_N_a = wolf_init_Na,
    wolf_p_det = runif(n_years, 0.6, 0.95),
    
    # elk regression coefficients
    beta0_calfSurv = qlogis(0.22),
    beta1_calfSurv_wolfN = 0,
    beta2_calfSurv_wintPPT = 0,
    beta3_calfSurv_grizN = 0,
    beta4_calfSurv_elkN = 0,
    
    beta0_yaSurv = qlogis(0.90),
    beta1_yaSurv_wolfN = 0,
    beta2_yaSurv_wintPPT = 0,
    beta3_yaSurv_grizN = 0,
    beta4_yaSurv_elkN = 0,
    
    beta0_oaSurv = qlogis(0.80),
    beta1_oaSurv_wolfN = 0,
    beta2_oaSurv_wintPPT = 0,
    beta3_oaSurv_grizN = 0,
    beta4_oaSurv_elkN = 0,
    
    sigma_calf = 0.1,
    sigma_ya = 0.1,
    sigma_oa = 0.1,
    
    eps_elk_s_c = c(0, rep(0, n_years - 1)),
    eps_elk_s_ya = c(0, rep(0, n_years - 1)),
    eps_elk_s_oa = c(0, rep(0, n_years - 1)),
    
    # wolf regression coefficients
    beta0_wpupSurv = qlogis(0.50),
    beta1_wpupSurv_elkN = 0,
    beta2_wpupSurv_bisonN = 0,
    beta3_wpupSurv_wolfN = 0,
    
    beta0_wadSurv = qlogis(0.90),
    beta1_wadSurv_elkN = 0,
    beta2_wadSurv_bisonN = 0,
    beta3_wadSurv_wolfN = 0,
    
    sigma_wpup = 0.1,
    sigma_wad = 0.1,
    
    eps_wolf_s_p = c(0, rep(0, n_years - 1)),
    eps_wolf_s_a = c(0, rep(0, n_years - 1)),
    
    # first-year survival values
    logit_wolf_s_p = c(qlogis(0.5), rep(NA, n_years - 1)),
    logit_wolf_s_a = c(qlogis(0.9), rep(NA, n_years - 1)),
    logit_elk_s_c = c(qlogis(0.22), rep(NA, n_years - 1)),
    logit_elk_s_ya = c(qlogis(0.90), rep(NA, n_years - 1)),
    logit_elk_s_oa = c(qlogis(0.80), rep(NA, n_years - 1))
  )
}

run_one_chain <- function(chain_id,
                          icm_code,
                          icm_data,
                          icm_constants,
                          icm_params,
                          make_icm_inits,
                          ni,
                          nb,
                          th) {
  
  library(nimble)
  library(coda)
  
  set.seed(1000 + chain_id)
  
  fit <- nimbleMCMC(
    code = icm_code,
    data = icm_data,
    constants = icm_constants,
    inits = make_icm_inits(),
    monitors = icm_params,
    niter = ni,
    nburnin = nb,
    thin = th,
    nchains = 1,
    samplesAsCodaMCMC = TRUE,
    summary = FALSE
  )
  
  fit
}

start_time <- Sys.time()

cl <- makeCluster(nc)

clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  NULL
})

clusterExport(
  cl,
  varlist = c(
    "icm_code",
    "icm_data",
    "icm_constants",
    "icm_params",
    "make_icm_inits",
    "run_one_chain",
    "ni",
    "nb",
    "th",
    "n_years",
    "elk_init_N1y",
    "elk_init_Nya",
    "elk_init_Noa",
    "wolf_summer_pups",
    "wolf_init_Np",
    "wolf_init_Np_bio",
    "wolf_init_Na"
  ),
  envir = environment()
)

# optional test: make sure workers can build inits
clusterEvalQ(cl, {
  tmp <- make_icm_inits()
  names(tmp)[1:5]
})

chain_samples <- parLapply(
  cl,
  X = 1:nc,
  fun = function(chain_id) {
    run_one_chain(
      chain_id = chain_id,
      icm_code = icm_code,
      icm_data = icm_data,
      icm_constants = icm_constants,
      icm_params = icm_params,
      make_icm_inits = make_icm_inits,
      ni = ni,
      nb = nb,
      th = th
    )
  }
)

stopCluster(cl)

# combine chains into one mcmc.list
icm_samples <- mcmc.list(chain_samples)

# original combined matrix used only to find bad columns
post_mat <- as.matrix(icm_samples)

bad_cols <- which(
  apply(post_mat, 2, function(x) any(!is.finite(x)))
)

good_cols <- colnames(post_mat)[-bad_cols]

chain_samples_clean <- lapply(chain_samples, function(ch) {
  ch_mat <- as.matrix(ch)
  mcmc(ch_mat[, good_cols, drop = FALSE])
})

icm_samples_clean <- mcmc.list(chain_samples_clean)

icm_summary <- MCMCsummary(icm_samples_clean)

post_mat_clean <- as.matrix(icm_samples_clean)

end_time <- Sys.time()
run_time <- end_time - start_time

print(paste0('Model runtime: ', 
             round(run_time, 2), 
             ' ', 
             units(run_time)))

# save outputs
# stop("The following line will overwrite data. Are you sure you would like to proceed?")
save(
  icm_samples,
  icm_summary,
  file = "data/outputs/ICM_parallel_output_2026-05-19.RData"
)

################################################################################
################################################################################
################################################################################