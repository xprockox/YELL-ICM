### Integrated Community Model (ICM)
### Combines elk IPM + wolf IPM in one NIMBLE model
### Last updated: July 23, 2026

################################################################################
############################ Packages and settings #############################
################################################################################

library(tidyverse)
library(MCMCvis)
library(nimble)
library(coda)
library(cowplot)
library(popbio)
library(parallel)

wolf_range <- "full" # can be "NR" or "full"
grizzly_range <- "elk_mcp" # can be "NR", "park", "elk_mcp", or "full"
drop_regression_years <- 1995:2001
# drop intial years 95-01 because 
# (a) years 1995-1997 were heavily influenced by humans (esp. wolves) 
# (b) MODIS data only goes back to 2001
# (c) elk collar data only goes back to 2001
# (c) MODIS data is lagged in regressions, so the 2002 regressions rely on 2001 data

################################################################################
################################## Helpers #####################################
################################################################################

# reconstruct matrix identifying when a collared individual was first seen
rebuild_first_seen <- function(y_mat) {
  apply(y_mat, 1, function(row) {
    first <- which(row == 1)[1]
    if (is.na(first)) ncol(y_mat) else first
  }) %>%
    as.integer()
}

# trims matrices to shared years
trim_capture_object <- function(z_mat, y_mat, class1_mat, class2_mat, years) {
  idx <- match(years, colnames(z_mat))
  
  list(
    z = z_mat[, idx, drop = FALSE],
    y = y_mat[, idx, drop = FALSE],
    is_class1 = class1_mat[, idx, drop = FALSE],
    is_class2 = class2_mat[, idx, drop = FALSE]
  )
}

# convert capture history matrices to M-array format
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
  
  list(
    m_array = m_array,
    releases = rowSums(m_array),
    first_seen = first_seen_sub,
    keep = keep,
    n_individuals = sum(keep)
  )
}

# extract M-arrays
extract_marray_inputs <- function(y_mat, class_mat, class_name) {
  out <- make_class_marray(y_mat, class_mat, class_name)
  list(marray = out$m_array, releases = out$releases)
}

################################################################################
################################ Load elk data #################################
################################################################################

load("data/elk_adultSurvival_cjsMatrices.rData")
elk_z <- z
elk_y <- y
elk_is_class1 <- is_class1
elk_is_class2 <- is_class2
rm(z, y, is_class1, is_class2, first_seen)

colnames(elk_is_class1) <- colnames(elk_z)
colnames(elk_is_class2) <- colnames(elk_z)

elk_dat_n <- read.csv("data/elk_abundanceEstimates_stages.csv")
elk_dat_fec <- read.csv("data/elk_fecundity.csv")

elk_dat_n$n_female <- elk_dat_n$n_cow + (elk_dat_n$n_calf / 2)
elk_shared_years <- intersect(elk_dat_n$year, as.numeric(colnames(elk_z)))

################################################################################
################################ Load wolf data ################################
################################################################################

load("data/wolf_adultSurvival_cjsMatrices.rData")
wolf_z <- z
wolf_y <- y
wolf_is_class1 <- is_class1
wolf_is_class2 <- is_class2
rm(z, y, is_class1, is_class2, first_seen)

colnames(wolf_is_class1) <- colnames(wolf_z)
colnames(wolf_is_class2) <- colnames(wolf_z)

wolf_pop <- switch(
  wolf_range,
  NR = read.csv("data/wolf_nr_pop.csv"),
  full = read.csv("data/wolf_full_park_pop.csv"),
  stop("wolf_range must be 'NR' or 'full'")
)

wolf_shared_years <- intersect(wolf_pop$seasonal.year, as.numeric(colnames(wolf_z)))
wolf_pop <- wolf_pop %>%
  filter(seasonal.year %in% wolf_shared_years)

################################################################################
######################## Align elk and wolf to common years ####################
################################################################################

community_years <- intersect(elk_shared_years, wolf_shared_years)

regression_start_year <- max(drop_regression_years) + 1
regression_start_idx <- match(regression_start_year, community_years)

if (is.na(regression_start_idx)) {
  stop("regression_start_year was not found in community_years.")
}

elk_dat_n <- elk_dat_n %>%
  filter(year %in% community_years) %>%
  arrange(match(year, community_years))

elk_dat_fec <- elk_dat_fec %>%
  filter(year %in% community_years) %>%
  arrange(match(year, community_years))

wolf_pop <- wolf_pop %>%
  filter(seasonal.year %in% community_years) %>%
  arrange(match(seasonal.year, community_years))

elk_cap <- trim_capture_object(
  elk_z, elk_y, elk_is_class1, elk_is_class2, community_years
)

wolf_cap <- trim_capture_object(
  wolf_z, wolf_y, wolf_is_class1, wolf_is_class2, community_years
)

elk_z <- elk_cap$z
elk_y <- elk_cap$y
elk_is_class1 <- elk_cap$is_class1
elk_is_class2 <- elk_cap$is_class2
elk_first_seen <- rebuild_first_seen(elk_y)

wolf_z <- wolf_cap$z
wolf_y <- wolf_cap$y
wolf_is_class1 <- wolf_cap$is_class1
wolf_is_class2 <- wolf_cap$is_class2
wolf_first_seen <- rebuild_first_seen(wolf_y)

n_years <- nrow(elk_dat_n)
if (n_years != nrow(wolf_pop)) {
  stop("Elk and wolf year counts do not match after alignment.")
}

################################################################################
############################### Build M-arrays #################################
################################################################################

elk_ya <- extract_marray_inputs(elk_y, elk_is_class1, "young_adult")
elk_oa <- extract_marray_inputs(elk_y, elk_is_class2, "old_adult")

wolf_p <- extract_marray_inputs(wolf_y, wolf_is_class1, "pup")
wolf_a <- extract_marray_inputs(wolf_y, wolf_is_class2, "adult")

################################################################################
############################## Load covariates #################################
################################################################################

# import environmental covar data
env_covars <- read.csv("data/covariates/environmental_covariates_annual.csv")

# import bison obs & cull data
bison <- read.csv("data/covariates/yellowstone_bison_culls_harvests_1970_2023.csv") %>%
  select(Year, Max_Bison_North, Total_Culls_Harvests) %>%
  rename(
    year = Year,
    NR_Bison = Max_Bison_North,
    total_cull_harvest = Total_Culls_Harvests
  ) %>%
  mutate(year = year + 1,
         total_cull_harvest = readr::parse_number(total_cull_harvest),
         NR_Bison = readr::parse_number(NR_Bison))

bison <- tibble(year = community_years) %>%
  left_join(bison, by = "year")

# import grizzly data using conditional logic
if (grizzly_range %in% c("NR", "park", "elk_mcp")) {
  
  grizzly <- read.csv("data/covariates/annual_grizzly_estimates_by_polygon.csv")
  
  polygon_keep <- switch(
    grizzly_range,
    NR = "Northern Range",
    park = "Yellowstone National Park",
    elk_mcp = "Elk GPS MCP"
  )
  
  grizzly <- grizzly %>%
    filter(polygon_name == polygon_keep) %>%
    rename(griz_N = estimated_grizzlies) %>%
    select(year, griz_N)
  
} else if (grizzly_range == "full") {
  
  grizzly <- read.csv("data/covariates/IGBST_abundances_1995-2024.csv") %>%
    rename(
      griz_N = Mean,
      year = Year
    ) %>%
    select(year, griz_N)
  
} else {
  
  stop("grizzly_range must be 'NR', 'park', 'elk_mcp', or 'full'")
  
}

# grizzly data only goes to 2020, so we need to add some NAs at the end to have the model estimate years after
grizzly <- tibble(year = community_years) %>%
  left_join(grizzly, by = "year")

# import annual elk harvest
elk_harvest <- read.csv("data/covariates/annual_elk_harvest.csv") %>%
  rename(
    elk_harvest = total_female_harvested
  )

elk_harvest <- tibble(year = community_years) %>%
  left_join(elk_harvest, by = "year")

# import cougar data
cougars <- read.csv('data/allSpp_Abundances.csv') %>%
  select(Year, Cougars) %>%
  rename(
    year = Year,
    cougar_N = Cougars
  )

cougars <- tibble(year = community_years) %>%
  left_join(cougars, by = "year")


# combine all covars into same df
covars <- env_covars %>%
  left_join(bison, by = "year") %>%
  left_join(grizzly, by = "year") %>%
  left_join(elk_harvest, by = "year") %>%
  left_join(cougars, by = 'year') %>%
  filter(year %in% community_years) %>%
  arrange(match(year, community_years))

# standardize all covars
zscore_na <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

covars_std <- covars %>%
  mutate(across(-c(year, winter_severity), zscore_na)) # winter severity is already standardized

################################################################################
############################### NIMBLE model ###################################
################################################################################

icm_code <- nimbleCode({
  
  ##############################################################
  ########## ------------- ELK SUBMODEL ------------- ##########
  ##############################################################
  
  # prior for probability of a young adult being 13 y.o. 
  for (t in 1:n_years) {
    logit(elk_p_13[t]) ~ dnorm(qlogis(0.15), 1 / 0.5^2)
  }
  
  # prior for fecundity
  for (t in 1:(n_years - 1)) {
    elk_f_ya[t] ~ dbeta(1, 1) # bound 0 - 1; assumes cows have max 1 calf
    elk_f_oa[t] ~ dbeta(1, 1)
  }
  
  # obs. error
  elk_sigma_obs_female ~ dunif(0.05, 2)
  elk_tau_obs_female <- 1 / (elk_sigma_obs_female^2)
  
  # initial values
  elk_N_1y[1] ~ dunif(0, 5000)
  elk_N_ya[1] ~ dunif(0, 15000)
  elk_N_oa[1] ~ dunif(0, 5000)
  
  elk_N_female[1] <- elk_N_1y[1] + elk_N_ya[1] + elk_N_oa[1]
  elk_obs_female[1] ~ dlnorm(log(elk_N_female[1] + 1e-6), elk_tau_obs_female)
  
  # state-space abundance
  for (t in 1:(n_years - 1)) {
    elk_mu_1y[t + 1] <- # yearlings in year t + 1 are comprised of...
      elk_f_ya[t] * elk_s_c[t] * elk_N_ya[t] + # the number of calves born from young adults which survived 
      elk_f_oa[t] * elk_s_c[t] * elk_N_oa[t] # plus the number of calves born from old adults which survived
    
    elk_mu_ya[t + 1] <- # young adults in year t + 1 are comprised of...
      elk_s_ya[t] * elk_N_1y[t] + # yearlings that survived
      elk_s_ya[t] * (1 - elk_p_13[t]) * elk_N_ya[t] #- # plus young adults that survived and did not transition to old adults
      #elk_ya_harvest[t] # minus the number of young adults harvested
    
    elk_mu_oa[t + 1] <- # old adults in year t + 1 are comprised of...
      elk_s_ya[t] * elk_p_13[t] * elk_N_ya[t] + # young adults that survived and transitioned
      elk_s_oa[t] * elk_N_oa[t] #- # plus old adults that survived
      #elk_oa_harvest[t] # minus the number of old adults harvested
    
    # demographic stochasticity (actual numbers differ from expected)
    elk_N_1y[t + 1] ~ dpois(max(1e-6, elk_mu_1y[t + 1]))
    elk_N_ya[t + 1] ~ dpois(max(1e-6, elk_mu_ya[t + 1]))
    elk_N_oa[t + 1] ~ dpois(max(1e-6, elk_mu_oa[t + 1]))
    
    # total number of females accounting for observation error
    elk_N_female[t + 1] <- elk_N_1y[t + 1] + elk_N_ya[t + 1] + elk_N_oa[t + 1]
    elk_obs_female[t + 1] ~ dlnorm(log(elk_N_female[t + 1] + 1e-6), elk_tau_obs_female)
  }
  
  # prior for detection probability in CJS
  for (t in 1:n_years) {
    elk_p_det[t] ~ dunif(0, 1)
  }
  
  # CJS using M-arrays
  for (r in 1:(n_years - 1)) { # release occasions
    for (j in 1:(n_years - 1)) { # recaptures
      
      # ------ YOUNG ADULTS ------
      elk_marr_prob_ya[r, j] <-
        
        # the first recapture occurs immediately next year
        equals(j, r) * # this is 1 when j = r
        elk_s_ya[r] * # elk survived one year
        elk_p_det[r + 1] + # and was detected
        
        # OR the first recapture occurs later than next year
        step(j - r - 0.5) * # this is 1 when j > r
        prod(elk_s_ya[r:j]) * # elk survived from years r:j
        prod(1 - elk_p_det[(r + 1):j]) * # was not detected during that time
        elk_p_det[j + 1] # but was eventually detected
      
      # ------ OLD ADULTS ------
      elk_marr_prob_oa[r, j] <-
        
        # the first recapture occurs immediately next year
        equals(j, r) * 
        elk_s_oa[r] * elk_p_det[r + 1] +
        step(j - r - 0.5) *
        
        # OR the first recapture occurs later than next year
        prod(elk_s_oa[r:j]) *
        prod(1 - elk_p_det[(r + 1):j]) *
        elk_p_det[j + 1]
    }
    
    # this is the probability that an individual is never seen again after release
    elk_marr_prob_ya[r, n_years] <- 1 - sum(elk_marr_prob_ya[r, 1:(n_years - 1)])
    elk_marr_prob_oa[r, n_years] <- 1 - sum(elk_marr_prob_oa[r, 1:(n_years - 1)])
    
    # likelihood is multinomial, because the animal could be seen in 1st recapture, 2nd, 3rd...etc.
    elk_marray_ya[r, 1:n_years] ~ dmulti(elk_marr_prob_ya[r, 1:n_years], elk_rel_ya[r])
    elk_marray_oa[r, 1:n_years] ~ dmulti(elk_marr_prob_oa[r, 1:n_years], elk_rel_oa[r])
  }
  
  # calf survival and fecundity
  for (t in 1:(n_years - 1)) {
    elk_CCR_prob_young[t] <- elk_f_ya[t] * elk_s_c[t]
    elk_CCR_prob_old[t] <- elk_f_oa[t] * elk_s_c[t]
    
    elk_CCR_c_fromYoungCows[t] ~ dbin(elk_CCR_prob_young[t], elk_CCR_cow_youngadult[t])
    elk_CCR_c_fromOldCows[t] ~ dbin(elk_CCR_prob_old[t], elk_CCR_cow_oldadult[t])
    
    elk_young_num_preg[t] ~ dbin(elk_f_ya[t], elk_young_num_capt[t])
    elk_old_num_preg[t] ~ dbin(elk_f_oa[t], elk_old_num_capt[t])
    
    elk_CCR_c[t] <- elk_CCR_c_fromYoungCows[t] + elk_CCR_c_fromOldCows[t]
  }
  
  # proportion of young adults that are 13 y.o.
  for (t in 1:n_years) {
    elk_harvested_13yo[t] ~ dbin(elk_p_13[t], elk_harvested_ya[t])
  }
  
  ###############################################################
  ########## ------------- WOLF SUBMODEL ------------- ##########
  ###############################################################
  
  # prior for fecundity
  for (t in 1:(n_years-1)) {
    wolf_f[t] ~ dgamma(2, 2)
  }
  
  # obs. error
  wolf_sigma_obs ~ dunif(0.05, 2)
  wolf_tau_obs <- 1 / (wolf_sigma_obs^2)
  
  # initial values (using true numbers from reintroduction)
  wolf_N_a[1] <- 8
  wolf_N_p[1] <- 6
  wolf_N_p_sum[1] <- 0 # no summer pups the first year (not reintroduced yet)
  
  # pup survival / introduction
  for (t in 2:n_years) {
    
    # number of pups produced biologically (not introduced) is the number of pups in summer that survived
    wolf_N_p_bio[t] ~ dbin(wolf_s_p[t - 1], wolf_N_p_sum[t])
    
    # this just adds the 9 pups introduced in 1996
    wolf_N_p[t] <- wolf_N_p_bio[t] + equals(t, 2) * 9
  }
  
  # state-space abundance
  for (t in 1:(n_years - 1)) {
    
    wolf_mu_a[t + 1] <- # expected number of adults is comprised of...
      wolf_s_a[t] * wolf_N_a[t] + # the number of adults surviving from the previous year 
      wolf_N_p[t] + # plus the number of pups surviving from the previous year
      equals(t, 1) * 8 # and if t = 1 (year 1996), there's an explicit addition of 8 additional adults introduced
    
    # demographic stochasticity
    wolf_N_a[t + 1] ~ dpois(max(1e-6, wolf_mu_a[t + 1]))
    
    # then the number of pups born in the summer of t + 1 is the number of adults * fecundity
    wolf_N_p_sum[t + 1] ~ dpois(max(1e-6, wolf_f[t] * wolf_N_a[t]))
  }
  
  # observation processes
  for (t in 1:n_years) {
    
    # deterministic calculation of total number of wolves
    wolf_N_tot[t] <- wolf_N_a[t] + wolf_N_p[t]
    
    # summer pup observation error
    wolf_obs_p_sum[t] ~ dpois(max(1e-6, wolf_N_p_sum[t]))
    
    # observation error for the stage-specific counts
    wolf_obs_tot[t] ~ dlnorm(log(wolf_N_tot[t] + 1e-6), wolf_tau_obs)
    wolf_obs_p[t] ~ dlnorm(log(wolf_N_p[t] + 1e-6), wolf_tau_obs)
    wolf_obs_a[t] ~ dlnorm(log(wolf_N_a[t] + 1e-6), wolf_tau_obs)
    
    # prior on detection probability (assumed same for pups + adults)
    wolf_p_det[t] ~ dunif(0, 1)
  }
  
  # CJS using M-arrays
  for (r in 1:(n_years - 1)) { # release occasions
    for (j in 1:(n_years - 1)) { # recaptures
      
      # -------- PUPS -------- 
      wolf_marr_prob_p[r, j] <-
        
        # first recapture is the very next occasion after first capture
        equals(j, r) * # this is 1 when j = r
        wolf_s_p[r] *  # pup survived
        wolf_p_det[r + 1] + # pup was detected
        
        # first recapture is at a later occasion
        step(j - r - 0.5) *  # this is 1 when j > r
        wolf_s_p[r] * # pup survived first interval as a pup
        prod(wolf_s_a[(r + 1):j]) * # then (as an adult by that time) survived for subsequent years
        prod(1 - wolf_p_det[(r + 1):j]) * # was not detected during those years
        wolf_p_det[j + 1] # but then was eventually detected
      
      # -------- ADULTS -------- 
      wolf_marr_prob_a[r, j] <-
        
        # first recapture is the very next occasion after first capture
        equals(j, r) * # this is 1 when j = r
        wolf_s_a[r] * wolf_p_det[r + 1] + # adult survived
        
        # first recapture is at a later occasion
        step(j - r - 0.5) * # this is 1 when j > r
        prod(wolf_s_a[r:j]) * # adult survived for however many years later
        prod(1 - wolf_p_det[(r + 1):j]) * # but was not detected during those years
        wolf_p_det[j + 1] # but was eventually detected
    }
    
    # individuals that either die or are never again detected after last capture
    wolf_marr_prob_p[r, n_years] <- 1 - sum(wolf_marr_prob_p[r, 1:(n_years - 1)])
    wolf_marr_prob_a[r, n_years] <- 1 - sum(wolf_marr_prob_a[r, 1:(n_years - 1)])
    
    # multinomial distribution
    wolf_marray_p[r, 1:n_years] ~ dmulti(wolf_marr_prob_p[r, 1:n_years], wolf_rel_p[r])
    wolf_marray_a[r, 1:n_years] ~ dmulti(wolf_marr_prob_a[r, 1:n_years], wolf_rel_a[r])
  }
  
  ###############################################################
  ######### ------------- GRIZZLY SUBMODEL ------------ #########
  ###############################################################
  
  # latent grizzly abundance the first year (norm dist with mean = overall mean obs), 
  # just helps model get started
  griz_logN[1] ~ dnorm(griz_logN_init_mean, 1 / 0.5^2)  
  
  # 
  ##
  ###
  # since our elk abundance estimates are the pre-birth-pulse estimates,
  # we don't actually estimate the number of "calves" (instead, we do yearlings).
  # therefore, in order to estimate the number of calves that grizzlies would prey on,
  # we need to know how many total were born.
  # previously, we estimated how many would be expected to be born from young & old adults
  # separately. so now, we just derive total as the sum (deterministically)

  # total calves born in each year (available to be eaten by grizzlies)
  for (t in 1:(n_years - 1)) {
    elk_calves_born[t] <- elk_f_ya[t] * elk_N_ya[t] +
      elk_f_oa[t] * elk_N_oa[t]
  }
  ###
  ##
  #
  
  for (t in 2:n_years) {
    
    # for now, expected abundance in year t depends on abundance in year t - 1 + some effect of elk calves
    griz_mu[t] <-
      griz_logN[t - 1] +
      beta1_griz_elkCalves * log(elk_calves_born[t - 1] + 1e-6) # + 
    #   beta2_griz_x2 * griz_x2_std[t - 1] +
    #   beta3_griz_x3 * griz_x3_std[t - 1] # leaving these to make it easier to add more covariates as needed

    # demographic stochasticity
    griz_logN[t] ~ dnorm(griz_mu[t], griz_tau_proc)
  }
  
  # convert latent log-abundance back to natural scale
  for (t in 1:n_years) {
    griz_N[t] <- exp(griz_logN[t])
  }
  
  # observation error prior
  griz_sigma_obs ~ dunif(0.05, 2)
  griz_tau_obs <- 1 / (griz_sigma_obs^2)
  
  for (t in 1:n_years) {
    griz_obs[t] ~ dlnorm(griz_logN[t], griz_tau_obs)
  }
  
  # process error prior
  griz_sigma_proc ~ dunif(0.01, 1)
  griz_tau_proc <- 1 / (griz_sigma_proc^2)
  
  # standardize grizzly abundance for use in elk regressions
  for (t in 1:n_years) {
    griz_N_std[t] <- (griz_N[t] - griz_N_mean) / griz_N_sd
  }
  
  # priors for grizzly process covariates (uncomment as added)
  beta1_griz_elkCalves ~ dnorm(0, 1 / 1^2)
  # beta2_griz_x2 ~ dnorm(0, 1 / 0.3^2)
  # beta3_griz_x3 ~ dnorm(0, 1 / 0.3^2)
  
  ###############################################################
  ######### -------------- BISON SUBMODEL ------------- #########
  ###############################################################

  # latent bison abundance the first year (norm dist with mean = overall mean obs), 
  # just helps model get started
  bison_logN[1] ~ dnorm(bison_logN_init_mean, 1 / 0.5^2)  
  
  # state-space 
  for (t in 2:n_years) {
    
    # for now, expected abundance in year t depends on abundance in year t - 1 + some effect of bison harvest/cull
    bison_mu[t] <-
      bison_logN[t - 1] +
      beta1_bison_cull * log(bison_culled[t - 1] + 1e-6)# + 
    #   beta2_bison_x2 * bison_x2_std[t - 1] +
    #   beta3_bison_x3 * bison_x3_std[t - 1] # leaving these to make it easier to add more covariates as needed
    
    bison_logN[t] ~ dnorm(bison_mu[t], bison_tau_proc)
  }
  
  # convert latent log-abundance back to natural scale
  for (t in 1:n_years) {
    bison_N[t] <- exp(bison_logN[t])
  }
  
  # observation error prior
  bison_sigma_obs ~ dunif(0.05, 2)
  bison_tau_obs <- 1 / (bison_sigma_obs^2)
  
  for (t in 1:n_years) {
    bison_obs[t] ~ dlnorm(bison_logN[t], bison_tau_obs)
  }
  
  # process error prior
  bison_sigma_proc ~ dunif(0.01, 1)
  bison_tau_proc <- 1 / (bison_sigma_proc^2)
  
  # standardize bison abundance for use in wolf regressions
  for (t in 1:n_years) {
    bison_N_std[t] <- (bison_N[t] - bison_N_mean) / bison_N_sd
  }
  
  # priors for bison process covariates (uncomment as added)
  beta1_bison_cull ~ dnorm(0, 1 / 1^2)
  # beta2_bison_x2 ~ dnorm(0, 1 / 0.3^2)
  # beta3_bison_x3 ~ dnorm(0, 1 / 0.3^2)
  
  ###############################################################
  ######### ------------- COUGAR SUBMODEL ------------- #########
  ###############################################################
  
  # latent cougar abundance the first year (norm dist with mean = overall mean obs), 
  # just helps model get started
  cougar_logN[1] ~ dnorm(cougar_logN_init_mean, 1 / 0.5^2)  
  
  # state-space 
  for (t in 2:n_years) {
    
    # for now, expected abundance in year t depends on abundance in year t - 1 + some effect of elk total abundance
    cougar_mu[t] <-
      cougar_logN[t - 1] +
      beta1_cougar_elk * log(elk_N_female[t - 1] + 1e-6)# + 
    #   beta2_cougar_x2 * cougar_x2_std[t - 1] +
    #   beta3_cougar_x3 * cougar_x3_std[t - 1] # leaving these to make it easier to add more covariates as needed
    
    cougar_logN[t] ~ dnorm(cougar_mu[t], cougar_tau_proc)
  }
  
  # convert latent log-abundance back to natural scale
  for (t in 1:n_years) {
    cougar_N[t] <- exp(cougar_logN[t])
  }
  
  # observation error prior
  cougar_sigma_obs ~ dunif(0.05, 2)
  cougar_tau_obs <- 1 / (cougar_sigma_obs^2)
  
  for (t in 1:n_years) {
    cougar_obs[t] ~ dlnorm(cougar_logN[t], cougar_tau_obs)
  }
  
  # process error prior
  cougar_sigma_proc ~ dunif(0.01, 1)
  cougar_tau_proc <- 1 / (cougar_sigma_proc^2)
  
  # standardize cougar abundance for use in wolf regressions
  for (t in 1:n_years) {
    cougar_N_std[t] <- (cougar_N[t] - cougar_N_mean) / cougar_N_sd
  }
  
  # priors for cougar process covariates (uncomment as added)
  beta1_cougar_elk ~ dnorm(0, 1 / 1^2)
  # beta2_cougar_x2 ~ dnorm(0, 1 / 0.3^2)
  # beta3_cougar_x3 ~ dnorm(0, 1 / 0.3^2)
  
  ###############################################################
  ###### ------- ELK AND WOLF SURVIVAL REGRESSIONS ------- ######
  ###############################################################
  
  # ****** Elk regression priors ****** 
  
  # calves
  beta0_calfSurv ~ dnorm(qlogis(0.22), 1 / 0.3^2)
  beta1_calfSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_calfSurv_winterSeverity ~ dnorm(0, 1 / 0.3^2)
  beta3_calfSurv_grizN ~ dnorm(0, 1 / 0.3^2)
  beta4_calfSurv_cougarN ~ dnorm(0, 1 / 0.3^2)
  beta5_calfSurv_elkN ~ dnorm(0, 1 / 0.3^2)  

  # young adults
  beta0_yaSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2)
  beta1_yaSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_yaSurv_winterSeverity ~ dnorm(0, 1 / 0.3^2)
  beta4_yaSurv_cougarN ~ dnorm(0, 1 / 0.3^2)
  beta5_yaSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  beta6_yaSurv_annualNpp ~ dnorm(0, 1 / 0.3^2)
  beta7_yaSurv_browndown ~ dnorm(0, 1 / 0.3^2)
  beta8_yaSurv_pdsi ~ dnorm(0, 1 / 0.3^2)
  
  # old adults
  beta0_oaSurv ~ dnorm(qlogis(0.80), 1 / 0.3^2)
  beta1_oaSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  beta2_oaSurv_winterSeverity ~ dnorm(0, 1 / 0.3^2)
  beta4_oaSurv_cougarN ~ dnorm(0, 1 / 0.3^2)
  beta5_oaSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  beta6_oaSurv_annualNpp ~ dnorm(0, 1 / 0.3^2)
  beta7_oaSurv_browndown ~ dnorm(0, 1 / 0.3^2)
  beta8_oaSurv_pdsi ~ dnorm(0, 1 / 0.3^2)
  
  # error
  sigma_calf ~ dunif(0, 0.5)
  sigma_ya ~ dunif(0, 0.3)
  sigma_oa ~ dunif(0, 0.3)
  tau_calf <- 1 / (sigma_calf^2)
  tau_ya <- 1 / (sigma_ya^2)
  tau_oa <- 1 / (sigma_oa^2)
  
  # ****** Wolf regression priors ****** 
  
  # pups
  beta0_wpupSurv ~ dnorm(qlogis(0.5), 1 / 0.3^2)
  beta1_wpupSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  beta2_wpupSurv_bisonN ~ dnorm(0, 1 / 0.3^2)
  beta3_wpupSurv_wolfN ~ dnorm(0, 1 / 0.3^2)
  
  # adults
  beta0_wadSurv ~ dnorm(qlogis(0.90), 1 / 0.3^2)
  beta1_wadSurv_elkN ~ dnorm(0, 1 / 0.3^2)
  beta2_wadSurv_bisonN ~ dnorm(0, 1 / 0.3^2)
  beta3_wadSurv_wolfN ~ dnorm(0, 1 / 0.3^2)

  # error
  sigma_wpup ~ dunif(0, 0.5)
  sigma_wad ~ dunif(0, 0.3)
  tau_wpup <- 1 / (sigma_wpup^2)
  tau_wad <- 1 / (sigma_wad^2)
  
  # ****** Pre-regression survival priors ****** 
  
  # we wanted to exclude some of the years (1995-1997) for our regressions, 
  # because the first few years of wolf reintroduction were kind of a shitshow,
  # and there was a lot of human influence that can't really be accounted for. 
  # SO, pre-regression years get direct priors rather than regression structure
  
  for (t in 1:(regression_start_idx - 1)) {
    logit(elk_s_c[t]) ~ dnorm(qlogis(0.22), 1 / 0.5^2)
    logit(elk_s_ya[t]) ~ dnorm(qlogis(0.90), 1 / 0.5^2)
    logit(elk_s_oa[t]) ~ dnorm(qlogis(0.80), 1 / 0.5^2)
    logit(wolf_s_p[t]) ~ dnorm(qlogis(0.50), 1 / 0.5^2)
    logit(wolf_s_a[t]) ~ dnorm(qlogis(0.90), 1 / 0.5^2)
    
    eps_elk_s_c[t] <- 0
    eps_elk_s_ya[t] <- 0
    eps_elk_s_oa[t] <- 0
    eps_wolf_s_p[t] <- 0
    eps_wolf_s_a[t] <- 0
  }
  
  # ****** Standardizing abundances ****** 
  
  # standardize wolf and elk abundance to be used as predictors in each other's models
  for (t in 1:n_years) {
    wolf_N_tot_std[t] <- (wolf_N_tot[t] - wolf_tot_mean) / wolf_tot_sd
    elk_N_female_std[t] <- (elk_N_female[t] - elk_N_female_mean) / elk_N_female_sd
  }
  
  # ****** Regressions ****** 

  for (t in regression_start_idx:n_years) {
    
    # elk error priors
    eps_elk_s_c[t] ~ dnorm(0, tau_calf)
    eps_elk_s_ya[t] ~ dnorm(0, tau_ya)
    eps_elk_s_oa[t] ~ dnorm(0, tau_oa)
    
    # wolf error priors
    eps_wolf_s_p[t] ~ dnorm(0, tau_wpup)
    eps_wolf_s_a[t] ~ dnorm(0, tau_wad)
    
    # elk calf regression
    logit(elk_s_c[t]) <- # calf survival depends on....
      beta0_calfSurv +
      beta1_calfSurv_wolfN * wolf_N_tot_std[t - 1] + # wolf abundance
      beta2_calfSurv_winterSeverity * winterSeverity[t] + # winter precipitation
      beta3_calfSurv_grizN * griz_N_std[t - 1] + # grizzly abundance
      beta4_calfSurv_cougarN * cougar_N_std[t - 1] + # cougar abundance
      beta5_calfSurv_elkN * elk_N_female_std[t - 1] + # density dependence
      eps_elk_s_c[t] # error
    
    # elk young adult regression
    logit(elk_s_ya[t]) <- # young adult survival depends on...
      beta0_yaSurv +
      beta1_yaSurv_wolfN * wolf_N_tot_std[t - 1] + # wolf abundance
      beta2_yaSurv_winterSeverity * winterSeverity[t] + # winter precipitation
      beta4_yaSurv_cougarN * cougar_N_std[t - 1] + # cougar abundance
      beta5_yaSurv_elkN * elk_N_female_std[t - 1] + # density dependence
      beta6_yaSurv_annualNpp * annualNpp_std[t - 1] + # vegetation productivity
      beta7_yaSurv_browndown * browndown_std[t - 1] + # brown-down date (growing season length)
      beta8_yaSurv_pdsi * pdsi_std[t - 1] + # drought index
      eps_elk_s_ya[t] # error
    
    # elk old adult regression
    logit(elk_s_oa[t]) <- # old adult survival depends on...
      beta0_oaSurv +
      beta1_oaSurv_wolfN * wolf_N_tot_std[t - 1] + # wolf abundance
      beta2_oaSurv_winterSeverity * winterSeverity[t] + # winter precipitation
      beta4_oaSurv_cougarN * cougar_N_std[t - 1] + # cougar abundance
      beta5_oaSurv_elkN * elk_N_female_std[t - 1] + # density dependence
      beta6_oaSurv_annualNpp * annualNpp_std[t - 1] + # vegetation productivity
      beta7_oaSurv_browndown * browndown_std[t - 1] + # brown-down date (growing season length)
      beta8_oaSurv_pdsi * pdsi_std[t - 1] + # drought index
      eps_elk_s_oa[t] # error
    
    # wolf pup regression
    logit(wolf_s_p[t]) <- # wolf pup survival depends on...
      beta0_wpupSurv +
      beta1_wpupSurv_elkN * elk_N_female_std[t - 1] + # elk abundance
      beta2_wpupSurv_bisonN * bison_N_std[t - 1] + # bison abundance
      beta3_wpupSurv_wolfN * wolf_N_tot_std[t - 1] + # density dependence
      eps_wolf_s_p[t]
    
    # wolf adult regression
    logit(wolf_s_a[t]) <- # adult survival depends on...
      beta0_wadSurv +
      beta1_wadSurv_elkN * elk_N_female_std[t - 1] + # elk abundance
      beta2_wadSurv_bisonN * bison_N_std[t - 1] + # bison abundance
      beta3_wadSurv_wolfN * wolf_N_tot_std[t - 1] + # density dependence
      eps_wolf_s_a[t]
  }
})

################################################################################
############################ Constants, data, inits ############################
################################################################################

# constants
icm_constants <- list(
  # for the mechanics of the model
  n_years = n_years,
  regression_start_idx = regression_start_idx,
  # wolves
  wolf_tot_mean = mean(wolf_pop$total_abundance, na.rm = TRUE),
  wolf_tot_sd = sd(wolf_pop$total_abundance, na.rm = TRUE),
  # elk
  elk_N_female_mean = mean(elk_dat_n$n_female, na.rm = TRUE),
  elk_N_female_sd = sd(elk_dat_n$n_female, na.rm = TRUE),
  # grizzlies
  griz_N_mean = mean(grizzly$griz_N, na.rm = TRUE),
  griz_N_sd = sd(grizzly$griz_N, na.rm = TRUE),
  griz_logN_init_mean = log(pmax(1, grizzly$griz_N[which(!is.na(grizzly$griz_N))[1]])),
  # bison
  bison_N_mean = mean(bison$NR_Bison, na.rm = TRUE),
  bison_N_sd = sd(bison$NR_Bison, na.rm = TRUE),
  bison_logN_init_mean = log(bison$NR_Bison[which(!is.na(bison$NR_Bison))[1]]),
  bison_culled = bison$total_cull_harvest,
  # cougars
  cougar_N_mean = mean(cougars$cougar_N, na.rm = TRUE),
  cougar_N_sd = sd(cougars$cougar_N, na.rm = TRUE),
  cougar_logN_init_mean = log(pmax(1, cougars$cougar_N[which(!is.na(cougars$cougar_N))[1]]))
)

# data
icm_data <- list(
  # elk
  elk_obs_female = elk_dat_n$n_female,
  elk_marray_ya = elk_ya$marray,
  elk_rel_ya = elk_ya$releases,
  elk_marray_oa = elk_oa$marray,
  elk_rel_oa = elk_oa$releases,
  elk_young_num_preg = elk_dat_fec$young_num_preg,
  elk_young_num_capt = elk_dat_fec$young_num_capt,
  elk_old_num_preg = elk_dat_fec$old_num_preg,
  elk_old_num_capt = elk_dat_fec$old_num_capt,
  elk_CCR_cow_youngadult = elk_dat_fec$n_cows_young,
  elk_CCR_cow_oldadult = elk_dat_fec$n_cows_old,
  elk_harvested_13yo = elk_dat_fec$harvested_age13,
  elk_harvested_ya = elk_dat_fec$harvested_total,
  # wolves
  wolf_obs_tot = wolf_pop$total_abundance,
  wolf_obs_p_sum = wolf_pop$summer_pups,
  wolf_obs_p = wolf_pop$dec_pups,
  wolf_obs_a = wolf_pop$dec_adults,
  wolf_marray_p = wolf_p$marray,
  wolf_rel_p = wolf_p$releases,
  wolf_marray_a = wolf_a$marray,
  wolf_rel_a = wolf_a$releases,
  # covars
  winterSeverity = covars_std$winter_severity,
  annualNpp_std = covars_std$annual_npp,
  browndown_std = covars_std$browndown_onset_greenness_min,
  pdsi_std = covars_std$summer_avg_pdsi,
  bison_obs = bison$NR_Bison,
  griz_obs = grizzly$griz_N,
  cougar_obs = cougars$cougar_N
  # elk_ya_harvest = covars$age_2_13,
  # elk_oa_harvest = covars$age_14_plus
)

# initial values
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
if (n_years >= 2) wolf_init_Np_bio[2] <- max(0, wolf_init_Np[2] - 9)
if (n_years >= 3) wolf_init_Np_bio[3:n_years] <- wolf_init_Np[3:n_years]
wolf_init_Np_bio <- pmin(wolf_init_Np_bio, pmax(0, round(wolf_summer_pups)))
wolf_init_Np_bio[1] <- 0

wolf_init_Na <- pmax(1, round(wolf_init_Ntot - wolf_init_Np))

make_icm_inits <- function() {
  
  bison_init_logN <- log(pmax(1, ifelse(
    is.na(bison$NR_Bison),
    mean(bison$NR_Bison, na.rm = TRUE),
    bison$NR_Bison
  )))
  
  griz_init_logN <- log(pmax(1, ifelse(
    is.na(grizzly$griz_N),
    mean(grizzly$griz_N, na.rm = TRUE),
    grizzly$griz_N
  )))
  
  cougar_init_logN <- log(pmax(1, ifelse(
    is.na(cougars$cougar_N),
    mean(cougars$cougar_N, na.rm = TRUE),
    cougars$cougar_N
  )))
  
  griz_obs_init <- ifelse(
    is.na(grizzly$griz_N),
    mean(grizzly$griz_N, na.rm = TRUE),
    grizzly$griz_N
  )
  
  bison_obs_init <- ifelse(
    is.na(bison$NR_Bison),
    mean(bison$NR_Bison, na.rm = TRUE),
    bison$NR_Bison
  )
  
  list(
    # elk demography
    elk_p_13 = rep(0.15, n_years),
    elk_f_ya = rep(0.76, n_years - 1),
    elk_f_oa = rep(0.64, n_years - 1),
    
    # elk observation error
    elk_sigma_obs_female = 0.30,
    elk_p_det = runif(n_years, 0.6, 0.95),
    
    # elk abundances
    elk_N_1y = elk_init_N1y,
    elk_N_ya = elk_init_Nya,
    elk_N_oa = elk_init_Noa,
    
    # wolf demography
    wolf_f = rep(1.0, n_years - 1),
    
    # wolf observation error
    wolf_sigma_obs = 0.2,
    wolf_p_det = runif(n_years, 0.6, 0.95),
    
    # wolf abundances
    wolf_N_p_sum = pmax(1, round(wolf_summer_pups)),
    wolf_N_p = wolf_init_Np,
    wolf_N_p_bio = wolf_init_Np_bio,
    wolf_N_a = wolf_init_Na,
    
    # grizzly abundance
    griz_logN = griz_init_logN,
    griz_sigma_obs = 0.2,
    griz_sigma_proc = 0.1,
    beta1_griz_elkCalves = 0,
    
    # bison abundances and culls
    bison_logN = bison_init_logN,
    bison_sigma_obs = 0.2,
    bison_sigma_proc = 0.1,
    beta1_bison_cull = 0,
    
    # cougar abundances
    cougar_logN = cougar_init_logN,
    cougar_sigma_obs = 0.2,
    cougar_sigma_proc = 0.1,
    beta1_cougar_elk = 0,
    
    # elk survival covariates
    beta0_calfSurv = qlogis(0.22),
    beta1_calfSurv_wolfN = 0,
    beta2_calfSurv_winterSeverity = 0,
    beta3_calfSurv_grizN = 0,
    beta4_calfSurv_cougarN = 0,
    beta5_calfSurv_elkN = 0,
    
    beta0_yaSurv = qlogis(0.90),
    beta1_yaSurv_wolfN = 0,
    beta2_yaSurv_winterSeverity = 0,
    beta4_yaSurv_cougarN = 0,
    beta5_yaSurv_elkN = 0,
    beta6_yaSurv_annualNpp = 0, 
    beta7_yaSurv_browndown = 0, 
    beta8_yaSurv_pdsi = 0,
    
    beta0_oaSurv = qlogis(0.80),
    beta1_oaSurv_wolfN = 0,
    beta2_oaSurv_winterSeverity = 0,
    beta4_oaSurv_cougarN = 0,
    beta5_oaSurv_elkN = 0,
    beta6_oaSurv_annualNpp = 0, 
    beta7_oaSurv_browndown = 0, 
    beta8_oaSurv_pdsi = 0,
    
    sigma_calf = 0.1,
    sigma_ya = 0.1,
    sigma_oa = 0.1,
    eps_elk_s_c = c(0, rep(0, n_years - 1)),
    eps_elk_s_ya = c(0, rep(0, n_years - 1)),
    eps_elk_s_oa = c(0, rep(0, n_years - 1)),
    
    # wolf survival covariates
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
    eps_wolf_s_a = c(0, rep(0, n_years - 1))
  )
}

# construct the inits for use in WAIC calculations later on
waic_inits <- make_icm_inits()

# params to monitor
icm_params <- c(
  # elk demography
  "elk_s_c", "elk_s_ya", "elk_s_oa", "elk_p_13",
  "elk_f_ya", "elk_f_oa", "elk_p_det",
  # elk abundances
  "elk_N_1y", "elk_N_ya", "elk_N_oa", "elk_N_female",
  # wolf demography
  "wolf_s_p", "wolf_s_a", "wolf_f", "wolf_p_det",
  # wolf abundances
  "wolf_N_p_sum", "wolf_N_p", "wolf_N_a", "wolf_N_tot",
  # grizzly abundances and vars for state-space model
  "griz_N", "griz_logN", "griz_sigma_obs", "griz_sigma_proc", "beta1_griz_elkCalves",
  "griz_mu", "elk_calves_born", "griz_N_std",
  # bison abundances and vars for state-space model
  "bison_N", "bison_logN", "bison_mu",
  "bison_sigma_obs", "bison_sigma_proc",
  "beta1_bison_cull",
  # cougar abundances and vars for state-space
  "cougar_N", "cougar_logN", "cougar_mu",
  "cougar_sigma_obs", "cougar_sigma_proc",
  "beta1_cougar_elk", "cougar_N_std",
  
  # elk survival regression covariates
  "beta0_calfSurv", 
  "beta1_calfSurv_wolfN", 
  "beta2_calfSurv_winterSeverity", 
  "beta3_calfSurv_grizN",
  "beta4_calfSurv_cougarN",
  "beta5_calfSurv_elkN",
  
  "beta0_yaSurv", 
  "beta1_yaSurv_wolfN", 
  "beta2_yaSurv_winterSeverity", 
  "beta4_yaSurv_cougarN",
  "beta5_yaSurv_elkN",
  'beta6_yaSurv_annualNpp', 
  'beta7_yaSurv_browndown',
  'beta8_yaSurv_pdsi',
  
  "beta0_oaSurv", 
  "beta1_oaSurv_wolfN", 
  "beta2_oaSurv_winterSeverity", 
  "beta4_oaSurv_cougarN",
  "beta5_oaSurv_elkN",
  'beta6_oaSurv_annualNpp', 
  'beta7_oaSurv_browndown',
  'beta8_oaSurv_pdsi',
  
  # elk errors
  "sigma_calf", "sigma_ya", "sigma_oa", "eps_elk_s_c", "eps_elk_s_ya", "eps_elk_s_oa",
  
  # wolf survival regression covariates
  "beta0_wpupSurv", 
  "beta1_wpupSurv_elkN", 
  "beta2_wpupSurv_bisonN", 
  "beta3_wpupSurv_wolfN",
  
  "beta0_wadSurv", 
  "beta1_wadSurv_elkN", 
  "beta2_wadSurv_bisonN", 
  "beta3_wadSurv_wolfN",
  
  # wolf errors
  "sigma_wpup", "sigma_wad", "eps_wolf_s_p", "eps_wolf_s_a",
  
  # all remaining parameters must be monitored to retroactively calculate WAIC
  'logit_elk_p_13', 'elk_sigma_obs_female', 'wolf_sigma_obs', 'logit_elk_s_ya', 
  'logit_elk_s_oa', 'logit_wolf_s_p', 'logit_wolf_s_a', 'wolf_N_p_bio'
)

################################################################################
############################ Parallel model fitting ############################
################################################################################

# function to run a single MCMC chain
# this is the unit of work that will be sent to each CPU core
run_one_chain <- function(chain_id,
                          icm_code,
                          icm_data,
                          icm_constants,
                          icm_params,
                          make_icm_inits,
                          ni,
                          nb,
                          th) {
  
  # required packages must be reloaded on each core 
  # because each parallel worker starts a fresh R session
  library(nimble)
  library(coda)
  
  # each chain gets a different random seed
  set.seed(1000 + chain_id)
  
  # running a single chain NIMBLE MCMC
  nimbleMCMC(
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
}

# MCMC settings
set.seed(17)
nc <- 3
ni <- 200000
nb <- 40000
th <- 4

# stamp start time to calculate total runtime after run
start_time <- Sys.time()

# make parallel cluster with one worker per chain
cl <- makeCluster(nc)

# load required packages on every worker in the cluster
clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  NULL
})

# export all objects needed by the workers
# these are copied from the main R session to each worker session
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
    "wolf_init_Na",
    "grizzly",
    "bison",
    "cougars"
  ),
  envir = environment()
)

# run each chain in parallel
chain_samples <- parLapply(cl, 1:nc, function(chain_id) {
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
})

# shut down the cluster after all chains are complete
# this is important so background worker processes don't remain open
stopCluster(cl)

# stamp model end time and calculate total run time
end_time <- Sys.time()
run_time <- end_time - start_time
print(paste0("Model runtime: ", round(run_time, 2), " ", units(run_time)))

# combine list of individual chain outputs into one mcmc.list object
# (req. for MCMCsummary() to work later)
icm_samples <- mcmc.list(chain_samples)

# convert all posterior draws across chains into one matrix
post_mat <- as.matrix(icm_samples)

# keep only columns that are finite for every posterior draw 
# (some end up being NA which breaks summary functions later)
good_cols <- colnames(post_mat)[apply(post_mat, 2, function(x) all(is.finite(x)))]

# rebuild a clean mcmc.list using only the finite parameters
icm_samples_clean <- mcmc.list(lapply(chain_samples, function(ch) {
  mcmc(as.matrix(ch)[, good_cols, drop = FALSE])
}))

# summarize posterior distributions across the cleaned chains
# returns means, credible intervals, Rhat, and ESS
icm_summary <- MCMCsummary(icm_samples_clean)

# write data
# stop('The following line will overwrite data. Are you sure you would like to proceed?')
save(
  
  # posterior draws
  icm_samples,
  icm_samples_clean,
  
  # objects required to reconstruct the model for offline WAIC
  icm_code,
  icm_data,
  icm_constants,
  waic_inits,
  
  # useful supporting output
  icm_summary,
  icm_params,
  community_years,
  elk_dat_n,
  wolf_pop,
  grizzly,
  bison,
  cougars,
  covars,
  run_time,
  drop_regression_years,
  n_years,
  
  file = "data/outputs/ICM_parallel_output_2026-07-08.RData"
  
)

################################################################################
################################################################################
################################################################################