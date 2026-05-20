### Elk IPM Main Script
### Last updated: Nov. 14, 2025
### Contact: xprockox@gmail.com

############################################################################################
################# --------------- PACKAGES AND SET-UP ---------------- #####################
############################################################################################
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(MCMCvis)
library(ggplot2)
library(nimble)
library(coda)

############################################################################################
#################### --------------- DATA IMPORT ---------------- ##########################
############################################################################################
load('data/intermediate/adultSurvival_cjsMatrices.rData')
dat_n <- read.csv('data/intermediate/abundanceEstimates_stages.csv')
dat_fec <- read.csv('data/intermediate/fecundity.csv')

## 2025-10-30: Found an issue where there are latent state histories that are incorrect:
## e.g., one individual in z could be 1, 1, 1, 1, NA, NA, NA, 1, 1
## but if we know they're alive at the 8th and 9th timesteps, those NAs should also be 1s

z_fixed <- z # make a copy to work on

# fix the copy
for (i in 1:nrow(z)) {
  first_det <- which(z[i, ] == 1)[1]
  last_det <- max(which(z[i, ] == 1))
  
  if (!is.na(first_det) && !is.na(last_det) && first_det < last_det) {
    # Fill all NAs between first and last detection with 1
    z_fixed[i, first_det:last_det] <- 1
  }
}

# Use the fixed version
z <- z_fixed
rm(z_fixed)

# we have different years in the CJS and abundance data, so for now let's
# only use years that are shared across both datasets
shared_years <- intersect(as.numeric(dat_n$year),as.numeric(colnames(z)))
dat_n <- dat_n[dat_n$year %in% shared_years,]

colnames(is_class1) <- colnames(z)
colnames(is_class2) <- colnames(z)

z <- z[,which(colnames(z) %in% shared_years)]
y <- y[,which(colnames(y) %in% shared_years)]
is_class1 <- is_class1[,which(colnames(is_class1) %in% shared_years)]
is_class2 <- is_class2[,which(colnames(is_class2) %in% shared_years)]

# the same is true for dat_fec (has more years), so we need to remove years that aren't shared
dat_fec <- dat_fec[dat_fec$year %in% shared_years,]

### Nov. 5, 2025: The total number of elk (n_total) shouldn't be used because this doesn't exclude bulls
# instead, we want to approximate it as the total number of cows + half the total number of calves
dat_n$n_female <- dat_n$n_cow + (dat_n$n_calf/2)

############################################################################################
#################### --------------- NIMBLE CODE ---------------- ##########################
############################################################################################
# define years
n_years <- length(dat_n$n_calf)   

elk_ipm <- nimbleCode({
  
  ## -----------------------------
  ## (1) STATE-SPACE IPM 
  ## -----------------------------
  # priors 
  for (t in 1:n_years) {
    logit(s_c[t]) ~ dnorm(qlogis(0.22), 1 / 0.5^2)   # prior on calf survival
    logit(s_ya[t]) ~ dnorm(qlogis(0.90), 1 / 0.5^2)  # prior on young adult survival
    logit(s_oa[t]) ~ dnorm(qlogis(0.80), 1 / 0.5^2)  # prior on old adult survival
    logit(p_13[t]) ~ dnorm(qlogis(0.15), 1 / 0.5^2)  # prior on young-to-old growth
  }
  
  # observation error prior + derive precision from SD
  sigma_obs_female ~ dunif(0.05, 2)
  tau_obs_female <- 1 / (sigma_obs_female^2)
  
  # initial expected values of stage-specific abundances
  lambda_init_1y ~ dgamma(11.1, 0.00454)   # mean ≈ 2451 (dat_n$n_1y[1] = 2451)
  lambda_init_ya ~ dgamma(11.1, 0.00118)   # mean ≈ 9418 (dat_n$n_cow_youngadult[1] = 9418)
  lambda_init_oa ~ dgamma(11.1, 0.0111)    # mean ≈ 997 (dat_n$n_cow_oldadult[1] = 997)
  
  lambda_init_female ~ dgamma(11.1, 0.000863) # mean ≈ 12866 (sum of all three prev. values = 12866)
  
  # initial latent stage-specific abundances (from expected values)
  N_1y[1] ~ dpois(lambda_init_1y)
  N_ya[1] ~ dpois(lambda_init_ya)
  N_oa[1] ~ dpois(lambda_init_oa)
  
  N_female[1] ~ dpois(lambda_init_female)
  
  obs_female[1] ~ dlnorm(log(N_female[1] + 1e-6), tau_obs_female)
  
  # process model
  for (t in 1:(n_years-1)) {
    # expected values
    mu_1y[t+1] <- f_ya[t] * s_c[t] * N_ya[t] + f_oa[t] * s_c[t] * N_oa[t]
    mu_ya[t+1] <- s_ya[t] * N_1y[t] + s_ya[t] * (1 - p_13[t]) * N_ya[t]
    mu_oa[t+1] <- s_ya[t] * p_13[t] * N_ya[t] + s_oa[t] * N_oa[t]
    
    # latent states
    N_1y[t+1] ~ dpois(max(1e-6, mu_1y[t+1]))
    N_ya[t+1] ~ dpois(max(1e-9, mu_ya[t+1]))
    N_oa[t+1] ~ dpois(max(1e-9, mu_oa[t+1]))
    
    # derive latent total female abundance from latent stage-specific abundances 
    N_female[t+1] <- N_1y[t+1] + N_ya[t+1] + N_oa[t+1]
    
    # observations of total female abundance are related to latent total female abundance
    obs_female[t+1] ~ dlnorm(log(N_female[t+1] + 1e-6), tau_obs_female)
  }
  
  ## -----------------------------
  ## (2) CJS MODEL
  ## -----------------------------
  
  # detection prob
  for (t in 1:n_years){
    p[t] ~ dunif(0, 1)
  }
  
  
  for (i in 1:N) {
    # Time 1
    z[i, 1] ~ dbern(equals(1, first_seen[i]))
    
    # Times 2 onward
    for (t in 2:n_years) {
      # Ensure at least one class is active (add small constant)
      phi[i, t] <- is_class1[i, t-1] * s_ya[t-1] + 
        is_class2[i, t-1] * s_oa[t-1] +
        1e-10  # prevent exactly 0
      
      z[i, t] ~ dbern(
        equals(t, first_seen[i]) +
          step(t - first_seen[i] - 0.5) * (1 - equals(t, first_seen[i])) *
          z[i, t-1] * phi[i, t]
      )
    }
    
    # observations
    for (t in 1:n_years) {
      y[i, t] ~ dbern(p[t] * z[i, t])
    }
  }
  
  ## -----------------------------
  ## (3) CALF SURVIVAL 
  ## -----------------------------
  
  for (t in 1:(n_years-1)) {
    CCR_c_fromYoungCows[t] ~ dbin(f_ya[t] * s_c[t], CCR_cow_youngadult[t])
    CCR_c_fromOldCows[t] ~ dbin(f_oa[t] * s_c[t], CCR_cow_oldadult[t])
    CCR_c[t] <- CCR_c_fromYoungCows[t] + CCR_c_fromOldCows[t]
  }
  
  ## -----------------------------
  ## (4) FECUNDITY 
  ## -----------------------------

  for (t in 1:(n_years-1)){
    f_ya[t] ~ dbeta(1, 1)
    f_oa[t] ~ dbeta(1, 1)
    
    young_num_preg[t] ~ dbin(f_ya[t], young_num_capt[t])
    old_num_preg[t] ~ dbin(f_oa[t], old_num_capt[t])
  }
  
  ## -----------------------------
  ## (5) GROWTH 
  ## -----------------------------
  for (t in 1:(n_years)){
    harvested_13yo[t] ~ dbinom(p_13[t], harvested_ya[t])
  }
})

############################################################################################
########### --------------- MODEL SPECS, INITS, AND DATA SOURCES ---------------- ##########
############################################################################################
# constants & data 
N <- nrow(y)

elk_constants <- list(n_years = n_years, N = N)

elk_data <- list(
  # state-space
  obs_female = dat_n$n_female,  # female-only
  # CJS
  y = y,
  is_class1 = is_class1,
  is_class2 = is_class2,
  first_seen = first_seen,
  # fecundity
  young_num_preg = dat_fec$young_num_preg,
  young_num_capt = dat_fec$young_num_capt,
  old_num_preg = dat_fec$old_num_preg,
  old_num_capt = dat_fec$old_num_capt,
  # calf survival
  CCR_cow_youngadult = dat_fec$n_cows_young,
  CCR_cow_oldadult = dat_fec$n_cows_old,
  # growth
  harvested_13yo = dat_fec$harvested_age13,
  harvested_ya = dat_fec$harvested_total
)


## -----------------------------
## inits  
## -----------------------------

make_inits <- function() {
  # seed latent states near observations
  init_N1y <- ifelse(is.na(dat_n$n_calf), # is n_calf is NA for a given year,
                    pmax(1, round(mean(dat_n$n_calf, na.rm = TRUE))), # use the mean n_calf as initial, with na.rm
                    pmax(1, round(dat_n$n_calf))) # otherwise, use the "observed value"
  init_Nya <- ifelse(is.na(dat_n$n_cow_youngadult),
                    pmax(1, round(mean(dat_n$n_cow_youngadult, na.rm = TRUE))),
                    pmax(1, round(dat_n$n_cow_youngadult)))
  init_Noa <- ifelse(is.na(dat_n$n_cow_oldadult),
                    pmax(1, round(mean(dat_n$n_cow_oldadult, na.rm = TRUE))),
                    pmax(1, round(dat_n$n_cow_oldadult)))
  init_Nfemale <- ifelse(is.na(dat_n$n_female),
                    pmax(1, round(mean(dat_n$n_female, na.rm = TRUE))),
                    pmax(1, round(dat_n$n_female)))
  
  z_init <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in 1:nrow(y)) {
    detections <- which(y[i, ] == 1)
    if (length(detections) > 0) {
      first_det <- min(detections)
      last_det <- max(detections)
      # Set z=1 from first detection through last detection
      z_init[i, first_det:last_det] <- 1L
      # Everything before first detection stays NA (or set to 0)
      if (first_det > 1) {
        z_init[i, 1:(first_det-1)] <- 0L
      }
    }
  }
  
  list(
    # vital rates
    s_c = rep(0.22, n_years),
    s_ya = rep(0.90, n_years),
    s_oa = rep(0.80, n_years),
    p_13 = rep(0.15, n_years),
    
    # observation SDs
    sigma_obs_female = 0.30,
    
    # initial abundances
    lambda_init_1y = max(1, round(init_N1y[1])),
    lambda_init_ya = max(1, round(init_Nya[1])),
    lambda_init_oa = max(1, round(init_Noa[1])),
    lambda_init_female = max(1, round(init_Nfemale[1])),
    N_1y = pmax(1, init_N1y),
    N_ya = pmax(1, init_Nya),
    N_oa = pmax(1, init_Noa),
    N_female = pmax(1, init_Nfemale),
    
    # detection probability for cjs model
    p = runif(n_years, 0.6, 0.95),
    
    # initial z matrix for cjs model
    z = z_init,
    
    # fecundity initials
    f_ya = rep(0.76, n_years - 1),  # mean(dat_fec$young_prop_pregnant, na.rm=TRUE) = 0.76
    f_oa = rep(0.64, n_years - 1)  # mean(dat_fec$old_prop_pregnant, na.rm=TRUE) = 0.64
  )
}


## -----------------------------
## parameters to monitor
## -----------------------------
params <- c(
  # yearly vital rates (shared by IPM & CJS)
  "s_c", "s_ya","s_oa","p_13",
  "f_ya", 'f_oa',

  # detection (CJS)
  "p",
  
  # latent states (and female)
  "N_1y","N_ya","N_oa","N_female"
)

## -----------------------------
## model specs
## -----------------------------
set.seed(17)
nc <- 3
ni <- 1000000
nb <- 200000
th <- 4

############################################################################################
####################### --------------- RUN MODEL ---------------- #########################
############################################################################################

# run MCMC
elk_mod1 <- nimbleMCMC(
  code      = elk_ipm,
  data      = elk_data,
  constants = elk_constants,
  inits     = make_inits,
  monitors  = params,
  nchains   = nc,
  niter     = ni,
  nburnin   = nb,
  thin      = th,
  summary   = TRUE
)

# SAVE OUTPUT
# stop('The following line will overwrite data. Are you sure you would like to proceed?')
# save.image('environment_2026-02-17.RData')

# OR IMPORT PREVIOUSLY RUN MODEL TO WORK WITH RESULTS BEYOND HERE
load('data/results/elkIPM_environment_2026-02-17.RData')

## -----------------------------
## quick summary table
## -----------------------------
round(elk_mod1$summary$all.chains, 2)

############################################################################################

### for some reason, the following lines:

# MCMCtrace(elk_mod1,
#           pdf   = FALSE,
#           ind   = TRUE,
#           Rhat  = TRUE,
#           n.eff = TRUE,
#           params = 'all')
# 
# 
# MCMCsummary(elk_mod1,
#             params = 'all')

### produce the following error:

# Error in coda::mcmc.list(lapply(object, function(x) coda::mcmc(x))) : 
# Different start, end or thin values in each chain

# so let's pull the chains out one-by-one, trim them to be equal, and then stick them back together
############################################################################################

# -----------------------------
# A) Normalize chains' shapes
# -----------------------------

# 1) Each element of elk_mod1$samples is a chain; convert each of these to a plain matrix
mats <- lapply(elk_mod1$samples, function(ch) as.matrix(ch))

# 2) Trim all chains to the same number of iterations (using the minimum length of chain)
lens <- sapply(mats, nrow)
L <- min(lens)
mats_trim <- lapply(mats, function(M) tail(M, L))

# 3) Rebuild a consistent mcmc.list
mlist <- mcmc.list(lapply(mats_trim, function(M) mcmc(M, start = 1, end = L, thin = 1)))

# 4) create an alias (copy of the data) for clarity in the next step
ml <- mlist  

# -----------------------------------------
# B) Drop parameters that contain any NA/NaN
# -----------------------------------------

# Convert each chain back to a matrix for NA screening
mats <- lapply(ml, as.matrix)

# Identify columns (parameters) that have no NA/NaN in ANY chain.
# We compute, per chain, the set of columns with zero NA/NaN, then intersect across chains
keep_cols <- Reduce(intersect, lapply(mats, function(M) {
  colnames(M)[colSums(is.na(M) | is.nan(M)) == 0]
}))

# Build a "clean" mcmc.list that retains only the safe parameters (same columns in every chain)
# Also set consistent start/end/thin to keep MCMCvis happy.
ml_clean <- mcmc.list(lapply(mats, function(M) {
  mcmc(M[, keep_cols, drop = FALSE], start = 1, end = nrow(M), thin = 1)
}))

# -----------------------------------------
# C) Summarize and visualize with MCMCvis
# -----------------------------------------
round(MCMCsummary(ml_clean, params = 'all'),2)
MCMCtrace(ml_clean, pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, params = 'all')

############################################################################################

# Get summaries for all N parameters
N_summ <- MCMCsummary(ml_clean, params=c('N_1y', 'N_ya', 'N_oa', 'N_female'))   # extracts mean, sd, 2.5%, 50%, 97.5%

# Convert rownames to columns
N_summ <- N_summ %>%
  tibble::rownames_to_column("param") %>%
  mutate(
    stage = str_extract(param, "^N_[a-zA-Z0-9_]+"),  # Added 0-9 and _
    t     = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
  ) %>%
  select(stage, t, mean = mean, low = `2.5%`, high = `97.5%`) %>%
  arrange(stage, t)

years <- shared_years   # or however many years you have
N_summ <- N_summ %>% mutate(year = years[t])

N_summ <- N_summ %>%
  mutate(stage = recode(stage,
                        "N_1y"   = "Yearling",
                        "N_ya"   = "Young Adult",
                        "N_oa"   = "Old Adult",
                        "N_female" = "Total Females"))


dat_long <- dat_n %>%
  pivot_longer(
    cols = -c(X,year),          # everything except year becomes stage/value pairs
    names_to = "stage",
    values_to = "value"
  ) %>%
  mutate(stage = recode(stage,
                        n_calf  = "Yearling",
                        n_cow_youngadult = "Young Adult",
                        n_cow_oldadult   = "Old Adult",
                        n_female = "Total Females"))

dat_long$stage <- factor(dat_long$stage, levels = unique(N_summ$stage))

dat_long <- dat_long[dat_long$stage %in% c('Yearling', 'Young Adult', 'Old Adult', 'Total Females'),]

validation_plot <- ggplot(N_summ, aes(x = year, y = mean, group = stage)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill=stage), alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(data = dat_long[dat_long$stage=='Total Females',], aes(y = value), color = "red", size = 2) +
  geom_line(data = dat_long[dat_long$stage=='Total Females',], aes(y = value), color = "red", linetype = 2) +
  facet_wrap(~ stage, scales = "free_y") +
  theme_bw() +
  labs(x = "Year", y = "Abundance",
       title = "Posterior Population Estimates with Validation Data",
       subtitle = "Ribbon = 95% credible interval, Line = posterior mean, Red = observed")+
  theme(legend.position='none')

validation_plot

############################################################################################
### what about vital rates?

vrates <- MCMCsummary(ml_clean,
                      params = c("s_c","s_ya","s_oa","p_13","f_ya","f_oa")) %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  rename(mean = mean, low = `2.5%`, high = `97.5%`)

vrates2 <- vrates %>%
  mutate(
    year_index = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    rate = str_extract(param, "^[^\\[]+")  # remove bracket indices
  )

vrates2$rate <- factor(vrates2$rate,
                      levels = c("s_c","s_ya","s_oa","p_13","f_ya","f_oa"),
                      labels = c("Calf survival (s_c)",
                                 "Young Adult survival (s_ya)",
                                 "Old Adult survival (s_oa)",
                                 "Young→old transition (p_13)",
                                 "Fecundity (young) (f_y)",
                                 "Fecundity (old) (f_o)"))
vrates2_yearchunk1 <- rep(shared_years, 4)
vrates2_yearchunk2 <- rep(shared_years[-length(shared_years)], 2)
vrates2$year <- c(vrates2_yearchunk1, vrates2_yearchunk2)

vrate_plot <- ggplot(vrates2, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = rate), alpha = 0.2) +
  geom_line(size = 0.9) +
  facet_wrap(~ rate, scales = "free_y") +
  theme_minimal() +
  labs(x = "Year", y = "Estimated value",
       title = "Posterior Time-Varying Vital Rates (95% Credible Intervals)") + 
  theme(legend.position = "none")

vrate_plot

stop('The following line will overwrite data. Are you sure you would like to proceed?')

save.image('environment_2024-11-06.RData')

write.csv(vrates2, 'data/results/elk_vrates_2026-02-16.csv')
write.csv(N_summ, 'data/results/elk_N_byStages_2026-02-16.csv')
