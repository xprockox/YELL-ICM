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
nr_pop <- read.csv('data/intermediate/nr_pop.csv')
full_pop <- read.csv('data/intermediate/full_park_pop.csv')

# we have different years in the CJS and abundance data, so for now let's
# only use years that are shared across both datasets
shared_years <- intersect(as.numeric(nr_pop$seasonal.year),as.numeric(colnames(z)))
nr_pop <- nr_pop[nr_pop$seasonal.year %in% shared_years,]

colnames(is_class1) <- colnames(z)
colnames(is_class2) <- colnames(z)

z <- z[,which(colnames(z) %in% shared_years)]
y <- y[,which(colnames(y) %in% shared_years)]
is_class1 <- is_class1[,which(colnames(is_class1) %in% shared_years)]
is_class2 <- is_class2[,which(colnames(is_class2) %in% shared_years)]

# Rebuild first_seen to match the new years:
first_seen <- apply(y, 1, function(row) {
  first <- which(row == 1)[1]
  if (is.na(first)) return(ncol(y)) else return(first)
})
first_seen <- as.integer(first_seen)

############################################################################################
#################### --------------- NIMBLE CODE ---------------- ##########################
############################################################################################
# define years
n_years <- nrow(nr_pop)

wolf_ipm <- nimbleCode({
  
  ## -----------------------------
  ## (1) PRIORS
  ## -----------------------------
  for (t in 1:n_years) {
    logit(s_p[t]) ~ dnorm(qlogis(0.5), 1 / 0.5^2)
    logit(s_a[t]) ~ dnorm(qlogis(0.9), 1 / 0.5^2)
  }
  
  for (t in 1:(n_years - 1)) {
    f[t] ~ dbeta(1, 1)
  }
  
  sigma_obs ~ dunif(0.05, 2)
  tau_obs <- 1 / (sigma_obs^2)
  
  ## -----------------------------
  ## (2) INITIAL ADULT ABUNDANCE
  ## -----------------------------
  lambda_init_a ~ dgamma(10, 1)
  N_a[1] ~ dpois(lambda_init_a)
  
  ## -----------------------------
  ## (3) STATE-SPACE MODEL
  ## -----------------------------
  
  # For each year, winter pups come from summer pups surviving to December
  for (t in 1:n_years) {
    N_p[t] ~ dbin(s_p[t], obs_p_sum[t])
    N_tot[t] <- N_p[t] + N_a[t]
    obs_tot[t] ~ dlnorm(log(N_tot[t] + 1e-6), tau_obs)
  }
  
  # Adult dynamics and fecundity
  for (t in 1:(n_years - 1)) {
    
    # adults next year = surviving adults + surviving pups recruited to adult class
    mu_a[t + 1] <- N_p[t] + s_a[t] * N_a[t]
    N_a[t + 1] ~ dpois(max(1e-6, mu_a[t + 1]))
    
    # observed summer pup count next year
    lambda_p_sum[t + 1] <- f[t] * N_a[t]
    obs_p_sum[t + 1] ~ dpois(max(1e-6, lambda_p_sum[t + 1]))
  }
  
  ## -----------------------------
  ## (4) CJS MODEL
  ## -----------------------------
  for (t in 1:n_years) {
    p[t] ~ dunif(0, 1)
  }
  
  for (i in 1:N) {
    
    z[i, 1] ~ dbern(equals(1, first_seen[i]))
    
    for (t in 2:n_years) {
      phi[i, t] <- is_class1[i, t - 1] * s_p[t - 1] +
        is_class2[i, t - 1] * s_a[t - 1] +
        1e-10
      
      z[i, t] ~ dbern(
        equals(t, first_seen[i]) +
          step(t - first_seen[i] - 0.5) *
          (1 - equals(t, first_seen[i])) *
          z[i, t - 1] * phi[i, t]
      )
    }
    
    for (t in 1:n_years) {
      y[i, t] ~ dbern(p[t] * z[i, t])
    }
  }
})

############################################################################################
########### --------------- MODEL SPECS, INITS, AND DATA SOURCES ---------------- ##########
############################################################################################
# constants & data 
N <- nrow(y)

wolf_constants <- list(n_years = n_years, N = N)

wolf_data <- list(
  obs_tot = nr_pop$total_abundance, # total wolves
  obs_p_sum = nr_pop$summer_pups, # summer pups
  # survival stuff:
  y = y,
  is_class1 = is_class1,
  is_class2 = is_class2,
  first_seen = first_seen
)


## -----------------------------
## inits  
## -----------------------------
make_inits <- function() {
  
  init_Ntot <- ifelse(
    is.na(nr_pop$total_abundance),
    max(1, round(mean(nr_pop$total_abundance, na.rm = TRUE))),
    pmax(1, round(nr_pop$total_abundance))
  )
  
  init_Np <- ifelse(
    is.na(nr_pop$summer_pups),
    max(1, round(mean(nr_pop$summer_pups, na.rm = TRUE) * 0.5)),
    pmax(1, round(nr_pop$summer_pups * 0.5))
  )
  
  init_Na <- pmax(1, round(init_Ntot - init_Np))
  
  z_init <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
  
  for (i in 1:nrow(y)) {
    detections <- which(y[i, ] == 1)
    
    if (length(detections) > 0) {
      first_det <- min(detections)
      last_det <- max(detections)
      
      z_init[i, first_det:last_det] <- 1L
      
      if (first_det > 1) {
        z_init[i, 1:(first_det - 1)] <- 0L
      }
      
      if (last_det < ncol(y)) {
        z_init[i, (last_det + 1):ncol(y)] <- 0L
      }
    }
  }
  
  list(
    s_p = rep(0.5, n_years),
    s_a = rep(0.9, n_years),
    f = rep(0.5, n_years - 1),
    sigma_obs = 0.2,
    
    lambda_init_a = max(1, init_Na[1]),
    
    N_p = init_Np,
    N_a = init_Na,
    p = runif(n_years, 0.6, 0.95),
    z = z_init
  )
}

## -----------------------------
## parameters to monitor
## -----------------------------
params <- c(
  # yearly vital rates (shared by IPM & CJS)
  "s_p", "s_a", "f",

  # detection (CJS)
  "p",
  
  # latent states (and female)
  "N_p","N_a"
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
wolf_mod1 <- nimbleMCMC(
  code      = wolf_ipm,
  data      = wolf_data,
  constants = wolf_constants,
  inits     = make_inits,
  monitors  = params,
  nchains   = nc,
  niter     = ni,
  nburnin   = nb,
  thin      = th,
  summary   = TRUE
)

# SAVE OUTPUT
#stop('The following line will overwrite data. Are you sure you would like to proceed?')
# save.image('data/results/wolfIPM_environment_2026-03-25.RData')

# OR IMPORT PREVIOUSLY RUN MODEL TO WORK WITH RESULTS BEYOND HERE
load('data/results/wolfIPM_environment_2026-03-25.RData')

## -----------------------------
## quick summary table
## -----------------------------
round(wolf_mod1$summary$all.chains, 2)

############################################################################################
################### --------------- EXTRACT RESULTS ---------------- #######################
############################################################################################

# -----------------------------
# A) Normalize chains' shapes
# -----------------------------

# each element of wolf_mod1$samples is a chain; convert each to a plain matrix
mats <- lapply(wolf_mod1$samples, function(ch) as.matrix(ch))

# trim all chains to the same number of iterations
lens <- sapply(mats, nrow)
L <- min(lens)
mats_trim <- lapply(mats, function(M) tail(M, L))

# rebuild a consistent mcmc.list
mlist <- mcmc.list(lapply(mats_trim, function(M) mcmc(M, start = 1, end = L, thin = 1)))

# alias for clarity
ml <- mlist

# -----------------------------
# B) Drop parameters with any NA/NaN
# -----------------------------

mats <- lapply(ml, as.matrix)

keep_cols <- Reduce(intersect, lapply(mats, function(M) {
  colnames(M)[colSums(is.na(M) | is.nan(M)) == 0]
}))

ml_clean <- mcmc.list(lapply(mats, function(M) {
  mcmc(M[, keep_cols, drop = FALSE], start = 1, end = nrow(M), thin = 1)
}))

# combine all chains into one posterior matrix
post <- do.call(rbind, lapply(ml_clean, as.matrix))

############################################################################################
################## --------------- ABUNDANCE RESULTS ---------------- ######################
############################################################################################

# -----------------------------
# C) Summaries for N_p and N_a
# -----------------------------

N_summ <- MCMCsummary(ml_clean, params = c("N_p", "N_a")) %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  mutate(
    stage = str_extract(param, "^N_[a-zA-Z0-9_]+"),
    t = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))
  ) %>%
  select(stage, t, mean = mean, low = `2.5%`, high = `97.5%`) %>%
  arrange(stage, t)

# map model index to actual year
shared_years <- as.numeric(colnames(y))
N_summ <- N_summ %>%
  mutate(year = shared_years[t])

# relabel for plotting
N_summ <- N_summ %>%
  mutate(stage = recode(stage,
                        "N_p" = "December pups",
                        "N_a" = "December adults"))

# -----------------------------
# D) Compute N_tot from posterior draws
# -----------------------------
# do this from posterior samples, not by summing summary quantiles

N_a_cols <- grep("^N_a\\[", colnames(post), value = TRUE)
N_p_cols <- grep("^N_p\\[", colnames(post), value = TRUE)

# make sure columns are ordered by time index
extract_index <- function(x) as.integer(str_extract(x, "(?<=\\[)\\d+(?=\\])"))
N_a_cols <- N_a_cols[order(extract_index(N_a_cols))]
N_p_cols <- N_p_cols[order(extract_index(N_p_cols))]

N_tot_post <- post[, N_a_cols, drop = FALSE] + post[, N_p_cols, drop = FALSE]

N_tot_summ <- data.frame(
  t = seq_len(ncol(N_tot_post)),
  mean = apply(N_tot_post, 2, mean),
  low = apply(N_tot_post, 2, quantile, probs = 0.025),
  high = apply(N_tot_post, 2, quantile, probs = 0.975)
) %>%
  mutate(
    year = shared_years[t],
    stage = "December total"
  ) %>%
  select(stage, t, mean, low, high, year)

# combine all abundance summaries
N_all <- bind_rows(N_summ, N_tot_summ)

############################################################################################
################### --------------- OBSERVED DATA ----------------- ########################
############################################################################################

# make observed data long for plotting
obs_long <- nr_pop %>%
  transmute(
    year = seasonal.year,
    `December pups` = dec_pups,
    `December adults` = dec_adults,
    `December total` = total_abundance
  ) %>%
  pivot_longer(
    cols = -year,
    names_to = "stage",
    values_to = "observed"
  )

# match factor order
N_all$stage <- factor(N_all$stage,
                      levels = c("December pups", "December adults", "December total"))

obs_long$stage <- factor(obs_long$stage,
                         levels = c("December pups", "December adults", "December total"))

############################################################################################
#################### --------------- ABUNDANCE PLOT ---------------- #######################
############################################################################################

abundance_plot <- ggplot(N_all, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = stage), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_point(data = obs_long, aes(y = observed), color = "red", size = 1) +
  geom_line(data = obs_long, aes(y = observed), color = "red", linetype = 3) +
  facet_wrap(~ stage, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Abundance",
    title = "Posterior wolf abundance estimates with observed data",
    subtitle = "Black line = posterior mean, Ribbon = 95% credible interval, \nRed points/dashed line = observed data"
  ) +
  theme(legend.position = "none")

abundance_plot

############################################################################################
################### --------------- VITAL RATES ------------------- ########################
############################################################################################

# extract summaries for f, s_a, s_p
vrates <- MCMCsummary(ml_clean, params = c("f", "s_a", "s_p")) %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  rename(mean = mean, low = `2.5%`, high = `97.5%`)

vrates2 <- vrates %>%
  mutate(
    year_index = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
    rate = str_extract(param, "^[^\\[]+")
  )

# map to years
# s_a and s_p have length n_years
# f usually has length n_years - 1
vrates2 <- bind_rows(
  vrates2 %>%
    filter(rate %in% c("s_a", "s_p")) %>%
    mutate(year = shared_years[year_index]),
  
  vrates2 %>%
    filter(rate == "f") %>%
    mutate(year = shared_years[year_index])
)

# relabel
vrates2$rate <- factor(vrates2$rate,
                       levels = c("s_p", "s_a", "f"),
                       labels = c("Pup survival (s_p)",
                                  "Adult survival (s_a)",
                                  "Fecundity (f)"))

vital_rate_plot <- ggplot(vrates2, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = rate), alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ rate, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Estimate",
    title = "Posterior wolf vital rates",
    subtitle = "Black line = posterior mean, Ribbon = 95% credible interval"
  ) +
  theme(legend.position = "none")

vital_rate_plot

############################################################################################
################## --------------- OPTIONAL EXPORT ---------------- ########################
############################################################################################
stop('WARNING: The following lines will overwrite data. 
     Are you sure you would like to proceed?')

write.csv(N_all, "data/results/wolf_N_summary.csv", row.names = FALSE)
write.csv(vrates2, "data/results/wolf_vital_rates_summary.csv", row.names = FALSE)