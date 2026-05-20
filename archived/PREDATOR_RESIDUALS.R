### Predator residual calculations
### Last updated: Apr. 15, 2026
### xprockox@gmail.com

################################################################################
############------------------ Packages -------------------#####################
################################################################################

library(tidyverse)

################################################################################
############----------------- Data import -----------------#####################
################################################################################

dat <- read.csv('data/allSpp_Abundances.csv')

################################################################################
############---------------- Griz ~ Wolves ----------------#####################
################################################################################

griz_dat <- dat %>%
  select(Year, Wolves, Grizzly.Bears) %>%
  filter(Year %in% 1995:2025) %>%
  drop_na()

griz_lm <- lm(Grizzly.Bears ~ Wolves, data = griz_dat)

griz_dat <- griz_dat %>%
  mutate(
    grizzly_resid = resid(griz_lm)
  )

griz_dat

ggplot() +
  geom_point(data = griz_dat, aes(x = Wolves, y = Grizzly.Bears)) +
  geom_smooth(data = griz_dat, aes(x = Wolves, y = Grizzly.Bears), method = "lm") +
  geom_text(
    data = griz_dat,
    aes(x = Wolves, y = Grizzly.Bears, label = Year),
    nudge_y = 8,
    size = 3
  ) +
  theme_bw()

################################################################################
############----------- Cougars ~ Wolves + Griz -----------#####################
################################################################################

for (i in 1:nrow(dat)){
  dat$coug_lambda[i] <- dat$Cougars[i+1] / dat$Cougars[i] 
}

coug_lambda_sims <- rnorm(100000, 
                          mean = mean(dat$coug_lambda, na.rm = TRUE),
                          sd = sd(dat$coug_lambda, na.rm = TRUE))

dat$coug_lambda_sim <- sample(coug_lambda_sims, size = nrow(dat))

dat$coug_sim <- NA

for (i in 1:nrow(dat)){
  dat$coug_sim[i] <- ifelse(is.na(dat$Cougars)==TRUE, dat$Cougars[i-1] * dat$coug_lambda[i], dat$Cougars[i])
}


coug_dat <- dat %>%
  select(Year, Wolves, Cougars, Grizzly.Bears) %>%
  filter(Year %in% 1995:2025) %>%
  drop_na()



coug_lm <- lm(Cougars ~ Wolves + Grizzly.Bears, data = coug_dat)

coug_dat <- coug_dat %>%
  mutate(
    cougar_resid = resid(coug_lm)
  )

coug_dat

ggplot() +
  geom_point(data = coug_dat, aes(x = Wolves, y = Grizzly.Bears)) +
  geom_smooth(data = coug_dat, aes(x = Wolves, y = Grizzly.Bears), method = "lm") +
  geom_text(
    data = coug_dat,
    aes(x = Wolves, y = Grizzly.Bears, label = Year),
    nudge_y = 8,
    size = 3
  ) +
  theme_bw()