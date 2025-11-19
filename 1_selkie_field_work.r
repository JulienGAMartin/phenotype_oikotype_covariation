################################################################################
# Field work (simulating data) for selkie case study
# Author: Julien Martin
# Date: 18 November 2024
###############################################################################

# loading packages
library(mvtnorm)
library(dplyr)

# setting seed for reproducibility
set.seed(19812012)

###############################################################################
# simulated model
###############################################################################

# Phenotypic trait: aggressivity
# y = bm0 + um0 + (bmh + sm1) * h + bm1 * sex + e
# sd(e) = exp(bv0 + uv0 + bv1 * sex)

# sd(e) = exp(bv0 + uv0 + bv1 * sex)

# Oikotypic trait: temperature
# h = bm0 + um0 + e
# sd(e) = exp(bv0 + uv0)

###############################################################################

# sample size
n_ind <- 50
n_obs <- 10
n_tot <- n_ind * n_obs

# fixed effects parameters
## Phenotype
# mean aggressiveness
y_bm0 <- 0
# effect of temperature on aggressiveness
y_bmh <- 1
# mean variance in aggressiveness (log scale)
y_bv0 <- 0
# effect of sex on aggressiveness
y_bm1 <- 0.5
# effect of sex on variance in aggressiveness (log scale)
y_bv1 <- 0.5

## Oikotype
# mean temperature
h_bm0 <- 0
# mean variance in temperature (log scale)
h_bv0 <- 0


# random effects parameters
## individual variances
### mean aggressiveness
y_vm <- 0.3
### slope aggressiveness vs temperature
y_vs <- 0.1
### variance in aggressiveness
y_vv <- 0.2
### mean temperature
h_vm <- 0.3
### variance in temperature
h_vv <- 0.0001

### vector of individual standard deviations
sd_i <- diag(sqrt(c(y_vm, y_vs, y_vv, h_vm, h_vv)))

## correlation matrix
cor_m <- matrix(c(
  1, 0, -0.7, 0.7, 0,
  0, 1, 0, 0, 0,
  -0.7, 0, 1, -0.7, 0,
  0.7, 0, -0.7, 1, 0,
  0, 0, 0, 0, 1
), ncol = ncol(sd_i), byrow = TRUE)

## variance matrix
sigma <- sd_i %*% cor_m %*% sd_i

## individual values
ind <- rmvnorm(n_ind, rep(0, ncol(sigma)), sigma)

# putting everything together

dat <- data.frame(
  id = rep(seq(n_ind), each = n_obs),
  sex = rep(rbinom(n_ind, 1, 0.5), each = n_obs),
  y_um = rep(ind[, 1], each = n_obs),
  y_us = rep(ind[, 2], each = n_obs),
  y_uv = rep(ind[, 3], each = n_obs),
  h_um = rep(ind[, 4], each = n_obs),
  h_uv = rep(ind[, 5], each = n_obs)
)

dat <- dat %>%
  mutate(
    e_h = rnorm(n_tot, 0, exp(h_bv0 + h_uv)),
    h = h_bm0 + h_um + e_h,
    e_y = rnorm(n_tot, 0, exp(y_bv0 + y_uv + y_bv1 * sex)),
    y = y_bm0 + y_um + (y_bmh + y_us) * h + y_bm1 * sex + e_y
  ) %>%
  rename(
    temp = h,
    aggr = y
  ) %>%
  select(id, sex, aggr, temp)

write.csv(dat, "2_selkie_data.csv", row.names = FALSE)
