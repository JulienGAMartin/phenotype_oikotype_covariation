################################################################################
# Field work (simulating data) for selkie case study
# Author: Julien Martin
# Date: 18 November 2024
###############################################################################

# loading packages
library(brms)
library(mvtnorm)
library(posterior)
library(bayestestR)
library(gt)
library(tidyverse)

# reading the data
dat <- read.csv("2_selkie_data.csv")

###############################################################################
# model fitting
###############################################################################

y_cs <- bf(
  aggr ~ temp + sex + (1 + temp | a | id),
  sigma ~ 1 + sex + (1 | a | id)
)
h_cs <- bf(
  temp ~ 1 + (1 | a | id),
  sigma ~ 1 + (1 | a | id)
)

m_cs <- brm(y_cs + h_cs + set_rescor(FALSE),
  data = dat,
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000, thin = 4
)

saveRDS(m_cs, file = "m_cs.rds")

if (!exists("m_cs")) m_cs <- readRDS("m_cs.rds")

summary(m_cs)

post_sp <- as_draws_df(m_cs, variable = "^b_|^sd_|^sigma_|^cor_", regex = TRUE) %>%
  rename_with(~ str_replace(.x, "^sd_", "var_")) %>%
  mutate(across(starts_with("var_"), ~ .x^2, .names = "{.col}"))


param_summary <- summarise_draws(post_sp, mean, median, sd, mad, default_convergence_measures())

# Calculate 95% HPDI for all parameters
hpdi_95 <- hdi(post_sp, ci = 0.95)

param_out <- full_join(param_summary, hpdi_95, by = join_by(variable == Parameter)) %>%
  select(variable, median, CI_low, CI_high, rhat) %>%
  mutate(trait = case_when(
    str_detect(variable, "_sigma_y_") ~ "agr_sigma",
    str_detect(variable, "_sigma_h_") ~ "temp_sigma",
    str_detect(variable, "_y_") ~ "agr",
    str_detect(variable, "_h_") ~ "temp",
    TRUE ~ "other"
  )) %>%
  mutate(param = case_when(
    str_detect(variable, "var_") ~ "var",
    str_detect(variable, "cor_") ~ "vcor",
    TRUE ~ "fix"
  )) %>%
  arrange(param, trait) %>%
  mutate(
    Parameter = str_remove(
      variable,
      "b_y_|b_h_|b_sigma_y_|b_sigma_h_"
    ),
    Parameter = str_replace_all(Parameter, "h", "Temp"),
    Parameter = str_replace(Parameter, "sex", "Sex"),
    Parameter = str_replace_all(Parameter, "y", "Aggr"),
    Parameter = str_replace(Parameter, "var_id__", "V("),
    Parameter = str_replace(Parameter, "cor_id__", "cor("),
    Parameter = str_replace_all(Parameter, "sigma_(Aggr|Temp)_Intercept", "\\1 dispersion"),
    Parameter = str_replace_all(Parameter, "(Aggr|Temp)_Intercept", "\\1 mean"),
    Parameter = str_replace(Parameter, "Aggr_Temp", "Aggr slope Temp"),
    Parameter = str_replace(Parameter, "__", ", "),
    Parameter = if_else(param != "fix", paste0(Parameter, ")"), Parameter),
    model = rep(c("Aggressiveness (mean)", "Aggressiveness (dispersion)", "Temperature (mean)", "Temperature (dispersion)", "Variance components"), c(3, 2, 1, 1, 15))
  ) %>%
  select(trait, param, model, variable, Parameter, median, CI_low, CI_high, rhat) %>%
  as.data.frame() %>%
  print()

# Format and display param_out as a gt table
gt_table <- param_out %>%
  mutate(
    median = round(median, 3),
    CI_low = round(CI_low, 3),
    CI_high = round(CI_high, 3),
    rhat = round(rhat, 3),
    HPDI = paste0(CI_low, "/", CI_high)
  ) %>%
  select(model, Parameter, median, HPDI, rhat) %>%
  #    arrange(model, Parameter) %>%
  gt(groupname_col = "model") %>%
  tab_header(title = "Selkie aggressiveness and temperature model") %>%
  cols_label(
    model = "Parameter type",
    Parameter = "Parameter",
    median = "Median",
    HPDI = "HPDI (95%)",
    rhat = "Rhat"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  )

print(gt_table)

gtsave(gt_table, "param_summary_table.docx")
