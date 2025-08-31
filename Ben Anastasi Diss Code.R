#-------------------------------------------------------------------------------
#Packages
#-------------------------------------------------------------------------------
library(purrr)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(BMA)
library(tidyverse)
library(glmnet)
library(brms)
library(loo)
library(future)
#library(cmdstanr)
library(future.apply)

#Scenarios
#Toxicity 0.01, 0.3; 0.02, 0.2; 0.02, 0.2; 0.02, 0.4; 0.05, 0.3
#Efficacy 10, 0.9; 20, 0.7; 20, 0.9; 10, 1.1; 30, 0.9;
#Sensitivity analysis for sd of emax model

#-------------------------------------------------------------------------------
#Config
#-------------------------------------------------------------------------------
n_datasets   <- 50
n_cohorts    <- 6
cohort_size  <- 3
eff_threshold <- 70
tox_threshold <- 0.25
probE <- 0.1
probT <- 0.1

dose_grid <- data.frame(dose = c(0.35, 0.7, 1.05, 1.4)) %>%
  mutate(dose2 = dose^2)

all_results <- data.frame()
all_responses <- data.frame()
all_weights <- data.frame()

combined_weights <- data.frame()
combined_responses <- data.frame()


#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Baseline - 0.02, 0.3
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904, b1 = 2.952, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)

#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
      ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                         data.frame(cohort = cohort_i,
                                    lin = weights[1],
                                    quad = weights[2],
                                    emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
    .groups = "drop"
  )

combined_weights <- bind_rows(combined_weights, avg_weights)

write.csv(all_results, "obd_baseline.csv", row.names = TRUE)
write.csv(combined_responses, "responses_baseline.csv", row.names = TRUE)
write.csv(combined_weights, "weights_baseline.csv", row.names = TRUE)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 2 - 0.01, 0.3
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -5.844 , b1 = 3.569, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
  
#safe in case of na; used for limit testing       
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd

#is.na in case of na obd found; used for limit testing    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 3 - 0.05, 0.3
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -3.643 , b1 = 1.997, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 4 - 0.02, 0.2
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.953 , b1 = 2.660, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 5 - 0.02, 0.4
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.851 , b1 = 3.188, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

write.csv(all_results, "obd_tox_varied.csv", row.names = TRUE)
write.csv(combined_responses, "responses_tox_varied.csv", row.names = TRUE)
write.csv(combined_weights, "weights_tox_varied.csv", row.names = TRUE)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 6 (now onto efficacy) - 10, 0.9
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904, b1 = 2.952, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 10,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
    .groups = "drop"
  )

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 7 - 20, 0.7
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.7,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 8 - 20, 1.1
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 1.1,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 9 - 30, 0.9
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
#Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
#Inverse-logit function
  pr = 1/(1+exp(-z))
  
#Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 30,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
#Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
#Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
#Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
#Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
#To store OBD
  obd_hist <- numeric(n_cohorts)
#Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
#Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
#Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
#Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
#Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Output results
#-------------------------------------------------------------------------------

write.csv(all_results, "obd_eff_varied.csv", row.names = TRUE)
write.csv(combined_responses, "responses_eff_varied.csv", row.names = TRUE)
write.csv(combined_weights, "weights_eff_varied.csv", row.names = TRUE)


#-------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 10 - 0.01, 0.2
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -5.844 , b1 = 3.166, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    #safe in case of na; used for limit testing       
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    #is.na in case of na obd found; used for limit testing    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 11 - 0.01, 0.4
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -5.792 , b1 = 3.936, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 12 - 0.05, 0.2
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -3.643 , b1 = 1.545, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 13 - 0.05, 0.40
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -3.588 , b1 = 2.215, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

write.csv(all_results, "obd_tox_varied_2.csv", row.names = TRUE)
write.csv(combined_responses, "responses_tox_varied_2.csv", row.names = TRUE)
write.csv(combined_weights, "weights_tox_varied_2.csv", row.names = TRUE)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 14 (now onto efficacy) - 10, 0.7
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904, b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 10,
                    ec50 = 0.7,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
    .groups = "drop"
  )

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 15 - 10, 1.1
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 10,
                    ec50 = 1.1,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 16 - 30, 0.7
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 30,
                    ec50 = 0.7,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 17 - 30, 1.1
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 30,
                    ec50 = 1.1,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
  )  %>% 
  ungroup()

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#Output results
#-------------------------------------------------------------------------------

write.csv(all_results, "obd_tot_varied_2.csv", row.names = TRUE)
write.csv(combined_responses, "responses_tot_varied_2.csv", row.names = TRUE)
write.csv(combined_weights, "weights_tot_varied_2.csv", row.names = TRUE)


# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 18 - changing the model in toxicity to quadratic negative
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

n_datasets   <- 50
n_cohorts    <- 6
cohort_size  <- 3
eff_threshold <- 70
tox_threshold <- 0.25
probE <- 0.1
probT <- 0.1

dose_grid <- data.frame(dose = c(0.35, 0.7, 1.05, 1.4)) %>%
  mutate(dose2 = dose^2)

all_results <- data.frame()
all_responses <- data.frame()
all_weights <- data.frame()

combined_weights <- data.frame()
combined_responses <- data.frame()

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Baseline - 0.02, 0.3
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.963, b1 = 3.720, b2 = -.5, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose + b2*dose^2
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
    .groups = "drop"
  )

combined_weights <- bind_rows(combined_weights, avg_weights)

write.csv(all_results, "obd_baseline.csv", row.names = TRUE)
write.csv(combined_responses, "responses_baseline.csv", row.names = TRUE)
write.csv(combined_weights, "weights_baseline.csv", row.names = TRUE)


# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 19 - changing the model in toxicity to quadratic positive
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Baseline - 0.02, 0.3
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -3.887, b1 = 1.623, b2 = 0.5, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose + b2*dose^2
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
    .groups = "drop"
  )

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 20 - changing the model in toxicity to exponential
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Baseline - 0.02, 0.3
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -5.535, b1 = 1.155, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*exp(dose)
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------

all_weights <- do.call(rbind, lapply(results, function(x) x$model_weights))

#Calculate average weights per cohort
avg_weights <- all_weights %>%
  group_by(cohort) %>%
  summarise(
    avg_lin = mean(lin),
    avg_quad = mean(quad),
    avg_emax = mean(emax),
    .groups = "drop"
  )

combined_weights <- bind_rows(combined_weights, avg_weights)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 21 - changing the model in efficacy to gamma = .5
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Baseline - 0.02, 0.3
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904, b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 0.5,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)

# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 22 - changing the model in efficacy to gamma = 3
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Baseline - 0.02, 0.3
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904, b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 3,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


combined_responses <- bind_rows(combined_responses, avg_responses)


# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 23 - changing the model in efficacy to gamma = 5
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Data simulation
#-------------------------------------------------------------------------------
#Baseline - 0.02, 0.3
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904, b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 5,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_hist <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  #Creates data frame to store weights from BMA
  weights_hist <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Model weights
    loo_lin   <- loo(mod_linear)
    loo_quad  <- loo(mod_quad)
    loo_emax  <- loo(mod_emax)
    
    weights <- loo_model_weights(list(loo_lin, loo_quad, loo_emax), method = "stacking")
    
    weights_hist <- rbind(weights_hist, 
                          data.frame(cohort = cohort_i,
                                     lin = weights[1],
                                     quad = weights[2],
                                     emax = weights[3]))
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_preds <- weights[1]*pred_lin +
      weights[2]*pred_quad +
      weights[3]*pred_emax
    
    eff_probs <- colMeans(eff_preds < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(eff_preds),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_hist[cohort_i] <- common_obd
    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_hist, cohort_dose_preds = cohort_dose_preds,  model_weights = weights_hist))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_dose_preds))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox),
    .groups = "drop"
  )


write.csv(all_results, "obd_differentmodelstest.csv", row.names = TRUE)
write.csv(combined_responses, "responses_differentmodelstest.csv", row.names = TRUE)
write.csv(combined_weights, "weights_differentmodelstest.csv", row.names = TRUE)

#-------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Comparing baseline to single model results
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
                                       
#-------------------------------------------------------------------------------
#Config
#-------------------------------------------------------------------------------
n_datasets   <- 50
n_cohorts    <- 6
cohort_size  <- 3
eff_threshold <- 70
tox_threshold <- 0.25
probE <- 0.1
probT <- 0.1

dose_grid <- data.frame(dose = c(0.35, 0.7, 1.05, 1.4)) %>%
  mutate(dose2 = dose^2)

all_results <- data.frame()
all_responses <- data.frame()

combined_responses <- data.frame()

#-------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 24 - linear
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_linear <- brm(
    efficacy ~ dose, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(76.97, 10), class = "Intercept"),
      prior(normal(-7.67, 5),   class = "b"),
      prior(student_t(3, 0, 7), class = "sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_path <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_linear    <- update(mod_linear,    newdata = eff_df, recompile = FALSE)

    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_lin     <- posterior_epred(mod_linear, newdata = dose_grid)

    eff_probs <- colMeans(pred_lin < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(pred_lin),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    #safe in case of na; used for limit testing       
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_path[cohort_i] <- common_obd
    
    #is.na in case of na obd found; used for limit testing    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_path, cohort_dose_preds = cohort_dose_preds))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 25 - quadratic
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_quad <- brm(
    efficacy ~ dose + dose2, data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(78.61, 10),  class = "Intercept"),
      prior(normal(-14.76, 6),  class = "b", coef = "dose"),
      prior(normal(5.06, 4),    class = "b", coef = "dose2"),
      prior(student_t(3,0,7), class="sigma")
    ), 
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  #To store OBD
  obd_path <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_quad      <- update(mod_quad,      newdata = eff_df, recompile = FALSE)
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_quad    <- posterior_epred(mod_quad, newdata = dose_grid)
    
    eff_probs <- colMeans(pred_quad < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(pred_quad),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    #safe in case of na; used for limit testing       
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_path[cohort_i] <- common_obd
    
    #is.na in case of na obd found; used for limit testing    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_path, cohort_dose_preds = cohort_dose_preds))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)

#-------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#SIMULATION 26 - emax
# -------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------
#Toxicity simulation
tox_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    b0 = -4.904 , b1 = 2.952, sigma = 1) {
  #Repeat for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Assumed true model (could change this to b0 + b1*dose and then fit more complex models)
  z = b0 + b1*dose
  
  #Inverse-logit function
  pr = 1/(1+exp(-z))
  
  #Bernoulli response variable
  tox = rbinom(length(pr),1,pr)
  
  return(data.frame(toxicity = tox, dose = dose))
}

#Efficacy simulation
eff_sim <- function(nrep = 3,
                    dose_levels = c(0.35, 0.7, 1.05, 1.4),
                    n = 1,
                    e0 = 79.2,
                    emax = 20,
                    ec50 = 0.9,
                    sd_efficacy = 7) {
  
  #Repeat dose levels for 3 patients per dose
  dose <- rep(dose_levels, each = nrep)
  
  #Calculate mean efficacy using Hill (sigmoidal form) model
  mean_efficacy <- e0 - (emax * dose^n) / (ec50^n + dose^n)
  
  #Simulate efficacy values for each dose level
  eff <- rnorm(length(dose), 
               mean = mean_efficacy, 
               sd = sd_efficacy)
  
  #Create data frame
  eff_df <- data.frame(
    dose = dose,
    efficacy = eff)
  return(data.frame(dose = dose, efficacy = eff))
}

#-------------------------------------------------------------------------------
#Simulate a single iteration
#-------------------------------------------------------------------------------
simulate_dataset <- function() {
  
  eff_df <- eff_sim(nrep = 12) %>% mutate(dose2 = dose^2)
  tox_df <- tox_sim(nrep = 12)
  
  tox_model_log <- brm(
    toxicity ~ exp(dose), data = tox_df,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(-5.891, 2.5), class = "Intercept"),
      prior(normal(1.308, 2.5), class = "b")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  mod_emax <- brm(
    bf(efficacy ~ e0 - (emax * dose) / (ec50 + dose),
       e0 + emax + ec50 ~ 1, nl = TRUE),
    data = eff_df, family = gaussian(),
    prior = c(
      prior(normal(79.2,20), nlpar="e0"),
      prior(normal(19.2,20), nlpar="emax"),
      prior(uniform(0.35, 1.4), nlpar = "ec50"),
      prior(student_t(3,0,7), class="sigma")
    ),
    chains = 4, cores = 1, iter = 1000, refresh = 0
  )
  
  
  #To store OBD
  obd_path <- numeric(n_cohorts)
  #Creates data frame to store dose_summary and average over later
  cohort_dose_preds <- data.frame()
  
  for (cohort_i in 1:n_cohorts) {
    
    #Update models
    tox_model_log <- update(tox_model_log, newdata = tox_df, recompile = FALSE)
    mod_emax      <- update(mod_emax,      newdata = eff_df, recompile = FALSE)
    
    #Predictions
    pred_log_tox <- posterior_epred(tox_model_log, newdata = dose_grid)
    pred_emax    <- posterior_epred(mod_emax, newdata = dose_grid)
    
    eff_probs <- colMeans(pred_emax < eff_threshold)
    tox_probs <- colMeans(pred_log_tox < tox_threshold)
    
    dose_summary <- data.frame(
      cohort = cohort_i,
      dose = dose_grid$dose,
      eff_probs, tox_probs,
      mean_pred_eff = colMeans(pred_emax),
      mean_pred_tox = colMeans(pred_log_tox)
    )
    cohort_dose_preds <- rbind(cohort_dose_preds, dose_summary)
    
    #safe in case of na; used for limit testing       
    safe_doses <- dose_summary %>% filter(eff_probs > probE, tox_probs > probT)
    if (nrow(safe_doses) == 0) {
      common_obd <- NA
    } else {
      common_obd <- safe_doses$dose[which.min(safe_doses$mean_pred_eff)]
    }
    
    obd_path[cohort_i] <- common_obd
    
    #is.na in case of na obd found; used for limit testing    
    if (!is.na(common_obd)) {
      new_eff <- eff_sim(nrep = cohort_size, dose_levels = common_obd) %>% mutate(dose2 = dose^2)
      new_tox <- tox_sim(nrep = cohort_size, dose_levels = common_obd)
      eff_df <- bind_rows(eff_df, new_eff)
      tox_df <- bind_rows(tox_df, new_tox)
    }
  }
  
  return(list(obd = obd_path, cohort_dose_preds = cohort_dose_preds))
}

#-------------------------------------------------------------------------------
#Parallelism
#-------------------------------------------------------------------------------
plan(multisession, workers = 8)
start <- Sys.time()
results <- future_lapply(1:n_datasets, function(i) simulate_dataset(), future.seed = TRUE)
end <- Sys.time()
print(end - start)

#-------------------------------------------------------------------------------
#Results summary
#-------------------------------------------------------------------------------
obd_df <- do.call(rbind, lapply(seq_along(results), function(i) {
  data.frame(
    dataset = i,
    cohort = seq_along(results[[i]]$obd),
    dose = results[[i]]$obd,
    eff_threshold = eff_threshold,
    tox_threshold = tox_threshold
  )
}))

obd_frequency <- obd_df %>%
  count(cohort, dose, name = "count")

all_results <- bind_rows(all_results, obd_frequency)

#-------------------------------------------------------------------------------

#Combine all response predictions per cohort
all_responses <- do.call(rbind, lapply(results, function(x) x$cohort_))

#Average per dose per cohort across datasets
avg_responses <- all_responses %>%
  group_by(cohort, dose) %>%
  summarise(
    avg_eff = mean(mean_pred_eff),
    avg_tox = mean(mean_pred_tox)
  ) %>%
  ungroup()

combined_responses <- bind_rows(combined_responses, avg_responses)


write.csv(all_results, "obd_singlemodel.csv", row.names = TRUE)
write.csv(combined_responses, "responses_singlemodel.csv", row.names = TRUE)
                                       




