#### Libraries ####
library(tidyverse)
library(cmdstanr)

source("../utils/covid.R")
data_file <- "../data/data_preprocessed.csv"
data_file_rc_sep_reg <- "data_rc_separate_regions.csv"
stan_file <- "../models/sim_g.stan"
saved_models_folder <- "fitted-models/"


#### Data ####
df_pp <- read_csv(data_file)
df_pp_sep_reg <- read_csv(data_file_rc_sep_reg)

#### Model ####
main <- cmdstan_model(stan_file, cpp_options = list(stan_threads = T))
sample <- function(model = main, ...) {
  fit <- model$sample(
    data = stan_dat_omit(make_stan_dat(...)), # make_stan_dat generates the data that is provided to STAN
    seed = 12345,
    chains = 4,
    parallel_chains = 4,
    threads_per_chain = 2,
    refresh = 200,
    max_treedepth = 15,
    adapt_delta = 0.9
  )
}

sample_and_save <- function(model = main, ...) {
  # Sample
  rc_fit <- sample(model, ...)
  
  # Save fitted model
  par_args <- list(...) 
  arg_names <- names(par_args)
  par_values <- as.character(round(unlist(par_args),2))
  rc_name <- paste(arg_names, par_values, collapse = "_", sep = "_")
  save_file <- paste0(saved_models_folder, "rc_", rc_name, ".rds")
  
  rc_fit$save_object(file = save_file)
}


#### Start of modeling phase: ModelStartCases ####
sample_and_save(ModelStartCases = 10)
sample_and_save(ModelStartCases = 50)

#### End of modeling phase: EpidemicEnd ####
sample_and_save(EpidemicEnd = 3 * 7)
sample_and_save(EpidemicEnd = 5 * 7)

#### Spline: tb ####
sample_and_save(ta = -1, tb = 1)
sample_and_save(tb = 1)
sample_and_save(tb = 5)

#### Prior for theta_m: theta_mix ####
sample_and_save(theta_mix = 0.3)
sample_and_save(theta_mix = 0.5)

#### Infection to reporting: p_IN ####
par_p_in_est <- expand.grid(mu_mu0 = c(eval(formals(p_in)$mu) - 0.25, eval(formals(p_in)$mu), eval(formals(p_in)$mu) + 0.25),
                            sigma_beta0 = c(eval(formals(p_in)$sigma) - 0.2, eval(formals(p_in)$sigma), eval(formals(p_in)$sigma) + 0.2))
par_p_in_est$mu_sigma0 <- 0.5
par_p_in_est$sigma_alpha0 <- 2

for (i in 1:nrow(par_p_in_est)) {
  sample_and_save(p_in_mu_mu0 = par_p_in_est$mu_mu0[i], p_in_mu_sigma0 = par_p_in_est$mu_sigma0[i],
                  p_in_sigma_alpha0 = par_p_in_est$sigma_alpha0[i], p_in_sigma_beta0 = par_p_in_est$sigma_beta0[i])
}


#### Generation time: p_G ####
par_p_g <- expand.grid(par1 = c(eval(formals(p_g)$p1) - 0.5, eval(formals(p_g)$p1), eval(formals(p_g)$p1) + 0.5),
                       par2 = c(eval(formals(p_g)$p2) - 1, eval(formals(p_g)$p2), eval(formals(p_g)$p2) + 1))
for (i in 1:nrow(par_p_g)) {
  sample_and_save(par1_p_g = par_p_g$par1[i], par2_p_g = par_p_g$par2[i])
}


#### LOO Country ####
for (cid in 1:20) {
  df_pp_sub <- df_pp %>%
    dplyr::filter(country_id != cid)
  main_sub_fit <- main$sample(
    data = stan_dat_omit(make_stan_dat(dat = df_pp_sub)),
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    threads_per_chain = 2,
    refresh = 200,
    max_treedepth = 15,
    adapt_delta = 0.9
  )

  main_sub_fit$save_object(file = paste0(saved_models_folder, "rc_loo_", cid, ".rds"))
}

# ### Separate regions ####
rc_dat_sep_reg_looAUSA1 <- sample(dat = df_pp_sep_reg %>% dplyr::filter(country != "Western Australia"))
rc_dat_sep_reg_looAUSA1$save_object(paste0(saved_models_folder, "rc_separate_regions_loo_WA.rds"))
rc_dat_sep_reg_looAUSA2 <- sample(dat = df_pp_sep_reg %>% dplyr::filter(country != "Eastern Australia"))
rc_dat_sep_reg_looAUSA2$save_object(paste0(saved_models_folder, "rc_separate_regions_loo_EA.rds"))
