prep_file <- "utils/covid.R"
data_file <- "data/data_preprocessed.csv"
stan_file <- "models/sim_g.stan"
save_file <- "fitted-models/base.rds"

library(tidyverse)
source(prep_file)

df_pp <- read_csv(data_file) 

stanDat_main <- stan_dat_omit(make_stan_dat(df_pp))

# Fit model
main <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = T))
main_fit <- main$sample(
  data =  stanDat_main,
  seed = 12345,
  chains = 4,
  iter_warmup = 1000, 
  iter_sampling = 1000,
  parallel_chains = 4,
  threads_per_chain = 2, 
  refresh = 200,
  max_treedepth = 15,
  adapt_delta = 0.9
)
main_fit$save_object(file = save_file)

