require(posterior)
require(cmdstanr)
require(tidybayes)
require(bayesplot)
require(LaplacesDemon)
require(loo)

# To support tidybayes gather_draws for cmdstanr
tidy_draws.CmdStanMCMC <- function(model, ...) { return(as_draws_df(model$draws())) }

# LOO for ragged matrices (padded matrices that are originally of unequal length0)
loo.ragged_mat <- function(fit) {
  fit$draws("log_lik") %>% 
    apply(3, na.omit) %>%
    Filter(f = length) %>%
    abind::abind(along = 3) %>%
    loo
}