functions {
  // Discretization of the cumulative lognormal distribution function
  real diff_lnorm(real x, real mu, real sigma) {
    if (x == 0) {
      return lognormal_cdf(0.5, mu, sigma);
    } else {
      return lognormal_cdf(x+0.5, mu, sigma) - lognormal_cdf(x-0.5, mu, sigma);
    }
  }
  // Mixture distribution of a Halfnormal for negative and Uniform distributioon for positive values
  real mixtureHalfNormalUniform_lpdf(real param, real w, real a, real b) {
    if ( param < 0 ) {
      return log(w) + normal_lpdf(param|0, w / (sqrt(2 * pi()) * (1 - w)));
    } else {
      return log1m(w) + uniform_lpdf(param|a, b);
    }
  }
}


data {
  int<lower=0> max_N1; // maximum number of observations (days in seeding + modelling phase) across countries
  int<lower=0> max_N2; // maximum number of observations (days before seeding phase + in seeding phase + modeling phase) across countries
  int<lower=1> J; // number of countries
  int<lower=1> EpiStart; // -EpiStart is the start of the seeding phase (default is 33 days prior to 100 cumulative cases)
  int<lower=1> beforeEpiStart; // 1 + the number of days before the seeding phase (default is 1+7)
  int<lower=1> M; // number of NPIs
  int<lower=1> N[J]; // number of days in the modeling phase for each country j
  int<lower=1> N1[J]; // number of days in the seeding + modeling phase for each country j
  int<lower=0> cum_N1[J+1]; // observations are concatenated to one vector, this is the index vector for the time series (including seeding phase)
  int<lower=0> new_cases[sum(N1)]; // documented number of new cases
  matrix[sum(N1),M] X; // feature matrix adjusted for population weights and time-delayed response function (sum_Rj p_rj * f(T_mjrt))
  real p_in_mu_mu0; // prior: location hyperparameter mu0 in mu^p_IN ~ Normal(mu0, sigma0)
  real p_in_mu_sigma0; // prior: scale hyperparameter sigma0 in mu^p_IN ~ Normal(mu0, sigma0)
  real p_in_sigma_alpha0; // prior: shape hyperparameter alpha0 in sigma^p_IN ~ Gamma(alpha0, beta0)
  real p_in_sigma_beta0; /// prior: inverse scale hyperparameter beta0 in sigma^p_IN ~ Gamma(alpha0, beta0)
  vector<lower=0,upper=1>[max_N2] p_g; // assumed generation time distribution in reversed order (rev(p_G))
  real<lower=0,upper=1> theta_mix; // weight for negative effects in the mixture prior for theta_m
}


parameters {
  real alpha_0; // exp(alpha_0) is the global transmission rate in the absence of NPIs 
  vector[J] alpha_j; // exp(alpha_0 + alpha_j) is the country-specific transmission rate in the absence of NPIs
  real<lower=0> tau; // between-country variation for alpha_j
  vector<upper=1>[M] theta; // NPI effects
  real<lower=0> xi_new_cases; // over-dispersion on the parameter 1 / sqrt(phi^N) for the number of new cases
  real<lower=0> xi_new_infections; // over-dispersion on the parameter 1 / sqrt(phi^I) for the number of new infections
  real<lower=0> new_infections_EpiStart[J]; // number of new infections in country j at the start of the seeding phase, i.e., at day t = -EpiStart
  real<lower=0> lambda; // expected number of new infections in country j at the start of the seeding phase, i.e., at day t = -EpiStart
  vector<lower=0>[sum(N1)] new_infections; // number of new infections (latent / unobserved)
  real mu_p_in; // log mean in p_IN ~ Lognormal(mu, sigma)
  real<lower=0> sigma_p_in; // log standard deviation in p_IN ~ Lognormal(mu, sigma)
}


transformed parameters {
  real<lower=0> phi_new_cases = (1. / xi_new_cases) ^ 2; // over-dispersion parameter for the number of new cases
  real<lower=0> phi_new_infections = (1. / xi_new_infections) ^ 2; // over-dispersion parameter for the number of new infections
  vector<lower=0>[sum(N1)] mu_new_cases; // expected number of new cases
  vector<lower=0>[sum(N1)] being_contagious; // number of contagious subjects (up to a normalizing constant)
  vector<lower=0>[sum(N1)] mu_new_infections; // expected number of new infections
  vector<lower=0>[sum(N1)] sigma2_new_infections; // variance of the expected number of new infections
  vector[sum(N1)] log_mu_new_infections;
  vector<lower=0>[max_N2] p_in; // probability distribution for the time from infection to reporting of a new case p_IN
  matrix<lower=0>[beforeEpiStart,J] initially_infected; // number of initially infected before the start of the seeding phase (including new_infections_EpiStart)


  // Compute discretized p_IN distribution
  for (s in 1:max_N2) { 
   p_in[s] = diff_lnorm(max_N2-s, mu_p_in, sigma_p_in);
  }

  for (j in 1:J) {
    // number of new infections at t = -Epistart
    initially_infected[beforeEpiStart,j] = new_infections_EpiStart[j];
    
    // initially infected before EpiStart
    for (i in 1:(beforeEpiStart-1)) {
      initially_infected[i,j] = initially_infected[beforeEpiStart,j] / exp(3.28 / 5.5 * (beforeEpiStart-i));
    }
    
    // number of contagious subjects, expected number of new infections and cases at day t=-EpiStart+1
    being_contagious[cum_N1[j]+1] = dot_product(initially_infected[1:beforeEpiStart,j], tail(p_g, beforeEpiStart));
    mu_new_infections[cum_N1[j]+1] = being_contagious[cum_N1[j]+1] * exp(alpha_0 + alpha_j[j]) * prod(1 - to_vector(X[cum_N1[j]+1]) .* theta); // estimate number of new infections at EpiStart
    sigma2_new_infections[cum_N1[j]+1] = mu_new_infections[cum_N1[j]+1] * (1 + mu_new_infections[cum_N1[j]+1] / phi_new_infections);
    mu_new_cases[cum_N1[j]+1] = dot_product(append_row(initially_infected[1:beforeEpiStart,j], mu_new_infections[cum_N1[j]+1]),  tail(p_in, beforeEpiStart+1));
    
    // number of contagious subjects (up to a normalizing constant), expected number of new infections and cases from day t=-EpiStart+2 onwards
    for (i in (cum_N1[j]+2):cum_N1[j+1]) {
      being_contagious[i] = dot_product(append_row(initially_infected[1:beforeEpiStart,j], new_infections[(cum_N1[j]+1):(i-1)]), tail(p_g, i-cum_N1[j]+beforeEpiStart-1));
      mu_new_infections[i] = being_contagious[i] * exp(alpha_0 + alpha_j[j]) * prod(1 - to_vector(X[i]) .* theta);
      sigma2_new_infections[i] = mu_new_infections[i] * (1 + mu_new_infections[i] / phi_new_infections);
      mu_new_cases[i] = dot_product(append_row(initially_infected[1:beforeEpiStart,j], new_infections[(cum_N1[j]+1):(i)]), tail(p_in, i-cum_N1[j]+beforeEpiStart));
    }
  }
}


model {
  // priors
  alpha_0 ~ student_t(7., 0., 10.0);
  alpha_j ~ normal(0., tau);
  tau ~ student_t(4., 0., 1.);
  for (m in 1:M) {
   theta[m] ~ mixtureHalfNormalUniform(theta_mix, 0., 1.);
  }
  xi_new_cases ~ normal(0., 1.);
  xi_new_infections ~ normal(0., 1.);
  new_infections_EpiStart ~ exponential(1. / lambda);
  lambda ~ exponential(1.);
  mu_p_in ~ normal(p_in_mu_mu0, p_in_mu_sigma0);
  sigma_p_in ~ gamma(p_in_sigma_alpha0, p_in_sigma_beta0);
  
  // likelihood
  for (j in 1:J) {
    new_infections[(cum_N1[j]+1):cum_N1[j+1]] ~ normal(mu_new_infections[(cum_N1[j]+1):cum_N1[j+1]], sqrt(sigma2_new_infections[(cum_N1[j]+1):cum_N1[j+1]]));
    new_cases[(cum_N1[j]+1+EpiStart):cum_N1[j+1]] ~ neg_binomial_2(mu_new_cases[(cum_N1[j]+1+EpiStart):cum_N1[j+1]], phi_new_cases);
  }
}

generated quantities {
  // transform the expected number of contagious subjects (up to a normalizing constant), number of new infections and cases into matrix form
  matrix<lower=0>[max_N1,J] E_new_cases_rep = rep_matrix(positive_infinity(), max_N1, J);
  matrix<lower=0>[max_N1,J] E_new_infections_rep = rep_matrix(positive_infinity(), max_N1, J);
  matrix<lower=0>[max_N1,J] E_being_contagious_rep = rep_matrix(positive_infinity(), max_N1, J);
  
  // log-likelihood
  matrix[max_N1,J] log_lik;
  for (j in 1:J) {
    for (i in 1:N[j]) {
      log_lik[i,j] = neg_binomial_2_log_lpmf(new_cases[cum_N1[j]+i+EpiStart] | log(mu_new_cases[cum_N1[j]+i+EpiStart]), phi_new_cases);
    }
    
    // compartment quantities
    for (i in 1:N1[j]) {
      E_new_cases_rep[i,j] = mu_new_cases[cum_N1[j]+i];
      E_new_infections_rep[i,j] = mu_new_infections[cum_N1[j]+i];
      E_being_contagious_rep[i,j] = being_contagious[cum_N1[j]+i];
    }
  }
}
