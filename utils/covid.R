require(tidyverse)
require(lubridate)

# Subset data based on n_cases at start and n_days after more than 50% of regions have implemented the last measure
start_end <- function(dat, npi_vars = c("schools", "borders", "events", "gatherings", "venues", "lockdown", "homeoffice"), n_cases = 100, n_days = 3 * 7) {
  end_dates <- dat %>%
    dplyr::select(all_of(c("country", "date", npi_vars))) %>%
    reshape2::melt(id.vars = c("country", "date")) %>%
    group_by(country, variable) %>%
    dplyr::filter(value >= max(value)) %>%
    arrange(date) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(country) %>%
    arrange(desc(date)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(date = ifelse(country == "Australia", "2020-04-02", as.character(date)),
           date = ifelse(country == "Canada", "2020-04-01", as.character(date)),
           date = as.Date(ifelse(country == "Eastern Australia", "2020-04-02", as.character(date)))) %>%
    mutate(end_date = date %m+% days(n_days)) %>%
    dplyr::select(country, end_date)
  
  dat <- dat %>%
    left_join(end_dates, by = "country") %>%
    dplyr::filter(cases >= n_cases) %>%
    dplyr::filter(date <= end_date) 
  
  return(dat)
} 


# First order spline when having population-weighted NPI measures
weighted_spline <- function(x, t0, t1) {
  # Measure not implemented at all
  if (all(x == 0)) {
    return(y_sum <- rep(0, length(x)))
  }
  # Measures fully implemented from the beginning
  if (all(x == 1)) {
    return(y_sum <- rep(0, length(x)))
  }
  
  # Change points
  cp <- numeric(length(x))
  cp[1] <- 0
  k <- 0
  for (i in 1:(length(x)-1)) {
    if (x[i+1] > x[i]) {
      k <- k + 1
      cp[i+1] <- k
    }
  }
  
  # t for all change points, i.e., t by regional aggregates
  cp_uni <- unique(cp)[-1]
  cp_imp <- ifelse(cp == 0, NA, cp)
  cp_imp <- zoo::na.locf(cp_imp, fromLast = F, na.rm = F)
  cp_imp <- ifelse(is.na(cp_imp), 0, cp_imp)
  t <- list()
  for (k in 1:length(cp_uni)) {
    t[[k]] <- numeric(length(x))
    m <- 0
    for (i in 1:length(x)) {
      if (cp_imp[i] < cp_uni[k]) {
        t[[k]][i] <- 0
      } else {
        m <- m + 1
        t[[k]][i] <- m
      }
    }
  }
  
  # For negative t0, adjust t
  if (t0 < 0) {
    t <- lapply(t, function(t_k) {
      t_k_adj <- t_k + t0
      t_k_adj <- dplyr::lead(t_k_adj, abs(t0)) 
      t_k_adj[is.na(t_k_adj)] <- max(t_k_adj, na.rm = T) + 1:(abs(t0))
      return(t_k_adj)
      })
  }
  
  # Splines 
  y <- list()
  for (k in 1:length(t)) {
    y[[k]] <- numeric(length(t[[k]]))
    for (i in 1:length(t[[k]])) {
      y[[k]][i] <- ifelse(t[[k]][i] <= t0, 0, ifelse(t[[k]][i] < t1, (t[[k]][i] - t0) / (t1 - t0), 1))
    }
  }
  
  # Population share
  x_uni <- unique(x)
  dx <- diff(x_uni)
  
  # Weight splines
  y_w <- list()
  for (k in 1:length(y)) {
    y_w[[k]] <- y[[k]] * dx[k]
  }
  
  # Aggregate
  y_sum <- numeric(length(x))
  for (i in 1:length(x)) {
    y_sum[i] <- sum(sapply(y_w, "[[", i))
  }
  
  return(y_sum)  
}

# Create stan data list
make_stan_dat <- function(dat = df_pp, npis = c("schools", "borders", "events", "gatherings", "venues", "lockdown", "homeoffice"), 
                          beforeEpiStart = 7+1, ModelStartCases = 100, EpidemicEnd = 28,
                          p_in_mu_mu0 = eval(formals(p_in)[["mu"]]), p_in_mu_sigma0 = 0.5, p_in_sigma_alpha0 = 2, p_in_sigma_beta0 = 2 / eval(formals(p_in)[["sigma"]]),
                          mu_p_in = eval(formals(p_in)[["mu"]]), sigma_p_in = eval(formals(p_in)[["sigma"]]), 
                          distr_p_g = "weibull", par1_p_g = eval(formals(p_g)[["p1"]]), par2_p_g = eval(formals(p_g)[["p2"]]), 
                          ta = 0, tb = 3,
                          theta_mix = 0.1) {
  # Delays
  p_g_vec <- rev(vp(200, F, p_g, distr_p_g, par1_p_g, par2_p_g)[-1])
  p_in_vec <- rev(vp(200, T, p_in, mu_p_in, sigma_p_in))
  
  # EpidemicEnd 
  if (is.null(EpidemicEnd)) {
    n_days_c <- dat %>%
      group_by(country) %>%
      summarize(n = n()) %>%
      ungroup()
    EpidemicEnd <- max(n_days_c$n)
  }
  
  # Subset study period
  dat_studyPeriod = dat %>%
    start_end(n_cases = 100, npi_vars = npis, n_days = EpidemicEnd)
  
  # EpiStart
  EpiStart <- floor(log(ModelStartCases, 3.28) * 5.5 + exp(mu_p_in + sigma_p_in^2 / 2)) + 1
  
  # Define dimensions and initialize matrices
  studyPeriod <- dat_studyPeriod %>%
    group_by(country) %>%
    summarize(N = n())
  
  max_N <- max(studyPeriod$N)
  max_N1 <- max_N + EpiStart
  max_N2 <- max_N + EpiStart + beforeEpiStart
  M <- length(npis)
  J <- length(unique(dat$country_id))
  N <- numeric(length = J)
  N1 <- numeric(length = J)
  cum_N1 <- numeric(length = J+1)
  cum_N1[[1]] <- 0
  dates <- list()
  Measures <- list()
  new_cases <- list()
  
  k <- 0
  for (j in unique(dat$country_id)) {
    k <- k + 1
    # Subset country
    dat_studyPeriod_j <- dat_studyPeriod %>%  
      dplyr::filter(country_id == j) 
    
    # Subset full data
    start <- dat_studyPeriod_j$date[1] %m-% days(EpiStart)
    end <- dat_studyPeriod_j$date[nrow(dat_studyPeriod_j)]
    dat_j <- dat %>% 
      dplyr::filter(country_id == j) %>%
      dplyr::filter(date >= start & date <= end) 
    Measures_j <- dplyr::select(dat_j, all_of(npis)) 
    new_cases_j <- dat_j$new_cases
    
    # Padding 
    pad_N <- (nrow(dat_studyPeriod_j) + EpiStart) - nrow(dat_j)
    if (pad_N > 0) {
      Measures_j <- bind_rows(
        data.frame(matrix(0, nrow = pad_N, ncol = M)) %>%
          set_names(npis),
        Measures_j
      )
      new_cases_j <- c(rep(0, pad_N), new_cases_j)
    }
    
    # Spline measures
    if (tb != 0) {
      Measures_j <- Measures_j %>% apply(2, weighted_spline, ta, tb)
    } 
    
    # Assign
    dates[[k]] <- seq(start, end, by = "day")
    N[k] <- nrow(dat_studyPeriod_j)
    N1[k] <- N[k] + EpiStart
    Measures[[k]] <- Measures_j
    new_cases[[k]] <- new_cases_j
  } 
  
  # Reduce
  cum_N1[2:(J+1)] <- cumsum(N1)
  Measures <- do.call(rbind, Measures)
  new_cases <- Reduce(c, new_cases)
  
  # Prepare stan data list
  stan_dat <- list(
    date = dates, J = J, M = M,
    EpiStart = EpiStart, beforeEpiStart = beforeEpiStart, N = N, N1 = N1, max_N1 = max_N1, max_N2 = max_N2, cum_N1 = cum_N1,
    new_cases = new_cases, X = Measures, 
    theta_mix = theta_mix, p_in = tail(p_in_vec, max_N2), p_g = tail(p_g_vec, max_N2),
    p_in_mu_mu0 = p_in_mu_mu0, p_in_mu_sigma0 = p_in_mu_sigma0, p_in_sigma_alpha0 = p_in_sigma_alpha0, p_in_sigma_beta0 = p_in_sigma_beta0
  )
  
  return(stan_dat)
}


# Omit variable from stan data list (e.g., there is an error when having dates)
stan_dat_omit <- function(stan_dat, var = "date") {
  stan_dat[[var]] <- NULL
  return(stan_dat)
}

# Get sigma for right shift of lognormal
get_sigma_p <- function(x, base_x, base_mu, base_sigma) {
  p_x <- dlnorm(base_x, base_mu, base_sigma)
  return( 1 / (x * p_x * sqrt(2 * pi)) )
}

# Get sigma for robustness check
rc_get_sigma_p <- function(mu, base_sigma) {
  base_mu <- log(mu[3])
  sigma <- numeric(length(mu))
  sigma[1] <- get_sigma_p(mu[1], mu[3], base_mu, base_sigma)
  sigma[2] <- sigma[1]
  sigma[3] <- base_sigma
  sigma[5] <- get_sigma_p(mu[5], mu[3], base_mu, base_sigma)
  sigma[4] <- sigma[5]
  return(sigma)
}

# Delays
# p_IN = Incubation period + symptom onset to reporting delay
# 1 From McAloon: Incubation period ~ Lognormal(1.63, 0.5)
# 2 From Cereda: Reporting delay ~ Gamma(1.88, 1/0.26)
mean_1 <- exp(1.63 + 0.5^2 / 2)
var_1 <- ( exp(0.5^2) - 1 ) * exp(2 * 1.63 + 0.5^2)
mean_2 <- 1.88 / 0.26
var_2 = 1.88 / 0.26^2
mean_pIN <- mean_1 + mean_2
var_pIN <- var_1 + var_2
mu_pIN <- log(mean_pIN^2 / sqrt(mean_pIN^2 + var_pIN))
sigma_pIN <- sqrt( log(1 + var_pIN / mean_pIN^2) )
# --> Lognormal(2.47, 0.45)
p_in <- function(x, mu = mu_pIN, sigma = sigma_pIN) { plnorm(x, mu, sigma) }

# From Ferretti: Weibull(3.2862, 6.1244)
p_g <- function(x, distr = "weibull", p1 = 3.2862, p2 = 6.1244) { 
  if (distr == "weibull") { pweibull(x, p1, p2) }
}

# Vectorize delay function
vp <- function(xT, from0, FUN, ...) {
  x <- seq(0,xT)
  y <- length(x)
  if (from0) {
    y[1] <- FUN(0.5, ...)
    y[2] <- FUN(1.5, ...) - FUN(0.5, ...)
  } else {
    y[1] <- 0
    y[2] <- FUN(1.5, ...)
  }
  for (i in 3:length(x)) {
    y[i] = FUN(x[i] + 0.5, ...) - FUN(x[i] - 0.5, ...)
  }
  return(y)
}


# Simulate data for delay
data_delay <- function(xT, from0, FUN, ...) {
  fun_args <- list(...)
  if (length(fun_args) == 0) {
    fun_args <- formals(FUN)[-1]
    arg_names <- names(formals(FUN)[-1])
    dat <- data.frame(x = 0:xT, value = vp(xT, from0, FUN, ...))
    return(dat)
  } else {
    l <- length(fun_args[[1]])
    p_x <- mapply(function(...) invoke(vp, xT = xT, from0 = from0, FUN = FUN,  ...), ...)
    p_par <- mapply(function(...) paste0("(", paste(sapply(list(...), round, 2), collapse = ", "), ")"), ...)
    dat <- data.frame(cbind(0:xT, p_x))
    colnames(dat) <- c("x", p_par)
    dat_long <- reshape2::melt(dat, "x") %>%
      mutate(par1 = stringi::stri_extract(variable, regex = "\\d.\\d{1,2}"),
             par2 = gsub(", ", "", stringi::stri_extract(variable, regex = ", \\d.\\d{1,2}")))
    return(dat_long)
  }
}


# Default plot for robustness checks
plot_rc <- function(..., color_levs, color_lab, default, xlab = "Estimated reduction in new infections (%)", xl = c(NA, .7), brks = seq(0,1,.25)) {
  dat <- do.call(rbind, ...) %>% 
    ungroup() %>%
    mutate(spec = rep(color_levs, each = 7 * 4000)) %>%
    mutate(spec = factor(spec, levels = color_levs)) %>% 
    mutate(def = rep(default, each = 7 * 4000)) %>%
    mutate(measure = as.factor(recode(as.character(measure), !!! npi_names))) %>%
    mutate(measure_ranking = measure) %>%
    mutate(measure_ranking = recode(measure_ranking, !!! npi_ranking))
  rc_pl <- ggplot(dat, aes(y = reorder(measure, measure_ranking), x = .value, color = spec, linetype = def)) +
    stat_pointintervalh(point_interval = ci_interval, 
                        position =  ggstance::position_dodgev(height = .75), 
                        shape = 21, point_size = 2, fill = "white") +
    geom_vline(aes(xintercept = 0), color = "red", linetype = "dotted") +
    scale_x_continuous(expand = c(0,0), limits = xl, labels = function(x) x * 100, breaks = brks)  +
    scale_color_viridis_d() +
    scale_linetype(guide = F) +
    guides(colour = guide_legend(override.aes = list(linetype = ifelse(default, "dashed", "solid")))) +
    labs(x = xlab,
         color = color_lab,
         linetype = color_lab) +
    theme_nature()  +
    theme(axis.title.y = element_blank(), legend.position = "bottom", legend.key.size =  unit(1, "cm")) 
}


# Plot model fit with compartments
plot_compartments <- function(ctry, stanDat, E_N, E_I, npi_labs = c("S", "B", "E", "G", "V", "H", "W"), pl_type = "pdf") {
  # Select country id
  j <- cid_country %>% dplyr::filter(country == ctry) %>% dplyr::select(country_id) %>% unlist
  date_j <- as.Date(stanDat$date[[j]])
  
  # Observed new cases
  observed_cases <- data.frame(
    .value = stanDat$new_cases[(stanDat$cum_N1[j]+1):stanDat$cum_N1[j+1]],
    .lower = NA,
    .upper = NA
    )
  
  # Expected new cases
  new_cases <- E_N %>%
    dplyr::filter(country == j) %>%
    dplyr::filter(day <= stanDat$N1[[j]]) %>% 
    ci_interval() %>%
    ungroup() %>%
    dplyr::select(.value, .lower, .upper)
  
  # Expected new infections
  new_infections <- E_I %>%
    dplyr::filter(country == j) %>%
    dplyr::filter(day <= stanDat$N1[[j]]) %>%
    ci_interval() %>%
    ungroup() %>%
    dplyr::select(.value, .lower, .upper)
  
  # Combine compartments
  if (pl_type == "tex") {
    comp_names <- c("$\\mu^I$", "$\\mu^N$", "Observed N")
  } else {
    comp_names <- c("EI", "EN", "N")
  }
  compartments <- rbind(new_infections, new_cases, observed_cases) %>%
    mutate(date = rep(date_j, times = 3),
           variable = factor(rep(comp_names, each = nrow(new_cases)), levels = comp_names)) %>%
    dplyr::select(date, variable, .value, .lower, .upper) 
  
  # Plot
  upper_y <- max(c(compartments$.upper, compartments$.value), na.rm = T)
  pl <- ggplot() +
    geom_line(data = compartments, mapping = aes(x = date, y = .value, color = variable)) +
    geom_ribbon(data = compartments, mapping = aes(x = date, ymin = .lower, ymax = .upper, fill = variable)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, upper_y*1.25)) +
    scale_x_date(expand = c(0,0), breaks = seq.Date(min(compartments$date), max(compartments$date, na.rm = T), by = "2 week"), 
                 date_labels ="%b %d") +
    labs(title = ctry, y = "Number of new cases/infections", color = "", fill = "")
  if (pl_type == "tex") {
    pl <- pl +
      scale_color_manual(values = alpha(c("seagreen3", "deepskyblue4", "black"), 1)) +
      scale_fill_manual(values = alpha(c("seagreen3", "deepskyblue4", "black"), .2)) 
  } else {
    pl <- pl + 
      scale_color_manual(values = alpha(c("EI" = "seagreen3", "EN" = "deepskyblue4",  "N" = "black"), 1), labels = expression(mu^I,mu^N,"Observed "*N)) +
      scale_fill_manual(values = alpha(c("EI" = "seagreen3", "EN" = "deepskyblue4",  "N" = "black"), .2), labels = expression(mu^I,mu^N,"Observed "*N)) 
  }
    pl <- pl +
    theme_nature() +
    theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())
    #guides(color = guide_legend(nrow = 2, byrow = T), fill = guide_legend(nrow = 2, byrow = T))
  
  # Add timing of NPIs
  Measures <- data.frame(stanDat$X[(stanDat$cum_N1[j]+1):stanDat$cum_N1[j+1], ])
  first_idx <- Measures %>% 
    summarize_all(function(x) date_j[which(x > 0)[1]]) %>% 
    # NPI is w.r.t. Ctilde, thus the drop occurs from day t-1 to t --> move date one day earlier, so that drop corresponds to intro of NPI
    mutate_all(function(x) x  %m-% days(1)) %>%
    gather %>%
    mutate(id = npi_labs) %>%
    group_by(value) %>%
    summarize(id = paste(id, collapse = "\n"),
              n = n()-1) %>%
    ungroup() %>%
    mutate()
  pl <- pl +
    geom_segment(data = first_idx, mapping = aes(x = value, xend = value, y = 0, yend = upper_y*(1.05-n*0.04)), 
                 linetype = "dashed", color = "red", alpha = .5) +
    geom_text(data = first_idx, mapping = aes(x = value, y = upper_y*(1.1-n*0.02), label = id), 
              size = 6 / cm(1), color = "red")
    
  
  # Add seeding and modeling phase
  seed_phase <- data.frame(xstart = date_j[1], yend = upper_y) %>%
    mutate(x = xstart %m+% days(floor(stanDat$EpiStart / 2)),
           xend = xstart %m+% days(stanDat$EpiStart))
  modeling_phase <- data.frame(xstart = date_j[1] %m+% days(stanDat$EpiStart), yend = upper_y) %>%
    mutate(xend = date_j[length(date_j)]) %>%
    mutate(x = seq.Date(from = xstart, to = xend, "days") %>% .[[round(length(.) / 2)]])
  pl <- pl +
    geom_segment(data = seed_phase, mapping = aes(x = xstart, xend = xend, y = yend*1.15, yend = yend*1.15), 
                 color = "grey50", arrow = arrow(ends = "both", angle = 10, length = unit(.25, "cm"), type = "closed")) +
    geom_text(data = seed_phase, mapping = aes(x = x, y = yend*1.2, label = ifelse(pl_type == "pdf", "Non-modeling phase", "Non-modeling")), 
              size = 8 / cm(1), color = "grey50") +
    geom_segment(data = modeling_phase, mapping = aes(x = xstart, xend = xend, y = yend*1.15, yend = yend*1.15), 
                 color = "grey50", arrow = arrow(ends = "both", angle = 10, length = unit(.25, "cm"), type = "closed")) +
    geom_text(data = modeling_phase, mapping = aes(x = x, y = yend*1.2, label = "Modeling phase"), 
              size = 8 / cm(1), color = "grey50")
  
  return(pl)
} 
