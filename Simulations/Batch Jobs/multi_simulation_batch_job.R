# get job id and set seed
jobid = as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(jobid)

# change-point packages
library(InspectChangepoint)
library(L2hdchange)
library(ecp)

# load MICH functions
source("~/MICH/mich.R")
source("~/MICH/mich_vector.R")
source("~/MICH/mich_matrix.R")
source("~/MICH/credible_sets.R")
source("~/MICH/priors.R")
source("~/MICH/utilities.R")
source("~/MICH/simulation_functions.R")
Rcpp::sourceCpp("~/MICH/mich.cpp")
Rcpp::sourceCpp("~/MICH/multi_mich.cpp")
Rcpp::sourceCpp("~/MICH/force_adj.cpp")

# significance level
level <- 0.9
delta <- 0.5
tol <- 1e-5

# simulation settings
T <- 250
min_space <- 10
settings <- expand.grid(L = c(5, 10, 20), d = c(10, 50, 100), p = c(0.1, 0.5, 1))
max_length <- log(T)^1.5 # detection threshold

# conditional detection window
window <- floor(min(sqrt(T) / 2, 15)) 

cols <- c("method", "L", "d", "p", "rho", "adapt", "time","L_est", "n_detected", 
          "n_covered", "avg_len", "hausdorff_1", 
          "hausdorff_2", "fpsle", "fnsle")

# initialize results matrix
result <- matrix(ncol = length(cols), nrow = 0)
result <- data.frame(result)
names(result) <- cols

for (j in 1:nrow(settings)) {
  print(j)
  # extract settings
  d <- settings$d[j]; L <- settings$L[j]; p <- settings$p[j];
    
  #### Identity #### 
  C <- 10
  rho <- 0
  adapt <- FALSE
  # generate data
  cp_data <- multi_simulation(T, L, d, p, rho, C, min_space, adapt)
  true_cp <- cp_data$changepoints
  
  #### ECP ####
  time <- system.time({fit <- e.divisive(cp_data$Y, sig.lvl = 1 - level, 
                                         min.size = min_space, alpha = 2,
                                         R=499)})[3]
  est_cp <- fit$estimates
  # estimate number of change-points
  L_est <- length(est_cp) - 2
  
  # store results
  result[nrow(result) + 1,"method"] <- "ecp"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), est_cp), 
                               hausdorff(est_cp, c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), est_cp), 
                               fnsle(c(1, true_cp, T+1), est_cp))
  
  #### Inspect ####
  penalty <- log(log(T)*d/2)
  while(TRUE) {
    time <- try(system.time({fit <- inspect(t(cp_data$Y), lambda = penalty)})[3])
    if(!"try-error" %in% class(time)) break
    penalty <- penalty + 1
  }
  
  est_cp <- fit$changepoints[,1] + 1
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # store results
  result[nrow(result) + 1,"method"] <- "Inspect"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### L2HDC ####
  time <- system.time({ts_l2_fit <- ts_hdchange(t(cp_data$Y)); fit <- hdchange(ts_l2_fit);})[3]
  est_cp <- fit$time_stamps + 1

  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # store results
  result[nrow(result) + 1,"method"] <- "L2HDC"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  #### MICH ####
  #### MICH Auto ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$Y, L_auto = TRUE, max_iter = Inf, tol = tol, restart = FALSE)})[3] 
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Auto Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$Y, L_auto = TRUE, max_iter = Inf, tol = tol, reverse = TRUE, restart = FALSE)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Auto rev"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Oracle ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$Y, L = L, max_iter = Inf, tol = tol)})[3] 
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Ora Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$Y, L = L, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora rev"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### Adaptive #### 
  C <- 25
  rho <- 0
  adapt <- TRUE
  # generate data
  cp_data <- multi_simulation(T, L, d, p, rho, C, min_space, adapt)
  true_cp <- cp_data$changepoints
  
  #### ECP ####
  time <- system.time({fit <- e.divisive(cp_data$Y, sig.lvl = 1 - level, 
                                         min.size = min_space, alpha = 2,
                                         R=499)})[3]
  est_cp <- fit$estimates
  # estimate number of change-points
  L_est <- length(est_cp) - 2
  
  # store results
  result[nrow(result) + 1,"method"] <- "ecp"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), est_cp), 
                               hausdorff(est_cp, c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), est_cp), 
                               fnsle(c(1, true_cp, T+1), est_cp))
  
  #### Inspect ####
  penalty <- log(log(T)*d/2)
  while(TRUE) {
    time <- try(system.time({fit <- inspect(t(cp_data$Y), lambda = penalty)})[3])
    if(!"try-error" %in% class(time)) break
    penalty <- penalty + 1
  }
  est_cp <- fit$changepoints[,1] + 1
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # store results
  result[nrow(result) + 1,"method"] <- "Inspect"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### L2HDC ####
  time <- system.time({ts_l2_fit <- ts_hdchange(t(cp_data$Y)); fit <- hdchange(ts_l2_fit);})[3]
  est_cp <- fit$time_stamps + 1
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # store results
  result[nrow(result) + 1,"method"] <- "L2HDC"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  #### MICH ####
  #### MICH Auto ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$Y, L_auto = TRUE, max_iter = Inf, tol = tol, restart = FALSE)})[3] 
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Auto Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$Y, L_auto = TRUE, max_iter = Inf, tol = tol, reverse = TRUE, restart = FALSE)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Auto rev"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Oracle ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$Y, L = L, max_iter = Inf, tol = tol)})[3] 
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Ora Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$Y, L = L, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora rev"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### Spatial - 0.7 #### 
  C <- 10
  rho <- 0.7
  adapt <- FALSE
  # generate data
  cp_data <- multi_simulation(T, L, d, p, rho, C, min_space, adapt)
  true_cp <- cp_data$changepoints
  
  #### ECP ####
  time <- system.time({fit <- e.divisive(cp_data$Y, sig.lvl = 1 - level, 
                                         min.size = min_space, alpha = 2,
                                         R=499)})[3]
  est_cp <- fit$estimates
  # estimate number of change-points
  L_est <- length(est_cp) - 2
  
  # store results
  result[nrow(result) + 1,"method"] <- "ecp"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), est_cp), 
                               hausdorff(est_cp, c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), est_cp), 
                               fnsle(c(1, true_cp, T+1), est_cp))
  
  #### Inspect ####
  penalty <- log(log(T)*d/2)
  while(TRUE) {
    time <- try(system.time({fit <- inspect(t(cp_data$Y), lambda = penalty)})[3])
    if(!("try-error" %in% class(time))) break
    penalty <- penalty + 1
  }
  est_cp <- fit$changepoints[,1] + 1
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # store results
  result[nrow(result) + 1,"method"] <- "Inspect"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### L2HDC ####
  time <- system.time({ts_l2_fit <- ts_hdchange(t(cp_data$Y)); fit <- hdchange(ts_l2_fit);})[3]
  est_cp <- fit$time_stamps + 1
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # store results
  result[nrow(result) + 1,"method"] <- "L2HDC"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, NA, NA, NA,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  #### MICH ####
  #### MICH Auto ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$Y, L_auto = TRUE, max_iter = Inf, tol = tol, restart = FALSE)})[3] 
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Auto Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$Y, L_auto = TRUE, max_iter = Inf, tol = tol, reverse = TRUE, restart = FALSE)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Auto rev"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Oracle ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$Y, L = L, max_iter = Inf, tol = tol)})[3] 
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
  
  #### MICH Ora Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$Y, L = L, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$L > 0) {
    post_sets <- mich_sets(fit$pi_bar, max_length, level)
    cs <- post_sets$sets
    est_cp <- post_sets$cp
  }    
  
  # estimate number of change-points
  L_est <- length(est_cp) 
  
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora rev"
  result[nrow(result),-1] <- c(L, d, p, rho, adapt, time, L_est, n_detected, 
                               n_covered, avg_len,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)))
}

write.csv(result, paste0("~/MICH/multi_results/result_",jobid,".csv"), row.names = FALSE)
