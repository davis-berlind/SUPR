# get job id and set seed
jobid = as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(jobid)

# change-point packages
cp_pkgs <- c("stepR", "changepoint", "mosum", "not", "nsp")
lapply(cp_pkgs, require, character.only = TRUE)

# load MICH functions
source("~/MICH/mich.R")
source("~/MICH/mich_vector.R")
source("~/MICH/mich_matrix.R")
source("~/MICH/credible_sets.R")
source("~/MICH/priors.R")
source("~/MICH/utilities.R")
source("~/MICH/simulation_functions.R")
Rcpp::sourceCpp("mich.cpp")
Rcpp::sourceCpp("multi_mich.cpp")
Rcpp::sourceCpp("force_adj.cpp")

# significance level
level <- 0.9
tol <- 1e-5
prior <- 1e-3
delta <- 0.5

# simulation settings
settings <- data.frame(T = c(rep(100, 2), rep(500, 2), rep(1000, 3)),
                       L = c(2, 5, 2, 5, 2, 5, 10),
                       min_space = c(15, 15, 30, 30, 50, 50, 50))

cols <- c("method", "T", "L", "min_space", "prior", "tol", "delta", "time", 
          "n_comp", "L_est", "n_detected", "n_covered", "avg_len", "mean_mse", 
          "var_mse", "hausdorff_1", "hausdorff_2", "fpsle", "fnsle", "overlap")

# initialize results matrix
result <- matrix(ncol = length(cols), nrow = 0)
result <- data.frame(result)
names(result) <- cols

for (j in 1:nrow(settings)) {
  print(j)
  # extract settings
  T <- settings$T[j]; L <- settings$L[j]; min_space <- settings$min_space[j];
  
  # conditional detection window
  window <- floor(min(sqrt(T) / 2, 15)) 
  
  # generate data
  cp_data <- hsmuce_simulation(T, L, 200, min_space)
  true_cp <- cp_data$changepoints
  
  #### MICH ####
  # detection threshold
  max_length <- log(T)^(1 + delta)
  
  #### MICH (J,0,0) AutoFast ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, restart = FALSE, tol = tol)})[3] 
  
  # calculate mean/precision signal MSEs
  mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
  var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$J > 0) {
    post_sets <- mich_sets(fit$meanvar_model$pi_bar, max_length, level)
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
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH AutoFast"
  result[nrow(result),-1] <- c(T, L, min_space, prior, tol, delta, 
                               time, fit$J, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH (J,0,0) AutoFast Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, 
                                       restart = FALSE, reverse = TRUE, tol = tol)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # calculate mean/precision signal MSEs
  mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
  var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$J > 0) {
    post_sets <- mich_sets(fit$meanvar_model$pi_bar, max_length, level)
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
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH AutoFast Rev"
  result[nrow(result),-1] <- c(T, L, min_space, prior, tol, delta, 
                               time, fit$J, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH (J,0,0) Oracle ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$y, J = L, max_iter = Inf, tol = tol)})[3] 
  
  # calculate mean/precision signal MSEs
  mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
  var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$J > 0) {
    post_sets <- mich_sets(fit$meanvar_model$pi_bar, max_length, level)
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
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora"
  result[nrow(result),-1] <- c(T, L, min_space, prior, tol, delta, 
                               time, fit$J, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH (J,0,0) Ora Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$y, J = L, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
  time <- max(time, time_rev)
  if (max(fit$elbo) < max(fit_rev$elbo)) {
    fit <- fit_rev
  }
  
  # calculate mean/precision signal MSEs
  mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
  var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
  
  # determine MICH credible sets and change-points
  cs <- numeric(0)
  est_cp <- numeric(0)
  if (fit$J > 0) {
    post_sets <- mich_sets(fit$meanvar_model$pi_bar, max_length, level)
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
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MICH Ora Rev"
  result[nrow(result),-1] <- c(T, L, min_space, prior, tol, delta, 
                               time, fit$J, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
}

write.csv(result, paste0("~/MICH/e5_results/result_",jobid,".csv"), row.names = FALSE)
