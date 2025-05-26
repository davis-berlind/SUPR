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
delta <- 1.1
tol <- 1e-8

# Threshold for NSP
wn003 <- sim_max_holder(100, 500, 0.03)
nsp_thresh <- as.numeric(stats::quantile(wn003, 0.9))

# simulation settings
settings <- data.frame(T = c(rep(100, 2), rep(500, 2), rep(1000, 3)),
                       L = c(2, 5, 2, 5, 2, 5, 10),
                       min_space = c(15, 15, 30, 30, 50, 50, 50))

cols <- c("method", "family", "T", "L", "min_space", "time","L_est", "n_detected", 
          "n_covered", "avg_len", "mean_mse", "var_mse", "hausdorff_1", 
          "hausdorff_2", "fpsle", "fnsle", "overlap")

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
  
  #### laplace #### 
  # generate data
  cp_data <- hsmuce_simulation(T, L, 200, min_space, family = "laplace")
  true_cp <- cp_data$changepoints
  
  #### NSP ####
  # fit model and time process
  time <- system.time({fit <- nsp_selfnorm(cp_data$y, x = matrix(1, nrow = T), 
                                           M = 1000, lambda = nsp_thresh)})[3] 
  
  # estimate number of change-points
  est_cp <- fit$intervals$midpoints
  L_est <- nrow(fit$intervals)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) fit$intervals[i, "starts"]:fit$intervals[i, "ends"])
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))]
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  # check if sets overlap
  if (L_est > 1) overlap <- any(sapply(seq_len(L_est-1), function(i) fit$intervals[i, "ends"] > fit$intervals[i+1, "starts"]))
  else overlap <- 0
  
  # store results
  result[nrow(result) + 1,"method"] <- "NSP"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, NA, NA, NA, NA, NA, NA, overlap)
  
  #### NOT ####
  # fit model and time process
  time <- system.time({fit = not(cp_data$y, contrast = "pcwsConstMeanVar")})[3] 
  # estimate change-point locations
  est_cp <- features(fit)$cpt + 1
  # estimate number of change-points
  L_est <- length(est_cp)    
  # calculate mean/precision signal MSEs
  mean_mse <- mean((predict(fit)[,1] - cp_data$mean_signal)^2)
  var_mse <- mean((predict(fit)[,2]^2  - cp_data$var_signal)^2)

  # store results
  result[nrow(result) + 1,"method"] <- "NOT"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, NA, NA, NA, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### MOSUM BUM ####
  # fit model and time process
  time <- system.time({fit = multiscale.bottomUp(cp_data$y, do.confint = TRUE, 
                                                 alpha = 0.1, level = 0.1, 
                                                 var.est.method = "mosum.min")})[3] 
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))

  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM BUM"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MOSUM LP ####
  time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                   alpha = 0.1, level = 0.1, 
                                                   var.est.method = "mosum.min")})[3])
  if("try-error" %in% class(time)) {
    time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                         alpha = 0.1, level = 0.1, 
                                                         var.est.method = "mosum")})[3])
  }
    
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))

  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM LP"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  #### PELT ####
  # fit model and time process
  time <- system.time({fit <- cpt.meanvar(cp_data$y, method = "PELT")})[3] 
  # estimate number of change-points
  L_est <- length(fit@cpts) - 1
  # estimate change-point locations
  est_cp <- fit@cpts[-(L_est+1)] + 1
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit@param.est$mean, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <-  mean((rep(fit@param.est$variance, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)

  # store results
  result[nrow(result) + 1,"method"] <- "PELT"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, NA, NA, NA, 
                               mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### H-SMUCE 0.1 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.1, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.1)"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### H-SMUCE 0.5 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.5, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.5)"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH ####
  # detection threshold
  max_length <- 2*log(T)^1.5
  
  #### MICH (J,0,0) Auto ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, tol = tol)})[3] 
  
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
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH (J,0,0) Auto Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
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
  result[nrow(result) + 1,"method"] <- "MICH Auto Rev"
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
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
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("laplace", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  
  #### t #### 
  # generate data
  cp_data <- hsmuce_simulation(T, L, 200, min_space, family = "t", df = 4)
  true_cp <- cp_data$changepoints
  
  #### NSP ####
  # fit model and time process
  time <- system.time({fit <- nsp_selfnorm(cp_data$y, x = matrix(1, nrow = T), 
                                           M = 1000, lambda = nsp_thresh)})[3] 
  
  # estimate number of change-points
  est_cp <- fit$intervals$midpoints
  L_est <- nrow(fit$intervals)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) fit$intervals[i, "starts"]:fit$intervals[i, "ends"])
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))]
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  # check if sets overlap
  if (L_est > 1) overlap <- any(sapply(seq_len(L_est-1), function(i) fit$intervals[i, "ends"] > fit$intervals[i+1, "starts"]))
  else overlap <- 0
  
  # store results
  result[nrow(result) + 1,"method"] <- "NSP"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, NA, NA, NA, NA, NA, NA, overlap)
  
  #### NOT ####
  # fit model and time process
  time <- system.time({fit = not(cp_data$y, contrast = "pcwsConstMeanVar")})[3] 
  # estimate change-point locations
  est_cp <- features(fit)$cpt + 1
  # estimate number of change-points
  L_est <- length(est_cp)    
  # calculate mean/precision signal MSEs
  mean_mse <- mean((predict(fit)[,1] - cp_data$mean_signal)^2)
  var_mse <- mean((predict(fit)[,2]^2  - cp_data$var_signal)^2)
  
  # store results
  result[nrow(result) + 1,"method"] <- "NOT"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, NA, NA, NA, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### MOSUM BUM ####
  # fit model and time process
  time <- system.time({fit = multiscale.bottomUp(cp_data$y, do.confint = TRUE, 
                                                 alpha = 0.1, level = 0.1, 
                                                 var.est.method = "mosum.min")})[3] 
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM BUM"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MOSUM LP ####
  time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                       alpha = 0.1, level = 0.1, 
                                                       var.est.method = "mosum.min")})[3])
  if("try-error" %in% class(time)) {
    time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                         alpha = 0.1, level = 0.1, 
                                                         var.est.method = "mosum")})[3])
  }
  
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM LP"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  #### PELT ####
  # fit model and time process
  time <- system.time({fit <- cpt.meanvar(cp_data$y, method = "PELT")})[3] 
  # estimate number of change-points
  L_est <- length(fit@cpts) - 1
  # estimate change-point locations
  est_cp <- fit@cpts[-(L_est+1)] + 1
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit@param.est$mean, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <-  mean((rep(fit@param.est$variance, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  
  # store results
  result[nrow(result) + 1,"method"] <- "PELT"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, NA, NA, NA, 
                               mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### H-SMUCE 0.1 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.1, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.1)"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### H-SMUCE 0.5 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.5, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.5)"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH ####
  # detection threshold
  max_length <- 2*log(T)^1.5
  
  #### MICH (J,0,0) Auto ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, tol = tol)})[3] 
  
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
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH (J,0,0) Auto Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
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
  result[nrow(result) + 1,"method"] <- "MICH Auto Rev"
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
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
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("t", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  
  #### MA 0.3 #### 
  # generate data
  cp_data <- hsmuce_simulation(T, L, 200, min_space, family = "MA", theta = 0.3)
  true_cp <- cp_data$changepoints
  
  #### NSP ####
  # fit model and time process
  time <- system.time({fit <- nsp_selfnorm(cp_data$y, x = matrix(1, nrow = T), 
                                           M = 1000, lambda = nsp_thresh)})[3] 
  
  # estimate number of change-points
  est_cp <- fit$intervals$midpoints
  L_est <- nrow(fit$intervals)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) fit$intervals[i, "starts"]:fit$intervals[i, "ends"])
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))]
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  # check if sets overlap
  if (L_est > 1) overlap <- any(sapply(seq_len(L_est-1), function(i) fit$intervals[i, "ends"] > fit$intervals[i+1, "starts"]))
  else overlap <- 0
  
  # store results
  result[nrow(result) + 1,"method"] <- "NSP"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, NA, NA, NA, NA, NA, NA, overlap)
  
  #### NOT ####
  # fit model and time process
  time <- system.time({fit = not(cp_data$y, contrast = "pcwsConstMeanVar")})[3] 
  # estimate change-point locations
  est_cp <- features(fit)$cpt + 1
  # estimate number of change-points
  L_est <- length(est_cp)    
  # calculate mean/precision signal MSEs
  mean_mse <- mean((predict(fit)[,1] - cp_data$mean_signal)^2)
  var_mse <- mean((predict(fit)[,2]^2  - cp_data$var_signal)^2)
  
  # store results
  result[nrow(result) + 1,"method"] <- "NOT"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, NA, NA, NA, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### MOSUM BUM ####
  # fit model and time process
  time <- system.time({fit = multiscale.bottomUp(cp_data$y, do.confint = TRUE, 
                                                 alpha = 0.1, level = 0.1, 
                                                 var.est.method = "mosum.min")})[3] 
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM BUM"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MOSUM LP ####
  time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                       alpha = 0.1, level = 0.1, 
                                                       var.est.method = "mosum.min")})[3])
  if("try-error" %in% class(time)) {
    time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                         alpha = 0.1, level = 0.1, 
                                                         var.est.method = "mosum")})[3])
  }
  
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM LP"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  #### PELT ####
  # fit model and time process
  time <- system.time({fit <- cpt.meanvar(cp_data$y, method = "PELT")})[3] 
  # estimate number of change-points
  L_est <- length(fit@cpts) - 1
  # estimate change-point locations
  est_cp <- fit@cpts[-(L_est+1)] + 1
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit@param.est$mean, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <-  mean((rep(fit@param.est$variance, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  
  # store results
  result[nrow(result) + 1,"method"] <- "PELT"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, NA, NA, NA, 
                               mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### H-SMUCE 0.1 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.1, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.1)"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### H-SMUCE 0.5 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.5, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.5)"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH ####
  # detection threshold
  max_length <- 2*log(T)^1.5
  
  #### MICH (J,0,0) Auto ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, tol = tol)})[3] 
  
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
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH (J,0,0) Auto Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
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
  result[nrow(result) + 1,"method"] <- "MICH Auto Rev"
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
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
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("MA 0.3", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  
  #### MA 0.7 #### 
  # generate data
  cp_data <- hsmuce_simulation(T, L, 200, min_space, family = "MA", theta = 0.7)
  true_cp <- cp_data$changepoints
  
  #### NSP ####
  # fit model and time process
  time <- system.time({fit <- nsp_selfnorm(cp_data$y, x = matrix(1, nrow = T), 
                                           M = 1000, lambda = nsp_thresh)})[3] 
  
  # estimate number of change-points
  est_cp <- fit$intervals$midpoints
  L_est <- nrow(fit$intervals)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) fit$intervals[i, "starts"]:fit$intervals[i, "ends"])
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))]
  n_detected <- length(detected)
  # number of credible sets that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average credible set length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
  # check if sets overlap
  if (L_est > 1) overlap <- any(sapply(seq_len(L_est-1), function(i) fit$intervals[i, "ends"] > fit$intervals[i+1, "starts"]))
  else overlap <- 0
  
  # store results
  result[nrow(result) + 1,"method"] <- "NSP"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, NA, NA, NA, NA, NA, NA, overlap)
  
  #### NOT ####
  # fit model and time process
  time <- system.time({fit = not(cp_data$y, contrast = "pcwsConstMeanVar")})[3] 
  # estimate change-point locations
  est_cp <- features(fit)$cpt + 1
  # estimate number of change-points
  L_est <- length(est_cp)    
  # calculate mean/precision signal MSEs
  mean_mse <- mean((predict(fit)[,1] - cp_data$mean_signal)^2)
  var_mse <- mean((predict(fit)[,2]^2  - cp_data$var_signal)^2)
  
  # store results
  result[nrow(result) + 1,"method"] <- "NOT"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, NA, NA, NA, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### MOSUM BUM ####
  # fit model and time process
  time <- system.time({fit = multiscale.bottomUp(cp_data$y, do.confint = TRUE, 
                                                 alpha = 0.1, level = 0.1, 
                                                 var.est.method = "mosum.min")})[3] 
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM BUM"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MOSUM LP ####
  time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                       alpha = 0.1, level = 0.1, 
                                                       var.est.method = "mosum.min")})[3])
  if("try-error" %in% class(time)) {
    time <- try(system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                         alpha = 0.1, level = 0.1, 
                                                         var.est.method = "mosum")})[3])
  }
  
  # estimate change-point locations
  est_cp <- fit$cpts + 1
  # estimate number of change-points
  L_est <- length(est_cp)
  # calculate mean/precision signal MSEs
  mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
  mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  # determine confidence intervals
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
  }
  # detected true change-points
  detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
  n_detected <- length(detected)
  # number of CIs that cover a detected change-point
  n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs)))
  # average CI length
  avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
  # distance between estimated change-points (preprocessed)
  chp_dist <- abs(outer(est_cp, est_cp, `-`))
  
  # store results
  result[nrow(result) + 1,"method"] <- "MOSUM LP"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  #### PELT ####
  # fit model and time process
  time <- system.time({fit <- cpt.meanvar(cp_data$y, method = "PELT")})[3] 
  # estimate number of change-points
  L_est <- length(fit@cpts) - 1
  # estimate change-point locations
  est_cp <- fit@cpts[-(L_est+1)] + 1
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit@param.est$mean, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
  var_mse <-  mean((rep(fit@param.est$variance, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
  
  # store results
  result[nrow(result) + 1,"method"] <- "PELT"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, NA, NA, NA, 
                               mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               NA)
  
  #### H-SMUCE 0.1 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.1, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.1)"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### H-SMUCE 0.5 ####
  # fit model and time process
  time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.5, 
                                      jumpint = TRUE, family = "hsmuce")})[3] 
  # calculate mean/precision signal MSEs
  mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
  hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
  var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
  # estimate change-point locations
  est_cp <- fit$leftEnd[-1]
  # estimate number of change-points
  L_est <- length(est_cp)
  # determine credible sets
  if (L_est == 0) {
    cs <- numeric(0)
  } else {
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
  }
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
  result[nrow(result) + 1,"method"] <- "HSMUCE (0.5)"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)), 
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH ####
  # detection threshold
  max_length <- log(T)^(1+delta)
  
  #### MICH (J,0,0) Auto ####
  # fit model and time process
  time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, tol = tol)})[3] 
  
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
  result[nrow(result) + 1,"method"] <- "MICH Auto"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
  #### MICH (J,0,0) Auto Rev ####
  # fit model and time process
  time_rev <- time
  time <- system.time({fit_rev <- mich(cp_data$y, J_auto = TRUE, max_iter = Inf, reverse = TRUE, tol = tol)})[3] 
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
  result[nrow(result) + 1,"method"] <- "MICH Auto Rev"
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
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
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
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
  result[nrow(result),-1] <- c("MA 0.7", T, L, min_space, time, L_est, n_detected, 
                               n_covered, avg_len, mean_mse, var_mse,
                               hausdorff(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               hausdorff(c(1, est_cp, T+1), c(1, true_cp, T+1)),
                               fpsle(c(1, true_cp, T+1), c(1, est_cp, T+1)), 
                               fnsle(c(1, true_cp, T+1), c(1, est_cp, T+1)),
                               any(chp_dist[lower.tri(chp_dist)] < min_space / 2))
  
}

write.csv(result, paste0("~/MICH/heavy_results/result_",jobid,".csv"), row.names = FALSE)
