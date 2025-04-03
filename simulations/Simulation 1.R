library(foreach)
library(doParallel)
library(dplyr)

# change-point packages
library(blipr)
library(stepR)
library(changepoint)
library(mosum)
library(not)

registerDoParallel(8)

N <- 1000 # number of simulation repetitions
tol <- 1e-3 # convergence criterion

# simulation settings
settings <- data.frame(T = c(rep(100, 2), rep(500, 4), rep(1000, 4)),
                       L = c(2, 5, 2, 2, 5, 5, 2, 2, 10, 10),
                       min_space = c(15, 15, 15, 30, 15, 30, 30, 50, 30, 50))

cols <- c("method", "T", "L", "min_space", "time", "L_est", "n_detected", "n_covered",
          "avg_len", "mean_mse", "var_mse", "hausdorff", "fpsle", "fnsle")

#### Simulation ####

results <- foreach(i = 1:N, .combine = rbind, .verbose = TRUE) %dopar% {
  # initialize results matrix
  result <- matrix(ncol = length(cols), nrow = 0)
  result <- data.frame(result)
  names(result) <- cols
  
  for (j in 1:nrow(settings)) {
    # extract settings
    T <- settings$T[j]; L <- settings$L[j]; min_space <- settings$min_space[j];
    
    # setting priors
    if (T == 100) {
      pi_j <- pi_j_100
      pi_l <- pi_l_100
      pi_k <- pi_k_100
    } else if (T == 500) {
      pi_j <- pi_j_500
      pi_l <- pi_l_500
      pi_k <- pi_k_500
    } else {
      pi_j <- pi_j_1000
      pi_l <- pi_l_1000
      pi_k <- pi_k_1000
    }
    
    # conditional detection window
    window <- floor(min(sqrt(T) / 2, 15)) 
    
    # generate data
    cp_data <- hsmuce_simulation(T, L, 200, min_space)
    true_cp <- cp_data$changepoints
    
    #### MICH Models ####
    local_rate <- round(T / sqrt(T * log(T))) 
    
    #### MICH (0,L,K) Localization Rate - Uniform Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     L = local_rate, K = local_rate,
                                     pi_l = rep(1 / T, T), pi_k = rep(1 / T, T),
                                     tol = tol, max_iter = Inf)})[3] 

    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$L_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$L_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (0,L,K) LR - Uni."
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (0,L,K) Localization Rate - Inverse Weight Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     L = local_rate, K = local_rate,
                                     pi_l = pi_l, pi_k = pi_k,
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$L_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$L_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (0,L,K) LR"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (0,L,K) Oracle - Uniform Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     L = L, K = L, 
                                     pi_l = rep(1 / T, T), pi_k = rep(1 / T, T),
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$L_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$L_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (0,L,K) Oracle - Uni."
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (0,L,K) Oracle - Inverse Weight Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     L = L, K = L, 
                                     pi_l = pi_l, pi_k = pi_k,
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$L_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$L_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (0,L,K) Oracle - IW"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (0,L,K) Auto - Uniform Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     L_auto = TRUE, K_auto = TRUE, 
                                     pi_l = rep(1 / T, T), pi_k = rep(1 / T, T),
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    if (fit$L == 0) cs <- numeric(0)
    else cs <- mich_sets(fit$L_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    if (is.null(fit$L_model$pi)) est_cp <- numeric(0)
    else est_cp <- unique(apply(fit$L_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (0,L,K) Auto - Uni."
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (0,L,K) Auto - Inverse Weight Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     L_auto = TRUE, K_auto = TRUE, 
                                     pi_l = pi_l, pi_k = pi_k,
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    if (fit$L == 0) cs <- numeric(0)
    else cs <- mich_sets(fit$L_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    if (is.null(fit$L_model$pi)) est_cp <- numeric(0)
    else est_cp <- unique(apply(fit$L_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (0,L,K) Auto - IW"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Localization Rate - Uniform Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     J = local_rate, pi_j = rep(1 / T, T),
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) LR - Uni."
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Localization Rate - Inverse Weight Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     J = local_rate, pi_j = pi_j,
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) LR - IW"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Oracle - Uniform Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J = L, pi_j = rep(1 / T, T),
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Oracle - Uni."
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Oracle - Inverse Weight Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J = L, pi_j = pi_j,
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Oracle - IW"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Auto - Uniform Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     J_auto = TRUE, pi_j = rep(1 / T, T),
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    if (fit$J == 0) cs <- numeric(0)
    else cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)    
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    if (is.null(fit$J_model$pi)) est_cp <- numeric(0)
    else est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Auto - Uni."
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Auto - Inverse Weight Prior ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, 
                                     J_auto = TRUE, pi_j = pi_j,
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    if (is.null(fit$J_model$pi)) cs <- numeric(0)
    else cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)    
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    if (is.null(fit$J_model$pi)) est_cp <- numeric(0)
    else est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Auto - IW"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Fast-Auto - Uniform Prior ####
    # fit model and time process
    time <- system.time({fit <- mich2(cp_data$y, 
                                     J_auto = TRUE, pi_j = rep(1 / T, T),
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    if (fit$J == 0) cs <- numeric(0)
    else cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)    
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    if (is.null(fit$J_model$pi)) est_cp <- numeric(0)
    else est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Fast-Auto - Uni."
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MICH (J,0,0) Fast-Auto - Inverse Weight Prior ####
    # fit model and time process
    time <- system.time({fit <- mich2(cp_data$y, 
                                     J_auto = TRUE, pi_j = pi_j,
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    if (fit$J == 0) cs <- numeric(0)
    else cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4)    
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    if (is.null(fit$J_model$pi)) est_cp <- numeric(0)
    else est_cp <- unique(apply(fit$J_model$pi, 2, which.max)) 
    # only keep points in extracted credible sets
    est_cp <- est_cp[est_cp %in% unlist(cs)] 
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Fast-Auto - IW"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### H-SMUCE ####
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
    # store results
    result[nrow(result) + 1,"method"] <- "H-SMUCE"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
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
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, NA, NA, NA, 
                                 mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    #### MOSUM ####
    
    # bottom up merging ####
    # fit model and time process
    time <- system.time({fit = multiscale.bottomUp(cp_data$y, do.confint = TRUE, 
                                                   level = 0.1, 
                                                   var.est.method="mosum.min")})[3] 
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
    # store results
    result[nrow(result) + 1,"method"] <- "MOSUM BU"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
    # local pruning ####
    time <- system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, 
                                                     level = 0.1, 
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
    # store results
    result[nrow(result) + 1,"method"] <- "MOSUM LP"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    
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
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, NA, NA, NA, 
                                 mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
  }
  result
}

saveRDS(results, "~/Desktop/simulation_1_results.rds")

other_results <- readRDS("~/Desktop/simulation_1_results.rds")

methods = c("MICH (J,0,0) Oracle - IW",
            "MICH (J,0,0) Fast-Auto - IW",
            "MICH (J,0,0) Auto - IW",
            "H-SMUCE", "PELT", "MOSUM BU", "NOT")
i=10
other_results %>% 
  filter(method %in% methods) %>% 
  filter(T == settings$T[i], L == settings$L[i], min_space == settings$min_space[i]) %>% 
  group_by(method, T, L, min_space) %>% 
  summarize("|L - L_hat|" = mean(abs(L - L_est)),
            "<= -2" = mean((L - L_est) <= -2),
            "= -1" = mean((L - L_est) == -1),
            "= 0" = mean((L - L_est) == 0),
            "= 1" = mean((L - L_est) == 1),
            ">= 2" = mean((L - L_est) >= 2),
            ci_length = sum(L_est * avg_len) / sum(L_est),
            coverage = sum(n_covered) / sum(n_detected),
            hausdorff = mean(hausdorff, na.rm = TRUE),
            # fpsle = mean(fpsle, na.rm = TRUE),
            # fnsle = mean(fnsle, na.rm = TRUE),
            # mean_mse = mean(mean_mse, na.rm = TRUE),
            # var_mse = mean(var_mse, na.rm = TRUE),
            time = mean(time, na.rm = TRUE)) %>% 
  mutate_if(is.numeric, ~format(round(.,2), nsmall = 2)) %>% 
  arrange(T, L, min_space) #%>%
  #write.csv("~/Desktop/results.csv", row.names = FALSE)
  
settings[i,]

            
            
            

