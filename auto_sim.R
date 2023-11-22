library(foreach)
library(doParallel)
library(dplyr)

registerDoParallel(8)

pi_j_100 <- prior_j(100, 0, 5e5, 1e-3, 1e-3, 1e-3)
pi_j_500 <- prior_j(500, 0, 5e5, 1e-3, 1e-3, 1e-3)
pi_j_1000 <- prior_j(1000, 0, 5e5, 1e-3, 1e-3, 1e-3)

N <- 1000 # number of simulation repetitions
tol <- 1e-3 # convergence criterion

# simulation settings
settings <- data.frame(T = c(rep(100, 2), rep(500, 4), rep(1000, 4)),
                       L = c(2, 5, 2, 2, 5, 5, 2, 2, 10, 10),
                       min_space = c(15, 15, 15, 30, 15, 30, 30, 50, 30, 50))

cols <- c("method", "T", "L", "min_space", "time", "L_est", "n_detected", "n_covered",
          "avg_len", "mean_mse", "var_mse", "hausdorff", "fpsle", "fnsle", "overlap", 
          "merge_prob_lms", "merge_prob_gms")

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
    if (T == 100) pi_j <- pi_j_100
    else if (T == 500) pi_j <- pi_j_500 
    else pi_j <- pi_j_1000
    
    # conditional detection window
    window <- floor(min(sqrt(T) / 2, 15)) 
    
    # generate data
    cp_data <- hsmuce_simulation(T, L, 200, min_space)
    true_cp <- cp_data$changepoints
    
    #### MICH (J,0,0) Auto ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- numeric(0)
    if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)    
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp_raw <- numeric(0)
    if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    est_cp <- numeric(0)
    if (L_est > 0) {
      for(i in 1:L_est) {
        set = cs[[i]]
        if (length(set) == 1) est_cp <- c(set, est_cp)
        else {
          if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
          else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
        }
      }
    }
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
    # distance between estimated change-points (preprocessed)
    chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # pairwise probability of equivalence
    lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    merge_prob_lms <- NA
    merge_prob_gms <- NA
    if (any(chp_dist < min_space & lower_mask)) {
      merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    }
    if (any(chp_dist >= min_space & lower_mask)) {
      merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    }
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Auto"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower_mask] < min_space),
                                 merge_prob_lms, merge_prob_gms)
    
    # #### MICH2 (J,0,0) Auto ####
    # # fit model and time process
    # time <- system.time({fit <- mich2(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
    #                                  tol = tol, max_iter = Inf)})[3] 
    # 
    # # calculate mean/precision signal MSEs
    # mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    # var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # # determine MICH credible sets
    # cs <- numeric(0)
    # if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)     
    # # estimate number of change-points
    # L_est <- length(cs) 
    # # extract MAP change-point locations
    # est_cp_raw <- numeric(0)
    # if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    # est_cp <- numeric(0)
    # if (L_est > 0) {
    #   for(i in 1:L_est) {
    #     set = cs[[i]]
    #     if (length(set) == 1) est_cp <- c(set, est_cp)
    #     else {
    #       if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
    #       else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
    #     }
    #   }
    # }
    # # detected true change-points
    # detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    # n_detected <- length(detected)
    # # number of credible sets that cover a detected change-point
    # n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # # average credible set length
    # avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # # distance between estimated change-points (preprocessed)
    # chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # # pairwise probability of equivalence
    # lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    # if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    # merge_prob_lms <- NA
    # merge_prob_gms <- NA
    # if (any(chp_dist < min_space & lower_mask)) {
    #   merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    # }
    # if (any(chp_dist >= min_space & lower_mask)) {
    #   merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    # }
    # # store results
    # result[nrow(result) + 1,"method"] <- "MICH2 (J,0,0) Auto"
    # result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
    #                              n_covered, avg_len, mean_mse, var_mse,
    #                              hausdorff(true_cp, est_cp, T), 
    #                              fpsle(true_cp, est_cp, T), 
    #                              fnsle(true_cp, est_cp, T),
    #                              any(chp_dist[lower_mask] < min_space),
    #                              merge_prob_lms, merge_prob_gms)
    # 
    # #### MICH3 (J,0,0) Auto ####
    # # fit model and time process
    # time <- system.time({fit <- mich3(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
    #                                  tol = tol, max_iter = Inf)})[3] 
    # 
    # # calculate mean/precision signal MSEs
    # mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    # var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # # determine MICH credible sets
    # cs <- numeric(0)
    # if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)     
    # # estimate number of change-points
    # L_est <- length(cs) 
    # # extract MAP change-point locations
    # est_cp_raw <- numeric(0)
    # if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    # est_cp <- numeric(0)
    # if (L_est > 0) {
    #   for(i in 1:L_est) {
    #     set = cs[[i]]
    #     if (length(set) == 1) est_cp <- c(set, est_cp)
    #     else {
    #       if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
    #       else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
    #     }
    #   }
    # }
    # # detected true change-points
    # detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    # n_detected <- length(detected)
    # # number of credible sets that cover a detected change-point
    # n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # # distance between estimated change-points (preprocessed)
    # chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # # pairwise probability of equivalence
    # lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    # if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    # merge_prob_lms <- NA
    # merge_prob_gms <- NA
    # if (any(chp_dist < min_space & lower_mask)) {
    #   merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    # }
    # if (any(chp_dist >= min_space & lower_mask)) {
    #   merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    # }
    # # store results
    # result[nrow(result) + 1,"method"] <- "MICH3 (J,0,0) Auto"
    # result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
    #                              n_covered, avg_len, mean_mse, var_mse,
    #                              hausdorff(true_cp, est_cp, T), 
    #                              fpsle(true_cp, est_cp, T), 
    #                              fnsle(true_cp, est_cp, T),
    #                              any(chp_dist[lower_mask] < min_space),
    #                              merge_prob_lms, merge_prob_gms)
    # 
    #### MICH4 (J,0,0) Auto ####
    # fit model and time process
    time <- system.time({fit <- mich4(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- numeric(0)
    if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)     
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp_raw <- numeric(0)
    if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    est_cp <- numeric(0)
    if (L_est > 0) {
      for(i in 1:L_est) {
        set = cs[[i]]
        if (length(set) == 1) est_cp <- c(set, est_cp)
        else {
          if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
          else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
        }
      }
    }
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # distance between estimated change-points (preprocessed)
    chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # pairwise probability of equivalence
    lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    merge_prob_lms <- NA
    merge_prob_gms <- NA
    if (any(chp_dist < min_space & lower_mask)) {
      merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    }
    if (any(chp_dist >= min_space & lower_mask)) {
      merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    }
    # store results
    result[nrow(result) + 1,"method"] <- "MICH4 (J,0,0) Auto"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower_mask] < min_space),
                                 merge_prob_lms, merge_prob_gms)
  }
  result
}


results %>% 
  #filter(overlap == FALSE) %>% 
  #filter(T == settings$T[i], L == settings$L[i], min_space == settings$min_space[i]) %>% 
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
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            #mean_mse = mean(mean_mse, na.rm = TRUE),
            #var_mse = mean(var_mse, na.rm = TRUE),
            overlap = mean(overlap),
            merge_prob_lms = mean(merge_prob_lms, na.rm = TRUE),
            merge_prob_gms = mean(merge_prob_gms, na.rm = TRUE),
            time = mean(time, na.rm = TRUE)) %>% 
  mutate_if(is.numeric, ~format(round(.,2), nsmall = 2)) %>% 
  arrange(method, T, L, min_space) #%>% 

results <- matrix(ncol = length(cols), nrow = 0)
results <- data.frame(results)
names(results) <- cols

for (k in 1:N) {
  # initialize results matrix
  result <- matrix(ncol = length(cols), nrow = 0)
  result <- data.frame(result)
  names(result) <- cols
  
  for (j in 1:nrow(settings)) {
    # extract settings
    T <- settings$T[j]; L <- settings$L[j]; min_space <- settings$min_space[j];
    
    # setting priors
    if (T == 100) pi_j <- pi_j_100
    else if (T == 500) pi_j <- pi_j_500 
    else pi_j <- pi_j_1000
    
    # conditional detection window
    window <- floor(min(sqrt(T) / 2, 15)) 
    
    # generate data
    cp_data <- hsmuce_simulation(T, L, 200, min_space)
    true_cp <- cp_data$changepoints
    
    #### MICH (J,0,0) Auto ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- numeric(0)
    if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)    
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp_raw <- numeric(0)
    if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    est_cp <- numeric(0)
    if (L_est > 0) {
      for(i in 1:L_est) {
        set = cs[[i]]
        if (length(set) == 1) est_cp <- c(set, est_cp)
        else {
          if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
          else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
        }
      }
    }
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x))))
    # distance between estimated change-points (preprocessed)
    chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # pairwise probability of equivalence
    lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    merge_prob_lms <- NA
    merge_prob_gms <- NA
    if (any(chp_dist < min_space & lower_mask)) {
      merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    }
    if (any(chp_dist >= min_space & lower_mask)) {
      merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    }
    # store results
    result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Auto"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower_mask] < min_space),
                                 merge_prob_lms, merge_prob_gms)
    
    # #### MICH2 (J,0,0) Auto ####
    # # fit model and time process
    # time <- system.time({fit <- mich2(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
    #                                  tol = tol, max_iter = Inf)})[3] 
    # 
    # # calculate mean/precision signal MSEs
    # mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    # var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # # determine MICH credible sets
    # cs <- numeric(0)
    # if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)     
    # # estimate number of change-points
    # L_est <- length(cs) 
    # # extract MAP change-point locations
    # est_cp_raw <- numeric(0)
    # if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    # est_cp <- numeric(0)
    # if (L_est > 0) {
    #   for(i in 1:L_est) {
    #     set = cs[[i]]
    #     if (length(set) == 1) est_cp <- c(set, est_cp)
    #     else {
    #       if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
    #       else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
    #     }
    #   }
    # }
    # # detected true change-points
    # detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    # n_detected <- length(detected)
    # # number of credible sets that cover a detected change-point
    # n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # # average credible set length
    # avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # # distance between estimated change-points (preprocessed)
    # chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # # pairwise probability of equivalence
    # lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    # if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    # merge_prob_lms <- NA
    # merge_prob_gms <- NA
    # if (any(chp_dist < min_space & lower_mask)) {
    #   merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    # }
    # if (any(chp_dist >= min_space & lower_mask)) {
    #   merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    # }
    # # store results
    # result[nrow(result) + 1,"method"] <- "MICH2 (J,0,0) Auto"
    # result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
    #                              n_covered, avg_len, mean_mse, var_mse,
    #                              hausdorff(true_cp, est_cp, T), 
    #                              fpsle(true_cp, est_cp, T), 
    #                              fnsle(true_cp, est_cp, T),
    #                              any(chp_dist[lower_mask] < min_space),
    #                              merge_prob_lms, merge_prob_gms)
    # 
    # #### MICH3 (J,0,0) Auto ####
    # # fit model and time process
    # time <- system.time({fit <- mich3(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
    #                                  tol = tol, max_iter = Inf)})[3] 
    # 
    # # calculate mean/precision signal MSEs
    # mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    # var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # # determine MICH credible sets
    # cs <- numeric(0)
    # if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)     
    # # estimate number of change-points
    # L_est <- length(cs) 
    # # extract MAP change-point locations
    # est_cp_raw <- numeric(0)
    # if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    # est_cp <- numeric(0)
    # if (L_est > 0) {
    #   for(i in 1:L_est) {
    #     set = cs[[i]]
    #     if (length(set) == 1) est_cp <- c(set, est_cp)
    #     else {
    #       if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
    #       else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
    #     }
    #   }
    # }
    # # detected true change-points
    # detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    # n_detected <- length(detected)
    # # number of credible sets that cover a detected change-point
    # n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # # distance between estimated change-points (preprocessed)
    # chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # # pairwise probability of equivalence
    # lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    # if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    # merge_prob_lms <- NA
    # merge_prob_gms <- NA
    # if (any(chp_dist < min_space & lower_mask)) {
    #   merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    # }
    # if (any(chp_dist >= min_space & lower_mask)) {
    #   merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    # }
    # # store results
    # result[nrow(result) + 1,"method"] <- "MICH3 (J,0,0) Auto"
    # result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
    #                              n_covered, avg_len, mean_mse, var_mse,
    #                              hausdorff(true_cp, est_cp, T), 
    #                              fpsle(true_cp, est_cp, T), 
    #                              fnsle(true_cp, est_cp, T),
    #                              any(chp_dist[lower_mask] < min_space),
    #                              merge_prob_lms, merge_prob_gms)
    # 
    #### MICH4 (J,0,0) Auto ####
    # fit model and time process
    time <- system.time({fit <- mich4(cp_data$y, J_auto = TRUE, pi_j = pi_j, 
                                      tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    # determine MICH credible sets
    cs <- numeric(0)
    if (fit$J > 0) cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = 1 / T)     
    # estimate number of change-points
    L_est <- length(cs) 
    # extract MAP change-point locations
    est_cp_raw <- numeric(0)
    if (L_est > 0 ) est_cp_raw <- apply(fit$J_model$pi, 2, which.max)
    est_cp <- numeric(0)
    if (L_est > 0) {
      for(i in 1:L_est) {
        set = cs[[i]]
        if (length(set) == 1) est_cp <- c(set, est_cp)
        else {
          if (L_est == 1) est_cp <- c(est_cp, which.max(fit$J_model$pi))
          else est_cp <- c(est_cp, set[which.max(fit$J_model$pi[set, which.max(apply(fit$J_model$pi[set,], 2, max))])])
        }
      }
    }
    # detected true change-points
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] 
    n_detected <- length(detected)
    # number of credible sets that cover a detected change-point
    n_covered <- ifelse(L_est == 0, 0, sum(detected %in% unlist(cs))) 
    # average credible set length
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) 
    # distance between estimated change-points (preprocessed)
    chp_dist <- abs(outer(est_cp_raw, est_cp_raw, `-`))
    # pairwise probability of equivalence
    lower_mask <- lower.tri(matrix(ncol = fit$J, nrow = fit$J))
    if (fit$J > 1) eq_prob <- t(fit$J_model$pi) %*% fit$J_model$pi
    merge_prob_lms <- NA
    merge_prob_gms <- NA
    if (any(chp_dist < min_space & lower_mask)) {
      merge_prob_lms <- mean(eq_prob[chp_dist < min_space & lower_mask])
    }
    if (any(chp_dist >= min_space & lower_mask)) {
      merge_prob_gms <- mean(eq_prob[chp_dist >= min_space & lower_mask])
    }
    # store results
    result[nrow(result) + 1,"method"] <- "MICH4 (J,0,0) Auto"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower_mask] < min_space),
                                 merge_prob_lms, merge_prob_gms)
  }
  results <- rbind(results, result)
  print(k)
}
 