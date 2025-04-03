library(foreach)
library(doParallel)
library(dplyr)

registerDoParallel(8)

N <- 1000 # number of simulation repetitions
tol <- 1e-4 # convergence criterion

# simulation settings
settings <- data.frame(T = c(rep(100, 2), rep(500, 4), rep(1000, 4)),
                       L = c(2, 5, 2, 2, 5, 5, 2, 2, 10, 10),
                       min_space = c(15, 15, 15, 30, 15, 30, 30, 50, 30, 50))

cols <- c("method", "T", "L", "min_space", "time", "n_comp", "L_est", "n_detected", "n_covered",
          "avg_len", "mean_mse", "var_mse", "hausdorff", "fpsle", "fnsle", "overlap")

#### Simulation ####
results <- foreach(i = 1:N, .errorhandling = "pass", .inorder = FALSE) %dopar% {
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
    # distance between estimated change-points (preprocessed)
    chp_dist <- abs(outer(est_cp, est_cp, `-`))
    # store results
    result[nrow(result) + 1,"method"] <- "HSMUCE"
    result[nrow(result),-1] <- c(T, L, min_space, time, NA, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Auto ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = FALSE, fit_scale = FALSE)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar)
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Auto Rev ####
    # fit model and time process
    rev <- FALSE
    rev_dex <- 1:T
    time <- system.time({fit_rev <- mich(cp_data$y[T:1], J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = TRUE, fit_scale = TRUE)})[3] 
    if (max(fit$elbo) < max(fit_rev$elbo)) {
      rev_dex <- T:1
      fit <- fit_rev
      rev <- TRUE
    }
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal[rev_dex])^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal[rev_dex])^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar[rev_dex,,drop=FALSE])
      cs <- post_sets$sets
      est_cp <- post_sets$cp
      if (rev) {
        cs <- lapply(1:length(cs), function(i) cs[[i]] + 1)
        est_cp <- est_cp+1
      }
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit_rev$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    
    #### MICH (J,0,0) Auto - 1####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = FALSE, fit_scale = FALSE, n_restart = 1)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar)
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
    result[nrow(result) + 1,"method"] <- "MICH Auto - 1"
    result[nrow(result),-1] <- c(T, L, min_space, time, fit$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Auto Rev - 1 ####
    # fit model and time process
    rev <- FALSE
    rev_dex <- 1:T
    time <- system.time({fit_rev <- mich(cp_data$y[T:1], J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = TRUE, fit_scale = TRUE, n_restart = 1)})[3] 
    if (max(fit$elbo) < max(fit_rev$elbo)) {
      rev_dex <- T:1
      fit <- fit_rev
      rev <- TRUE
    }
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal[rev_dex])^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal[rev_dex])^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar[rev_dex,,drop=FALSE])
      cs <- post_sets$sets
      est_cp <- post_sets$cp
      if (rev) {
        cs <- lapply(1:length(cs), function(i) cs[[i]] + 1)
        est_cp <- est_cp+1
      }
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
    result[nrow(result) + 1,"method"] <- "MICH Auto Rev - 1"
    result[nrow(result),-1] <- c(T, L, min_space, time, fit_rev$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) AutoFast ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = FALSE, fit_scale = FALSE, restart = FALSE)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar)
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    
    #### MICH (J,0,0) AutoFast Rev ####
    # fit model and time process
    rev <- FALSE
    rev_dex <- 1:T
    time <- system.time({fit_rev <- mich(cp_data$y[T:1], J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = TRUE, fit_scale = TRUE, restart = FALSE)})[3] 
    if (max(fit$elbo) < max(fit_rev$elbo)) {
      rev_dex <- T:1
      fit <- fit_rev
      rev <- TRUE
    }
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal[rev_dex])^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal[rev_dex])^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar[rev_dex,,drop=FALSE])
      cs <- post_sets$sets
      est_cp <- post_sets$cp
      if (rev) {
        cs <- lapply(1:length(cs), function(i) cs[[i]] + 1)
        est_cp <- est_cp+1
      }
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit_rev$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Oracle ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J = L, tol = tol, max_iter = Inf, fit_intercept = FALSE, fit_scale = FALSE)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar)
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Ora Rev ####
    # fit model and time process
    rev <- FALSE
    rev_dex <- 1:T
    time <- system.time({fit_rev <- mich(cp_data$y[T:1], J = L, tol = tol, max_iter = Inf, fit_intercept = TRUE, fit_scale = TRUE)})[3] 
    if (max(fit$elbo) < max(fit_rev$elbo)) {
      rev_dex <- T:1
      fit <- fit_rev
      rev <- TRUE
    }
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal[rev_dex])^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal[rev_dex])^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar[rev_dex,,drop=FALSE])
      cs <- post_sets$sets
      est_cp <- post_sets$cp
      if (rev) {
        cs <- lapply(1:length(cs), function(i) cs[[i]] + 1)
        est_cp <- est_cp+1
      }
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit_rev$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
  }
  result
}

#saveRDS(results, "~/Desktop/auto_results.rds")
results <- readRDS("~/Desktop/auto_results.rds")

i = 10
tmp %>% 
  #filter(overlap == FALSE) %>% 
  filter(T == settings$T[i], L == settings$L[i], min_space == settings$min_space[i]) %>% 
  group_by(method, T, L, min_space) %>% 
  #group_by(method, T, L, min_space, overlap) %>% 
  summarize(n_comp = mean(n_comp),
            hausdorff = mean(hausdorff, na.rm = TRUE),
            fpsle = mean(fpsle, na.rm = TRUE),
            fnsle = mean(fnsle, na.rm = TRUE),
            "|L - L_hat|" = mean(abs(L - L_est)),
            "<= -2" = mean((L - L_est) <= -2),
            "= -1" = mean((L - L_est) == -1),
            "= 0" = mean((L - L_est) == 0),
            "= 1" = mean((L - L_est) == 1),
            ">= 2" = mean((L - L_est) >= 2),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            #mean_mse = mean(mean_mse, na.rm = TRUE),
            #var_mse = mean(var_mse, na.rm = TRUE),
            overlap = mean(overlap),
            #merge_prob_lms = mean(merge_prob_lms, na.rm = TRUE),
            #merge_prob_gms = mean(merge_prob_gms, na.rm = TRUE),
            time = mean(time, na.rm = TRUE)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  arrange(method, T, L, min_space) #%>% 

results <- matrix(ncol = length(cols), nrow = 0)
results <- data.frame(results)
names(results) <- cols

for (k in 1:N) {
  print(paste0("iter: ", k))
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
    # distance between estimated change-points (preprocessed)
    chp_dist <- abs(outer(est_cp, est_cp, `-`))
    # store results
    result[nrow(result) + 1,"method"] <- "HSMUCE"
    result[nrow(result),-1] <- c(T, L, min_space, time, NA, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Auto ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = FALSE, fit_scale = FALSE)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar)
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Auto Rev ####
    # fit model and time process
    rev <- FALSE
    rev_dex <- 1:T
    time <- system.time({fit_rev <- mich(cp_data$y[T:1], J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = TRUE, fit_scale = TRUE)})[3] 
    if (max(fit$elbo) < max(fit_rev$elbo)) {
      rev_dex <- T:1
      fit <- fit_rev
      rev <- TRUE
    }
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal[rev_dex])^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal[rev_dex])^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar[rev_dex,,drop=FALSE])
      cs <- post_sets$sets
      est_cp <- post_sets$cp
      if (rev) {
        cs <- lapply(1:length(cs), function(i) cs[[i]] + 1)
        est_cp <- est_cp+1
      }
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit_rev$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) AutoFast ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = FALSE, fit_scale = FALSE, restart = FALSE)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal)^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar)
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    
    #### MICH (J,0,0) AutoFast Rev ####
    # fit model and time process
    rev <- FALSE
    rev_dex <- 1:T
    time <- system.time({fit_rev <- mich(cp_data$y[T:1], J_auto = TRUE, tol = tol, max_iter = Inf, fit_intercept = TRUE, fit_scale = TRUE, restart = FALSE)})[3] 
    if (max(fit$elbo) < max(fit_rev$elbo)) {
      rev_dex <- T:1
      fit <- fit_rev
      rev <- TRUE
    }
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal[rev_dex])^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal[rev_dex])^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar[rev_dex,,drop=FALSE])
      cs <- post_sets$sets
      est_cp <- post_sets$cp
      if (rev) {
        cs <- lapply(1:length(cs), function(i) cs[[i]] + 1)
        est_cp <- est_cp+1
      }
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit_rev$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Oracle ####
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J = L, tol = tol, max_iter = Inf, fit_intercept = FALSE, fit_scale = FALSE)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar)
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
    
    #### MICH (J,0,0) Ora Rev ####
    # fit model and time process
    rev <- FALSE
    rev_dex <- 1:T
    time <- system.time({fit_rev <- mich(cp_data$y[T:1], J = L, tol = tol, max_iter = Inf, fit_intercept = TRUE, fit_scale = TRUE)})[3] 
    if (max(fit$elbo) < max(fit_rev$elbo)) {
      rev_dex <- T:1
      fit <- fit_rev
      rev <- TRUE
    }
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal[rev_dex])^2)
    var_mse <- mean((1 / fit$lambda - cp_data$var_signal[rev_dex])^2)
    
    # determine MICH credible sets and change-points
    cs <- numeric(0)
    est_cp <- numeric(0)
    if (fit$J > 0) {
      post_sets <- mich_sets(fit$J_model$pi_bar[rev_dex,,drop=FALSE])
      cs <- post_sets$sets
      est_cp <- post_sets$cp
      if (rev) {
        cs <- lapply(1:length(cs), function(i) cs[[i]] + 1)
        est_cp <- est_cp+1
      }
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
    result[nrow(result),-1] <- c(T, L, min_space, time, fit_rev$J, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T),
                                 any(chp_dist[lower.tri(chp_dist)] < min_space))
  }
  results <- rbind(results, result)
}
 