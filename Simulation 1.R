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

N <- 1000
settings <- data.frame(T = c(rep(100, 2), rep(500, 4), rep(1000, 4)),
                       L = c(2, 5, 2, 2, 5, 5, 2, 2, 10, 10),
                       min_space = c(15, 15, 15, 30, 15, 30, 30, 50, 30, 50))

cols <- c("method", "T", "L", "min_space", "time", "L_est", "n_detected", "n_covered",
          "avg_len", "mean_mse", "var_mse", "hausdorff", "fpsle", "fnsle")

#### Simulation ####

results <- foreach(i = 1:N, .combine = rbind, .verbose = TRUE) %dopar% {
  result <- matrix(ncol = length(cols), nrow = 0)
  result <- data.frame(result)
  names(result) <- cols
  for (j in 1:nrow(settings)) {
    T <- settings$T[j]; L <- settings$L[j]; min_space <- settings$min_space[j];
    window <- floor(min(sqrt(T)/2, 15)) # detection window
    
    # generate data
    cp_data <- hsmuce_simulation(T, L, 200, min_space)
    true_cp <- cp_data$changepoints
    
    #### MICH I ####
    local_rate <- round(T / sqrt(T * log(T))) 
    time <- system.time({fit <- mich_i(cp_data$y, L = local_rate , K = local_rate, 
                                       tol = 0.001, fit.intercept = FALSE)})[3] # fit model and time process
    mean_mse <- mean((fit$beta - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    
    # MICH CS
    cs <- cred_set_prisca(fit$pi, max_length = T %/% 4)
    L_est <- length(cs) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(cs)] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(cs, function(x) any(detected %in% x)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH I"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    # BLiP CS
    blip_groups <- suppressWarnings(susie_groups(t(fit$pi), X = NULL, q = 0.1)) # generate blip candidate groups
    blip_sets <- BLiP(cand_groups = blip_groups) # generate blip sets
    L_est <- length(blip_sets) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(sapply(blip_sets, function(x) x$group))] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(blip_sets, function(x) any(detected %in% x$group)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(blip_sets, function(x) length(x$group)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH I + BLiP"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    #### MICH I oracle ####
    time <- system.time({fit <- mich_i(cp_data$y, L = L, K = L, tol = 0.001, fit.intercept = FALSE)})[3] # fit model and time process
    mean_mse <- mean((fit$beta - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    
    # MICH CS
    cs <- cred_set_prisca(fit$pi, max_length = T %/% 4)
    L_est <- length(cs) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(cs)] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(cs, function(x) any(detected %in% x)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH I Oracle"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    # BLiP CS
    blip_groups <- suppressWarnings(susie_groups(t(fit$pi), X = NULL, q = 0.1)) # generate blip candidate groups
    blip_sets <- BLiP(cand_groups = blip_groups) # generate blip sets
    L_est <- length(blip_sets) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(sapply(blip_sets, function(x) x$group))] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(blip_sets, function(x) any(detected %in% x$group)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(blip_sets, function(x) length(x$group)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH I + BLiP Oracle"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    #### MICH II ####
    time <- system.time({fit <- mich_ii(cp_data$y, L = local_rate, tol = 0.001, fit.intercept = FALSE)})[3] # fit model and time process
    mean_mse <- mean((fit$beta - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    
    # MICH CS
    cs <- cred_set_prisca(fit$pi, max_length = T %/% 4)
    L_est <- length(cs) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(cs)] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(cs, function(x) any(detected %in% x)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH II"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    # BLiP CS
    blip_groups <- suppressWarnings(susie_groups(t(fit$pi), X = NULL, q = 0.1)) # generate blip candidate groups
    blip_sets <- BLiP(cand_groups = blip_groups) # generate blip sets
    L_est <- length(blip_sets) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(sapply(blip_sets, function(x) x$group))] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(blip_sets, function(x) any(detected %in% x$group)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(blip_sets, function(x) length(x$group)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH II + BLiP"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    #### MICH II oracle ####
    time <- system.time({fit <- mich_ii(cp_data$y, L = L, tol = 0.001, fit.intercept = FALSE)})[3] # fit model and time process
    mean_mse <- mean((fit$beta - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    
    # MICH CS
    cs <- cred_set_prisca(fit$pi, max_length = T %/% 4)
    L_est <- length(cs) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(cs)] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(cs, function(x) any(detected %in% x)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH II Oracle"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    # BLiP CS
    blip_groups <- suppressWarnings(susie_groups(t(fit$pi), X = NULL, q = 0.1)) # generate blip candidate groups
    blip_sets <- BLiP(cand_groups = blip_groups) # generate blip sets
    L_est <- length(blip_sets) 
    est_cp <- unique(apply(fit$pi, 2, which.max)) # extract MAP change-point locations
    est_cp <- est_cp[est_cp %in% unlist(sapply(blip_sets, function(x) x$group))] # only keep points in extracted CS
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(blip_sets, function(x) any(detected %in% x$group)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(blip_sets, function(x) length(x$group)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MICH II + BLiP Oracle"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    #### H-SMUCE ####
    time <- system.time({fit <- stepFit(cp_data$y, alpha = 0.1, jumpint = TRUE, family = "hsmuce")})[3] # fit model and time process
    mean_mse <- mean((rep(fit$value, diff(c(fit$leftEnd,T+1))) - cp_data$mean_signal)^2)
    hsmuce_var <- sapply(1:length(fit$leftEnd), function(i) var(cp_data$y[fit$leftEnd[i]:fit$rightEnd[i]]))
    var_mse <-  mean((rep(hsmuce_var, diff(c(fit$leftEnd,T+1))) - cp_data$var_signal)^2)
    est_cp <- fit$leftEnd[-1]
    L_est <- length(est_cp)
    cs <- lapply(1:L_est, function(i) (fit$leftEndLeftBound[i+1]-1):fit$leftEndRightBound[i+1])
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(cs, function(x) any(detected %in% x)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) # average CS length
    result[nrow(result) + 1,"method"] <- "H-SMUCE"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    #### PELT ####
    time <- system.time({fit <- cpt.meanvar(cp_data$y, method = "PELT")})[3] # fit model and time process
    L_est <- length(fit@cpts) - 1
    est_cp <- fit@cpts[-(L_est+1)] + 1
    mean_mse <- mean((rep(fit@param.est$mean, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
    var_mse <-  mean((rep(fit@param.est$variance, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
    result[nrow(result) + 1,"method"] <- "PELT"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, NA, NA, NA, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    #### MOSUM ####
    
    # bottom up merging
    time <- system.time({fit = multiscale.bottomUp(cp_data$y, do.confint = TRUE, level = 0.1, var.est.method = "mosum.min")})[3] # fit model and time process
    est_cp <- fit$cpts + 1
    L_est <- length(est_cp)
    mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
    vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
    mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
    var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(cs, function(x) any(detected %in% x)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MOSUM BU"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    # local pruning 
    time <- system.time({fit = multiscale.localPrune(cp_data$y, do.confint = TRUE, level = 0.1, var.est.method = "mosum.min")})[3] # fit model and time process
    est_cp <- fit$cpts + 1
    L_est <- length(est_cp)
    mn <- sapply(1:(L_est + 1), function(i) mean(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
    vr <- sapply(1:(L_est + 1), function(i) var(cp_data$y[c(1, est_cp, T + 1)[i]:(c(1, est_cp, T + 1)[i+1]-1)]))
    mean_mse <- mean((rep(mn, diff(c(1, est_cp, T+1))) - cp_data$mean_signal)^2)
    var_mse <- mean((rep(vr, diff(c(1, est_cp, T+1))) - cp_data$var_signal)^2)
    cs <- lapply(1:L_est, function(i) (fit$ci$CI$unif.left[i] + 1):(fit$ci$CI$unif.right[i] + 1))
    detected <- true_cp[true_cp %in% unlist(lapply(est_cp, function(x) (x - window):(x + window)))] # detected true change-points
    n_detected <- length(detected)
    n_covered <- ifelse(L_est == 0, 0, sum(sapply(cs, function(x) any(detected %in% x)))) # number of CS that cover a detected change-point
    avg_len <- ifelse(L_est == 0, NA, mean(sapply(cs, function(x) length(x)))) # average CS length
    result[nrow(result) + 1,"method"] <- "MOSUM LP"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, n_detected, n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
    
    #### NOT ####
    time <- system.time({fit = not(cp_data$y, contrast = "pcwsConstMeanVar")})[3] # fit model and time process
    est_cp <- features(fit)$cpt + 1
    L_est <- length(est_cp)
    mean_mse <- mean((predict(fit)[,1] - cp_data$mean_signal)^2)
    var_mse <- mean((predict(fit)[,2]^2  - cp_data$var_signal)^2)
    result[nrow(result) + 1,"method"] <- "NOT"
    result[nrow(result),-1] <- c(T, L, min_space, time, L_est, NA, NA, NA, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), fpsle(true_cp, est_cp, T), fnsle(true_cp, est_cp, T))
  }
  result
}

saveRDS(results, "Desktop/results.rds")

results <- readRDS("Desktop/results.rds")
i=1

results %>% 
  filter(T == settings$T[i], L == settings$L[i], min_space == settings$min_space[i]) %>% 
  group_by(method, T, L, min_space) %>% 
  summarize("|L - L_hat|" = mean(abs(L - L_est)),
            "<= -2" = mean((L - L_est) <= -2),
            "= -1" = mean((L - L_est) == -1),
            "= 0" = mean((L - L_est) == 0),
            "= 1" = mean((L - L_est) == 1),
            ">= 2" = mean((L - L_est) >= 2),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est),
            coverage = sum(n_covered, na.rm = TRUE) / sum(n_detected, na.rm = TRUE),
            hausdorff = mean(hausdorff),
            fpsle = mean(fpsle),
            fnsle = mean(fpsle),
            mean_mse = mean(mean_mse),
            var_mse = mean(var_mse),
            time = mean(time))

            
            
            

