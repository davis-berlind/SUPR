library(foreach)
library(doParallel)
library(dplyr)

# simulation settings
settings <- data.frame(T = c(rep(100, 2), rep(500, 4), rep(1000, 4)),
                       L = c(2, 5, 2, 2, 5, 5, 2, 2, 10, 10),
                       min_space = c(15, 15, 15, 30, 15, 30, 30, 50, 30, 50))

cols <- c("method", "merge_prob", "T", "L", "min_space", "time", "L_est", 
          "n_detected", "n_covered", "avg_len", "mean_mse", "var_mse", 
          "hausdorff", "fpsle", "fnsle")

merge_probs <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.6, 0.7, 0.8, 0.9)
pi_100 = rowMeans(sapply(1:1e5, function(i) smscp(rnorm(100), 
                                                   lambda = rep(1,100),
                                                   tau = 0.001, 
                                                   u = 0.001, 
                                                   v = rep(0.001,100), 
                                                   log_pi = rep(0,100), 
                                                   B_r = 0)$pi))
pi_100 = (1 / pi_100) / sum(1 / pi_100)
pi_500 = rowMeans(sapply(1:1e5, function(i) smscp(rnorm(500), 
                                                   lambda = rep(1, 500),
                                                   tau = 0.001, 
                                                   u = 0.001, 
                                                   v = rep(0.001, 500), 
                                                   log_pi = rep(0,500), 
                                                   B_r = 0)$pi))
pi_500 = (1 / pi_500) / sum(1 / pi_500)
pi_1000 = rowMeans(sapply(1:1e5, function(i) smscp(rnorm(1000), 
                                                    lambda = rep(1,1000),
                                                    tau = 0.001, 
                                                    u = 0.001, 
                                                    v = rep(0.001,1000), 
                                                    log_pi = rep(0,1000), 
                                                    B_r = 0)$pi))
pi_1000 = (1 / pi_1000) / sum(1 / pi_1000)

registerDoParallel(8)

N <- 1000 # number of simulation repetitions
tol <- 1e-3 # convergence criterion

#### Simulation ####

results <- foreach(i = 1:N, .combine = rbind, .verbose = TRUE) %dopar% {
  # initialize results matrix
  result <- data.frame(character(0), sapply(1:(length(cols)-1), function(i) numeric(0)))
  names(result) <- cols
  
  for (j in 1:nrow(settings)) {
    # extract settings
    T <- settings$T[j]; L <- settings$L[j]; min_space <- settings$min_space[j];
    
    # conditional detection window
    window <- floor(min(sqrt(T) / 2, 15)) 
    
    # generate data
    cp_data <- hsmuce_simulation(T, L, 200, min_space)
    true_cp <- cp_data$changepoints
    
    #### MICH Models ####
    local_rate <- round(T / sqrt(T * log(T))) 
    
    #### MICH (J,0,0) Localization Rate ####
    if (T == 100) pi_j <- sapply(1:local_rate, function(i) pi_100)
    else if (T == 500) pi_j <- sapply(1:local_rate, function(i) pi_500)
    else pi_j <- sapply(1:local_rate, function(i) pi_1000)
    
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J = local_rate, pi_j = pi_j, 
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    for (mp in merge_probs) {
      # determine MICH credible sets
      cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = mp)
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
      result[nrow(result) + 1,"method"] <- "MICH (J,0,0) LR"
      result[nrow(result),-1] <- c(mp, T, L, min_space, time, L_est, n_detected, 
                                   n_covered, avg_len, mean_mse, var_mse,
                                   hausdorff(true_cp, est_cp, T), 
                                   fpsle(true_cp, est_cp, T), 
                                   fnsle(true_cp, est_cp, T))
    }
    
    #### MICH (0,J) Oracle ####
    if (T == 100) pi_j <- sapply(1:L, function(i) pi_100)
    else if (T == 500) pi_j <- sapply(1:L, function(i) pi_500)
    else pi_j <- sapply(1:L, function(i) pi_1000)
    # fit model and time process
    time <- system.time({fit <- mich(cp_data$y, J = L, pi_j = pi_j, 
                                     tol = tol, max_iter = Inf)})[3] 
    
    # calculate mean/precision signal MSEs
    mean_mse <- mean((fit$mu - cp_data$mean_signal)^2)
    var_mse <- mean((1/fit$lambda - cp_data$var_signal)^2)
    
    for (mp in merge_probs) {
      # determine MICH credible sets
      cs <- mich_sets(fit$J_model$pi, max_length = T %/% 4, merge_prob = mp)
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
      result[nrow(result) + 1,"method"] <- "MICH (J,0,0) Oracle"
      result[nrow(result),-1] <- c(mp, T, L, min_space, time, L_est, n_detected, 
                                 n_covered, avg_len, mean_mse, var_mse,
                                 hausdorff(true_cp, est_cp, T), 
                                 fpsle(true_cp, est_cp, T), 
                                 fnsle(true_cp, est_cp, T))
    }
    
  }
  result
}

saveRDS(results, "~/Desktop/merge_prob_results.rds")

results <- readRDS("Desktop/merge_prob_results.rds")


i=10
results %>% 
  filter(method == "MICH (J,0,0) LR", T == settings$T[i], L == settings$L[i], min_space == settings$min_space[i]) %>% 
  group_by(method, merge_prob, T, L, min_space) %>% 
  summarize("|L - L_hat|" = mean(abs(L - L_est)),
            "<= -2" = mean((L - L_est) <= -2),
            "= -1" = mean((L - L_est) == -1),
            "= 0" = mean((L - L_est) == 0),
            "= 1" = mean((L - L_est) == 1),
            ">= 2" = mean((L - L_est) >= 2),
            ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est),
            coverage = sum(n_covered, na.rm = TRUE) / sum(n_detected, na.rm = TRUE),
            hausdorff = mean(hausdorff, na.rm = TRUE),
            # fpsle = mean(fpsle, na.rm = TRUE),
            # fnsle = mean(fnsle, na.rm = TRUE),
            #mean_mse = mean(mean_mse, na.rm = TRUE),
            #var_mse = mean(var_mse, na.rm = TRUE),
            time = mean(time, na.rm = TRUE)) %>% 
  mutate_if(is.numeric, ~format(round(.,3), nsmall = 3)) %>% 
  arrange(T, L, min_space) %>% print(n = 120) 
# write.csv("~/Desktop/results.csv", row.names = FALSE)





