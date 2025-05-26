mich_vector <- function(y, fit_intercept, fit_scale,
                        J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
                        pi_j_weighted, pi_l_weighted, pi_k_weighted,
                        tol, verbose, max_iter, reverse,
                        detect, merge_level, merge_prob,
                        restart, n_restart, n_search, increment,
                        omega_j, u_j, v_j, pi_j,
                        omega_l, pi_l,
                        u_k, v_k, pi_k) {
  
  #### set up ####
  # Flag for new model 
  refit <- FALSE
  
  # calculate dimensions of y
  T <- length(y)
  
  merge_counter = log(T) %/% 2
  
  #### initialize mu_0 ####
  if (fit_intercept) {
    mu_0 <- mean(y[1:ceiling(log(T))])
  } else mu_0 <- 0.0
  
  #### initialize lambda_0 ####
  if (fit_scale) {
    lambda_0 <- 1 / var(y[1:ceiling(log(T))])
  } else lambda_0 <- 1.0
  
  #### initializing posterior parameters ####
  
  # mean components
  log_pi_l <- log(pi_l)
  pi_bar_l <- matrix(1/T, nrow = T, ncol = L)
  log_pi_bar_l <- matrix(0, nrow = T, ncol = L)
  b_bar_l <- matrix(0.0, nrow = T, ncol = L)
  omega_bar_l <- matrix(1.0, nrow = T, ncol = L)
  
  # variance components
  log_pi_k <- log(pi_k)
  pi_bar_k <- matrix(1/T, nrow = T, ncol = K)
  log_pi_bar_k <- matrix(0, nrow = T, ncol = K)
  u_bar_k <- u_k + (T-1:T+1) / 2
  lgamma_u_bar_k <- lgamma(u_bar_k)
  digamma_u_bar_k <- digamma(u_bar_k)
  v_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T-1:T+1) / 2
  
  # mean-variance components
  log_pi_j <- log(pi_j)
  pi_bar_j <- matrix(1/T, nrow = T, ncol = J)
  log_pi_bar_j <- matrix(0.0, nrow = T, ncol = J)
  b_bar_j <- matrix(0.0, nrow = T, ncol = J)
  omega_bar_j <- matrix(1.0, nrow = T, ncol = J)
  u_bar_j <- u_j + (T-1:T+1) / 2
  lgamma_u_bar_j <- lgamma(u_bar_j)
  digamma_u_bar_j <- digamma(u_bar_j)
  
  if (((J_auto & !L_auto & !K_auto) | J > 0) & (L + K == 0)) {
    if (verbose) print("Initializing mean components.")
    if (J_auto) pi_l <- sapply(1:max(1,J), function(i) pi_l[,1])
    else pi_l <- pi_j
    
    fit <- mich_vector(y, fit_intercept, fit_scale,
                       J = 0, L = J, K, 
                       J_auto = L_auto, L_auto = J_auto, K_auto,
                       J_max = 0, L_max = J_max, K_max,
                       pi_j_weighted, pi_l_weighted, pi_k_weighted,
                       tol, verbose, max_iter, reverse = FALSE,
                       detect, merge_level, merge_prob,
                       restart=FALSE, n_restart, n_search %/% 2, increment,
                       omega_j, u_j, v_j, pi_j,
                       omega_l, pi_l,
                       u_k, v_k, pi_k)
    
    J <- fit$L
    if (J > 0) {
      keep <- seq_len(J)
      if (J_auto) {
        cred_sets <- apply(fit$mean_model$pi_bar, 2, cred_set, level = merge_level, simplify = FALSE)
        keep <- sapply(cred_sets, length) <= detect
        J = sum(keep)
        pi_j <- sapply(1:max(1,J), function(i) pi_j[,1])
        log_pi_j <- sapply(1:max(1,J), function(i) log_pi_j[,1])
      }
      pi_bar_j <- fit$mean_model$pi_bar[, keep, drop = FALSE]
      log_pi_bar_j <- log(pi_bar_j)
      b_bar_j <- fit$mean_model$b_bar[, keep, drop = FALSE]
      omega_bar_j <- fit$mean_model$omega_bar[, keep, drop = FALSE]
      refit <- TRUE
    }
  }
  v_bar_j <- sapply(1:max(1,J), function(i) u_bar_j)
  
  #### fit model and merge ####
  merged <- FALSE
  merge_attempts <- 0
  while (!merged) {
    fit <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                    fit_intercept, fit_scale, refit,
                    max_iter, verbose = (verbose & sum(J_auto, L_auto, K_auto) == 0), tol, 
                    omega_j, u_j, v_j, log_pi_j, 
                    pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, 
                    u_bar_j, v_bar_j, lgamma_u_bar_j, digamma_u_bar_j,
                    omega_l, log_pi_l, pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                    u_k, v_k, log_pi_k, pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                    lgamma_u_bar_k, digamma_u_bar_k)
    
    merged <- TRUE
    refit <- TRUE
    
    if (L > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE
      
      # compute pairwise merge probabilities 
      merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
      diag(merge_prob_mat) <- 0
      
      mu_bar_l <- fit$mean_model$mu_bar
      mu2_bar_l <- fit$mean_model$mu2_bar
      merge_residual <- fit$residual
      merge_lambda <- fit$lambda
      merge_delta <- fit$delta
      
      while (L > 1 & any(merge_prob_mat[keep_mat] > (merge_prob * 1^merge_attempts))) {
        merged <- FALSE
        L <- L - 1
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]

        mu_bar_merge <- rowSums(mu_bar_l[,merge_dex])
        merge_residual <- merge_residual + mu_bar_merge
        merge_delta <- merge_delta + rowSums(mu_bar_l[,merge_dex]^2 - mu2_bar_l[,merge_dex]) 
        
        merge_fit <- mean_scp(merge_residual, merge_lambda, omega_l, log_pi_l[,merge_dex[1]])
        
        pi_bar_l[,merge_dex[1]] <- merge_fit$pi_bar
        log_pi_bar_l[,merge_dex[1]] <- merge_fit$log_pi_bar
        b_bar_l[,merge_dex[1]] <- merge_fit$b_bar
        omega_bar_l[,merge_dex[1]] <- merge_fit$omega_bar

        mu_bar_l[,merge_dex[1]] <- mu_bar_fn(merge_fit$b_bar, merge_fit$pi_bar)
        mu2_bar_l[,merge_dex[1]] <- mu2_bar_fn(merge_fit$b_bar, merge_fit$omega_bar,  merge_fit$pi_bar)

        merge_residual <- merge_residual - mu_bar_l[,merge_dex[1]]
        merge_delta <- merge_delta - mu_bar_l[,merge_dex[1]]^2 + mu2_bar_l[,merge_dex[1]] 

        if (length(cred_set(pi_bar_l[,merge_dex[1]], level = merge_level)) > detect) {
          keep_mat[merge_dex[1], ] <- FALSE
          keep_mat[, merge_dex[1]] <- FALSE
        }
        
        merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_l[,merge_dex[1]] %*% pi_bar_l[,-merge_dex])
        merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
        merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
        keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
        
        pi_bar_l <- pi_bar_l[,-merge_dex[2], drop=FALSE] 
        log_pi_bar_l <- log_pi_bar_l[,-merge_dex[2], drop=FALSE]
        b_bar_l <- b_bar_l[,-merge_dex[2], drop=FALSE]
        omega_bar_l <- omega_bar_l[,-merge_dex[2], drop=FALSE]
        if (L_auto) {
          pi_l <- pi_l[,-merge_dex[2], drop=FALSE]
          log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
        }        
        mu_bar_l <- mu_bar_l[,-merge_dex[2], drop=FALSE]
        mu2_bar_l <- mu2_bar_l[,-merge_dex[2], drop=FALSE]
      }
    }
    
    if (K > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = K, nrow = K)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE
      
      # compute pairwise merge probabilities 
      merge_prob_mat <- t(pi_bar_k) %*% pi_bar_k
      diag(merge_prob_mat) <- 0
      
      lambda_bar_k <- fit$var_model$lambda_bar
      merge_residual <- fit$residual
      merge_lambda <- fit$lambda
      merge_delta <- fit$delta
      
      while (K > 1 & any(merge_prob_mat[keep_mat] > (merge_prob * 1^merge_attempts))) {
        merged <- FALSE
        K <- K - 1
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
        
        merge_lambda <- merge_lambda / apply(lambda_bar_k[,merge_dex], 1, prod)

        v_merge <- v_k + revcumsum(0.5 * merge_lambda * merge_delta)
        log_pi_merge <- log_pi_k[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))
        
        merge_fit <- var_scp(merge_residual, merge_lambda, u_bar_k, lgamma_u_bar_k, 
                             v_merge, log_pi_merge - max(log_pi_merge))
        
        pi_bar_k[,merge_dex[1]] <- merge_fit$pi_bar
        log_pi_bar_l[,merge_dex[1]] <- merge_fit$log_pi_bar
        v_bar_k[,merge_dex[1]] <- merge_fit$v_bar
        
        lambda_bar_k[,merge_dex[1]] <- lambda_bar_fn(u_bar_k, merge_fit$v_bar,  merge_fit$pi_bar)
        
        merge_lambda <- merge_lambda * lambda_bar_k[,merge_dex[1]]

        if (length(cred_set(pi_bar_k[,merge_dex[1]], level = merge_level)) > detect) {
          keep_mat[merge_dex[1], ] <- FALSE
          keep_mat[, merge_dex[1]] <- FALSE
        }
        
        merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_k[,merge_dex[1]] %*% pi_bar_k[,-merge_dex])
        merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
        merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
        keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
        
        pi_bar_k <- pi_bar_k[,-merge_dex[2], drop=FALSE] 
        log_pi_bar_k <- log_pi_bar_k[,-merge_dex[2], drop=FALSE]
        v_bar_k <- v_bar_k[, -merge_dex[2], drop=FALSE]
        if (K_auto) {
          pi_k <- pi_k[,-merge_dex[2], drop=FALSE]
          log_pi_k <- log_pi_k[,-merge_dex[2], drop=FALSE]
        }
        lambda_bar_k <- lambda_bar_k[,-merge_dex[2], drop=FALSE]
      }
    }
    
    if (J > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = J, nrow = J)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE
      
      # compute pairwise merge probabilities 
      merge_prob_mat <- t(pi_bar_j) %*% pi_bar_j
      diag(merge_prob_mat) <- 0
      
      mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar
      mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar
      lambda_bar_j <- fit$meanvar_model$lambda_bar
      merge_residual <- fit$residual
      merge_lambda <- fit$lambda
      merge_delta <- fit$delta
      
      while (J > 1 & any(merge_prob_mat[keep_mat] > (merge_prob * 1^merge_attempts))) {
        merged <- FALSE
        J <- J - 1
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
        
        mu_bar_merge <- rowSums(mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
        merge_residual <- merge_residual + mu_bar_merge
        merge_lambda <- merge_lambda / apply(lambda_bar_j[,merge_dex], 1, prod)
        merge_delta <- merge_delta + rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2) 
        merge_delta <- merge_delta - rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
        
        v_merge <- v_j + revcumsum(0.5 * merge_lambda * merge_delta)
        log_pi_merge <- log_pi_j[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))
        
        merge_fit <- meanvar_scp(merge_residual, merge_lambda, omega_j, u_bar_j, 
                                 lgamma_u_bar_j, v_merge, log_pi_merge - max(log_pi_merge))
        
        pi_bar_j[,merge_dex[1]] <- merge_fit$pi_bar
        log_pi_bar_j[,merge_dex[1]] <- merge_fit$log_pi_bar
        b_bar_j[,merge_dex[1]] <- merge_fit$b_bar
        omega_bar_j[,merge_dex[1]] <- merge_fit$omega_bar
        v_bar_j[,merge_dex[1]] <- merge_fit$v_bar
        
        mu_lambda_bar_j[,merge_dex[1]] <- mu_lambda_fn(merge_fit$b_bar, u_bar_j,  merge_fit$v_bar,  merge_fit$pi_bar)
        mu2_lambda_bar_j[,merge_dex[1]] <- mu2_lambda_fn(merge_fit$b_bar, merge_fit$omega_bar, u_bar_j,  merge_fit$v_bar,  merge_fit$pi_bar)
        lambda_bar_j[,merge_dex[1]] <- lambda_bar_fn(u_bar_j, merge_fit$v_bar,  merge_fit$pi_bar)
        
        mu_bar_merge <- mu_lambda_bar_j[,merge_dex[1]] / lambda_bar_j[,merge_dex[1]]
        merge_residual <- merge_residual - mu_bar_merge
        merge_lambda <- merge_lambda * lambda_bar_j[,merge_dex[1]]
        merge_delta <- merge_delta - rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2) 
        merge_delta <- merge_delta + rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
        
        if (length(cred_set(pi_bar_j[,merge_dex[1]], level = merge_level)) > detect) {
          keep_mat[merge_dex[1], ] <- FALSE
          keep_mat[, merge_dex[1]] <- FALSE
        }
        
        merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_j[,merge_dex[1]] %*% pi_bar_j[,-merge_dex])
        merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
        merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
        keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
        
        pi_bar_j <- pi_bar_j[,-merge_dex[2], drop=FALSE] 
        log_pi_bar_j <- log_pi_bar_j[,-merge_dex[2], drop=FALSE]
        b_bar_j <- b_bar_j[,-merge_dex[2], drop=FALSE]
        omega_bar_j <- omega_bar_j[,-merge_dex[2], drop=FALSE]
        v_bar_j <- v_bar_j[, -merge_dex[2], drop=FALSE]
        if (J_auto) {
          pi_j <- pi_j[,-merge_dex[2], drop=FALSE]
          log_pi_j <- log_pi_j[,-merge_dex[2], drop=FALSE]
        }
        mu_lambda_bar_j <- mu_lambda_bar_j[,-merge_dex[2], drop=FALSE]
        mu2_lambda_bar_j <- mu2_lambda_bar_j[,-merge_dex[2], drop=FALSE]
        lambda_bar_j <- lambda_bar_j[,-merge_dex[2], drop=FALSE]
      }
    }
    merge_attempts <- merge_attempts + 1
    if (verbose & !merged) print(paste0("Merging components. Attempt: ", merge_attempts))
  }
  
  # if components were merged out use auto procedure to increase to desired JLK
  merge_flag <- (J < J_max & !J_auto) | (L < L_max & !L_auto) | (K < K_max & !K_auto)
  merge_elbo <- -Inf
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components
  
  #### auto procedures ####
  #### auto procedure with single component ####
  last_restart <- ifelse(restart, 2, Inf)
  
  if (sum(J_auto, L_auto, K_auto) == 1 | merge_flag) {
    
    refit <- (L > 0 | K > 0 | J > 0)
    counter <- n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo

    if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", 
                              elbo_new, "; Counter: ", counter))
    
    # continue search until n_search exhausted or max components exceeded
    while ((J < J_max | L < L_max | K < K_max) & counter > 0) {
      if (J_auto | J < J_max) {
        # increment dimension of parameters
        J <- J + increment
        if (J > 1 & J_auto) {
          pi_j <- cbind(matrix(pi_j[,1], nrow = T, ncol = increment), pi_j)
          log_pi_j <- cbind(matrix(log_pi_j[,1], nrow = T, ncol = increment), log_pi_j)
        }
        pi_bar_j <- cbind(matrix(1/T, nrow = T, ncol = increment), pi_bar_j)
        log_pi_bar_j <- cbind(matrix(0.0, nrow = T, ncol = increment), log_pi_bar_j)
        b_bar_j <- cbind(matrix(0.0, nrow = T, ncol = increment), b_bar_j)
        omega_bar_j <- cbind(matrix(1.0, nrow = T, ncol = increment), omega_bar_j)
        v_bar_j <- cbind(matrix(v_j + (T-1:T+1) / 2, nrow = T, ncol = increment), v_bar_j)
      }
      
      if (L_auto | L < L_max) {
        # increment dimension of parameters
        L <- L + increment
        if (L > 1 & L_auto) {
          pi_l <- cbind(matrix(pi_l[,1], nrow = T, ncol = increment), pi_l)
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
        }
        pi_bar_l <- cbind(matrix(1/T, nrow = T, ncol = increment), pi_bar_l)
        log_pi_bar_l <- cbind(matrix(0.0, nrow = T, ncol = increment), log_pi_bar_l)
        b_bar_l <- cbind(matrix(0.0, nrow = T, ncol = increment), b_bar_l)
        omega_bar_l <- cbind(matrix(1.0, nrow = T, ncol = increment), omega_bar_l)
      }
      
      if (K_auto | K < K_max) {
        # increment dimension of parameters
        K <- K + increment
        if (K > 1 & K_auto) {
          pi_k <- cbind(matrix(pi_k[,1], nrow = T, ncol = increment), pi_k)
          log_pi_k <- cbind(matrix(log_pi_k[,1], nrow = T, ncol = increment), log_pi_k)
        }
        pi_bar_k <- cbind(matrix(1/T, nrow = T, ncol = increment), pi_bar_k)
        log_pi_bar_k <- cbind(matrix(0.0, nrow = T, ncol = increment), log_pi_bar_k)
        v_bar_k <- cbind(matrix(v_k + (T-1:T+1) / 2, nrow = T, ncol = increment), v_bar_k)
      }
      
      fit_new <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                          fit_intercept, fit_scale, refit,
                          max_iter = max_iter, 
                          verbose = FALSE, tol, 
                          omega_j, u_j, v_j, log_pi_j, 
                          pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, 
                          u_bar_j, v_bar_j, lgamma_u_bar_j, digamma_u_bar_j,
                          omega_l, log_pi_l, pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                          u_k, v_k, log_pi_k, pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                          lgamma_u_bar_k, digamma_u_bar_k)

      # test if model improved or restart ####
      if (!refit) refit <- TRUE
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J,
                                "): ELBO = ", elbo_new, "; Counter: ", counter))
      
      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- n_search
        if (last_restart < Inf) restart <- TRUE 
      } else if (restart & J_auto & J - n_restart > last_restart) {  
        refit <- FALSE
        restart <- FALSE
        last_restart <- J
        n_restart <- n_restart + 1
        
        J <- J - increment
        pi_j <- sapply(1:J, function(i) pi_j[,1])
        log_pi_j <- sapply(1:J, function(i) log_pi_j[,1])
        
        fit_new <- mich_vector(y, fit_intercept, fit_scale,
                               J = 0, L = J, K, 
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               J_max = 0, L_max = J, K_max,
                               pi_j_weighted, pi_l_weighted, pi_k_weighted,
                               tol, verbose = FALSE, max_iter, reverse = FALSE,
                               detect, merge_level, merge_prob,
                               restart = FALSE, n_restart, n_search, increment,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l = pi_j,
                               u_k, v_k, pi_k)
        
        refit <- TRUE
        pi_bar_j <- fit_new$mean_model$pi_bar
        log_pi_bar_j <-log(pi_bar_j)
        b_bar_j <- fit_new$mean_model$b_bar
        omega_bar_j <- fit_new$mean_model$omega_bar
        v_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T-1:T+1) / 2
      } else if (restart & L_auto & L - n_restart > last_restart) {
        refit <- FALSE
        restart <- FALSE
        last_restart <- L
        n_restart <- n_restart + 1
        
        L <- L - increment
        pi_l <- sapply(1:L, function(i) pi_l[,1])
        log_pi_l <- sapply(1:L, function(i) log_pi_l[,1])
        pi_bar_l <- matrix(1/T, nrow = T, ncol = L)
        log_pi_bar_l <- matrix(0, nrow = T, ncol = L)
        b_bar_l <- matrix(0.0, nrow = T, ncol = L)
        omega_bar_l <- matrix(1.0, nrow = T, ncol = L)
      } else if (restart & K_auto & K - n_restart > last_restart) {
        refit <- FALSE
        restart <- FALSE
        last_restart <- K
        n_restart <- n_restart + 1
        
        K <- K - increment
        pi_k <- sapply(1:K, function(i) pi_k[,1])
        log_pi_k <- sapply(1:K, function(i) log_pi_k[,1])
        pi_bar_k <- matrix(1/T, nrow = T, ncol = K)
        log_pi_bar_k <- matrix(0, nrow = T, ncol = K)
        v_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T-1:T+1) / 2
      } else if (counter == 1 | (L == L_max & K == K_max & J == J_max)) {
        # merging ####
        merge_attempts <- 0
        while (TRUE) {
          merged <- TRUE  
          L <- fit$L; K <- fit$K; J <- fit$J;
          if (verbose) print(paste0("Merging at (L = ", L, ", K = ", K, ", J = ", J,
                                    "): Merge Counter: ", merge_counter))
          
          if (L > 1) {
            pi_bar_l <- fit$mean_model$pi_bar
            log_pi_bar_l <- log(pi_bar_l)
            b_bar_l <- fit$mean_model$b_bar
            omega_bar_l <- fit$mean_model$omega_bar
            
            # test which columns actually contain change-points
            cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = L, nrow = L)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE
            
            # compute pairwise merge probabilities 
            merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
            diag(merge_prob_mat) <- 0
            
            mu_bar_l <- fit$mean_model$mu_bar
            mu2_bar_l <- fit$mean_model$mu2_bar
            merge_residual <- fit$residual
            merge_lambda <- fit$lambda
            merge_delta <- fit$delta
            
            while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob * (1^merge_attempts))) {
              merged <- FALSE
              L <- L - 1
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
              
              mu_bar_merge <- rowSums(mu_bar_l[,merge_dex])
              merge_residual <- merge_residual + mu_bar_merge
              merge_delta <- merge_delta + rowSums(mu_bar_l[,merge_dex]^2 - mu2_bar_l[,merge_dex]) 
              
              merge_fit <- mean_scp(merge_residual, merge_lambda, omega_l, log_pi_l[,merge_dex[1]])
              
              pi_bar_l[,merge_dex[1]] <- merge_fit$pi_bar
              log_pi_bar_l[,merge_dex[1]] <- merge_fit$log_pi_bar
              b_bar_l[,merge_dex[1]] <- merge_fit$b_bar
              omega_bar_l[,merge_dex[1]] <- merge_fit$omega_bar
              
              mu_bar_l[,merge_dex[1]] <- mu_bar_fn(merge_fit$b_bar, merge_fit$pi_bar)
              mu2_bar_l[,merge_dex[1]] <- mu2_bar_fn(merge_fit$b_bar, merge_fit$omega_bar,  merge_fit$pi_bar)
              
              merge_residual <- merge_residual - mu_bar_l[,merge_dex[1]]
              merge_delta <- merge_delta - mu_bar_l[,merge_dex[1]]^2 + mu2_bar_l[,merge_dex[1]]
              
              if (length(cred_set(pi_bar_l[,merge_dex[1]], level = merge_level)) > detect) {
                keep_mat[merge_dex[1], ] <- FALSE
                keep_mat[, merge_dex[1]] <- FALSE
              }
              
              merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_l[,merge_dex[1]] %*% pi_bar_l[,-merge_dex])
              merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
              merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
              keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
              
              pi_bar_l <- pi_bar_l[,-merge_dex[2], drop=FALSE] 
              log_pi_bar_l <- log_pi_bar_l[,-merge_dex[2], drop=FALSE]
              b_bar_l <- b_bar_l[,-merge_dex[2], drop=FALSE]
              omega_bar_l <- omega_bar_l[,-merge_dex[2], drop=FALSE]
              if (L_auto) {
                pi_l <- pi_l[,-merge_dex[2], drop=FALSE]
                log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
              }            
              mu_bar_l <- mu_bar_l[,-merge_dex[2], drop=FALSE]
              mu2_bar_l <- mu2_bar_l[,-merge_dex[2], drop=FALSE]
            }
          }
          
          if (K > 1) {
            pi_bar_k <- fit$var_model$pi_bar
            log_pi_bar_k <- log(pi_bar_k)
            v_bar_k <- fit$mean_model$v_bar
            
            # test which columns actually contain change-points
            cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = K, nrow = K)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE
            
            # compute pairwise merge probabilities 
            merge_prob_mat <- t(pi_bar_k) %*% pi_bar_k
            diag(merge_prob_mat) <- 0
            
            lambda_bar_k <- fit$var_model$lambda_bar
            merge_residual <- fit$residual
            merge_lambda <- fit$lambda
            merge_delta <- fit$delta
            
            while (K > 1 & any(merge_prob_mat[keep_mat] > merge_prob *(1^merge_attempts))) {
              merge <- FALSE
              K <- K - 1
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
              
              merge_lambda <- merge_lambda / apply(lambda_bar_k[,merge_dex], 1, prod)
              
              v_merge <- v_k + revcumsum(0.5 * merge_lambda * merge_delta)
              log_pi_merge <- log_pi_k[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))
              
              merge_fit <- var_scp(merge_residual, merge_lambda, u_bar_k, lgamma_u_bar_k, 
                                   v_merge, log_pi_merge - max(log_pi_merge))
              
              pi_bar_k[,merge_dex[1]] <- merge_fit$pi_bar
              log_pi_bar_l[,merge_dex[1]] <- merge_fit$log_pi_bar
              v_bar_k[,merge_dex[1]] <- merge_fit$v_bar
              
              lambda_bar_k[,merge_dex[1]] <- lambda_bar_fn(u_bar_k, merge_fit$v_bar,  merge_fit$pi_bar)
              
              merge_lambda <- merge_lambda * lambda_bar_k[,merge_dex[1]]
              
              if (length(cred_set(pi_bar_k[,merge_dex[1]], level = merge_level)) > detect) {
                keep_mat[merge_dex[1], ] <- FALSE
                keep_mat[, merge_dex[1]] <- FALSE
              }
              
              merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_k[,merge_dex[1]] %*% pi_bar_k[,-merge_dex])
              merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
              merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
              keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
              
              pi_bar_k <- pi_bar_k[,-merge_dex[2], drop=FALSE] 
              log_pi_bar_k <- log_pi_bar_k[,-merge_dex[2], drop=FALSE]
              v_bar_k <- v_bar_k[, -merge_dex[2], drop=FALSE]
              if (K_auto) {
                pi_k <- pi_k[,-merge_dex[2], drop=FALSE]
                log_pi_k <- log_pi_k[,-merge_dex[2], drop=FALSE]
              }            
              lambda_bar_k <- lambda_bar_k[,-merge_dex[2], drop=FALSE]
            }
          }
          
          if (J > 1) {
            pi_bar_j <- fit$meanvar_model$pi_bar
            log_pi_bar_j <- log(pi_bar_j)
            b_bar_j <- fit$meanvar_model$b_bar
            omega_bar_j <- fit$meanvar_model$omega_bar
            v_bar_j <- fit$meanvar_model$v_bar
            
            # test which columns actually contain change-points
            cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
            keep <- sapply(cred_sets, length) <= detect
            keep_mat <- matrix(FALSE, ncol = J, nrow = J)
            keep_mat[keep, keep] <- TRUE
            diag(keep_mat) <- FALSE
            
            # compute pairwise merge probabilities 
            merge_prob_mat <- t(pi_bar_j) %*% pi_bar_j
            diag(merge_prob_mat) <- 0
            
            mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar
            mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar
            lambda_bar_j <- fit$meanvar_model$lambda_bar
            merge_residual <- fit$residual
            merge_lambda <- fit$lambda
            merge_delta <- fit$delta
            
            while (J > 1 & any(merge_prob_mat[keep_mat] > merge_prob * (1^merge_attempts))) {
              merged <- FALSE
              J <- J - 1
              merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
              
              mu_bar_merge <- rowSums(mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
              merge_residual <- merge_residual + mu_bar_merge
              merge_lambda <- merge_lambda / apply(lambda_bar_j[,merge_dex], 1, prod)
              merge_delta <- merge_delta + rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2) 
              merge_delta <- merge_delta - rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
              
              v_merge <- v_j + revcumsum(0.5 * merge_lambda * merge_delta)
              log_pi_merge <- log_pi_j[,merge_dex[1]] - cumsum(c(0,0.5 * merge_lambda[-T] * merge_delta[-T]))
              
              merge_fit <- meanvar_scp(merge_residual, merge_lambda, omega_j, u_bar_j, 
                                       lgamma_u_bar_j, v_merge, log_pi_merge - max(log_pi_merge))
              
              pi_bar_j[,merge_dex[1]] <- merge_fit$pi_bar
              log_pi_bar_j[,merge_dex[1]] <- merge_fit$log_pi_bar
              b_bar_j[,merge_dex[1]] <- merge_fit$b_bar
              omega_bar_j[,merge_dex[1]] <- merge_fit$omega_bar
              v_bar_j[,merge_dex[1]] <- merge_fit$v_bar
              
              mu_lambda_bar_j[,merge_dex[1]] <- mu_lambda_fn(merge_fit$b_bar, u_bar_j,  merge_fit$v_bar,  merge_fit$pi_bar)
              mu2_lambda_bar_j[,merge_dex[1]] <- mu2_lambda_fn(merge_fit$b_bar, merge_fit$omega_bar, u_bar_j,  merge_fit$v_bar,  merge_fit$pi_bar)
              lambda_bar_j[,merge_dex[1]] <- lambda_bar_fn(u_bar_j, merge_fit$v_bar,  merge_fit$pi_bar)
              
              mu_bar_merge <- mu_lambda_bar_j[,merge_dex[1]] / lambda_bar_j[,merge_dex[1]]
              merge_residual <- merge_residual - mu_bar_merge
              merge_lambda <- merge_lambda * lambda_bar_j[,merge_dex[1]]
              merge_delta <- merge_delta - rowSums((mu_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])^2) 
              merge_delta <- merge_delta + rowSums(mu2_lambda_bar_j[,merge_dex] / lambda_bar_j[,merge_dex])
              
              if (length(cred_set(pi_bar_j[,merge_dex[1]], level = merge_level)) > detect) {
                keep_mat[merge_dex[1], ] <- FALSE
                keep_mat[, merge_dex[1]] <- FALSE
              }
              
              merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_j[,merge_dex[1]] %*% pi_bar_j[,-merge_dex])
              merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
              merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
              keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
              
              pi_bar_j <- pi_bar_j[,-merge_dex[2], drop=FALSE] 
              log_pi_bar_j <- log_pi_bar_j[,-merge_dex[2], drop=FALSE]
              b_bar_j <- b_bar_j[,-merge_dex[2], drop=FALSE]
              omega_bar_j <- omega_bar_j[,-merge_dex[2], drop=FALSE]
              v_bar_j <- v_bar_j[, -merge_dex[2], drop=FALSE]
              if (J_auto) {
                pi_j <- pi_j[,-merge_dex[2], drop=FALSE]
                log_pi_j <- log_pi_j[,-merge_dex[2], drop=FALSE]
              }            
              mu_lambda_bar_j <- mu_lambda_bar_j[,-merge_dex[2], drop=FALSE]
              mu2_lambda_bar_j <- mu2_lambda_bar_j[,-merge_dex[2], drop=FALSE]
              lambda_bar_j <- lambda_bar_j[,-merge_dex[2], drop=FALSE]
            }
          }
          
          if (!merged) {
            fit <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                            fit_intercept, fit_scale, refit,
                            max_iter = max_iter, 
                            verbose = FALSE, tol, 
                            omega_j, u_j, v_j, log_pi_j, 
                            pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, 
                            u_bar_j, v_bar_j, lgamma_u_bar_j, digamma_u_bar_j,
                            omega_l, log_pi_l, pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                            u_k, v_k, log_pi_k, pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                            lgamma_u_bar_k, digamma_u_bar_k)
          } else {
            break
          }
          merge_attempts <- merge_attempts + 1
          if (verbose & !merged) print(paste0("Merging components. Attempt: ", merge_attempts))
        }
        elbo <- max(fit$elbo)
        if (elbo > merge_elbo) {
          merge_elbo <- elbo
          fit_merge <- fit
        }
        counter <- n_search
        merge_counter <- merge_counter - 1
        if (merge_attempts == 0 | merge_counter == 0) {
          fit <- fit_merge
          break
        }
      } else {
        counter <- counter - 1
      }
    }
  }
  
  #### auto procedure with multiple components ####
  if (sum(J_auto, L_auto, K_auto) > 1) {
    
    auto_done <- FALSE
    counter <- (J_auto + L_auto + K_auto) * n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo))
    
    while (TRUE) {
      elbo_new <- rep(-Inf, 3)
      if (J_auto) {
        J <- J + increment
        if (J > 1) {
          pi_j <- cbind(matrix(pi_j[,1], nrow = T, ncol = increment), pi_j)
          log_pi_j <- cbind(matrix(log_pi_j[,1], nrow = T, ncol = increment), log_pi_j)
        }
        
        # fit new model 
        fit_new <- mich_vector(y, fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               pi_j_weighted, pi_l_weighted, pi_k_weighted,
                               tol, verbose = FALSE, max_iter, reverse = FALSE,
                               detect, merge_level, merge_prob,
                               restart, n_restart, n_search, increment,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l,
                               u_k, v_k, pi_k)
        
        elbo_new[1] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new[1]))
        if (elbo_new[1] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[1]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        J <- J - increment
        if (J > 0) {
          pi_j <- sapply(1:J, function(i) pi_j[,1])
          log_pi_j <- sapply(1:J, function(i) log_pi_j[,1])
        }
      }
      
      if (L_auto) {
        L <- L + increment
        if (L > 1) {
          pi_l <- cbind(matrix(pi_l[,1], nrow = T, ncol = increment), pi_l)
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
        }
        
        # fit new model 
        fit_new <- mich_vector(y, fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               pi_j_weighted, pi_l_weighted, pi_k_weighted,
                               tol, verbose = FALSE, max_iter, reverse = FALSE,
                               detect, merge_level, merge_prob,
                               restart, n_restart, n_search, increment,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l,
                               u_k, v_k, pi_k)
        
        elbo_new[2] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new[2]))
        if (elbo_new[2] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[2]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        L <- L - increment
        if (L > 0) {
          pi_l <- sapply(1:L, function(i) pi_l[,1])
          log_pi_l <- sapply(1:L, function(i) log_pi_l[,1])
        }
      }
      
      if (K_auto) {
        K <- K + increment
        if (K > 1) {
          pi_k <- cbind(matrix(pi_k[,1], nrow = T, ncol = increment), pi_k)
          log_pi_k <- cbind(matrix(log_pi_k[,1], nrow = T, ncol = increment), log_pi_k)
        }
        
        # fit new model 
        fit_new <- mich_vector(y, fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               pi_j_weighted, pi_l_weighted, pi_k_weighted,
                               tol, verbose = FALSE, max_iter, reverse = FALSE,
                               detect, merge_level, merge_prob,
                               restart, n_restart, n_search,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l,
                               u_k, v_k, pi_k)
        
        elbo_new[3] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new[3]))
        if (elbo_new[3] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[3]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        K <- K - 1
        if (K > 0) {
          pi_k <- sapply(1:K, function(i) pi_k[,1])
          log_pi_k <- sapply(1:K, function(i) log_pi_k[,1])
        }
      }
      if (counter <= 0) {
        auto_done <- TRUE
        break
      }
      
      if (which.max(elbo_new) == 1) {
        J <- J + increment 
        if (J > 1) {
          pi_j <- cbind(matrix(pi_j[,1], nrow = T, ncol = increment), pi_j)
          log_pi_j <- cbind(matrix(log_pi_j[,1], nrow = T, ncol = increment), log_pi_j)
        }
      } else if (which.max(elbo_new) == 2) {
        L <- L + increment 
        if (L > 1) {
          pi_l <- cbind(matrix(pi_l[,1], nrow = T, ncol = increment), pi_l)
          log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
        }
      } else {
        K <- K + increment
        if (K > 1) {
          pi_k <- cbind(matrix(pi_k[,1], nrow = T, ncol = increment), pi_k)
          log_pi_k <- cbind(matrix(log_pi_k[,1], nrow = T, ncol = increment), log_pi_k)
        }
      }
      if (auto_done) break
    }
  }
  
  #### reversing model ####
  if (reverse){
    if(verbose) print("reversing model")
    
    L <- fit$L; K <- fit$K; J <- fit$J;
    
    # reversed residuals, variance, and intercepts
    r_tilde <- fit$residual[T:1] 
    lambda_bar <- fit$lambda[T:1]
    delta <- fit$delta[T:1]
    mu_0 <- fit$mu[T]
    lambda_0 <- lambda_bar[1]
    
    # don't reverse weighted priors
    if (!pi_l_weighted) {
      pi_l <- pi_l[T:1,1:max(1,L),drop = FALSE]
      log_pi_l <- log_pi_l[T:1,1:max(1,L),drop = FALSE]
    } else {
      pi_l <- sapply(1:max(1,L), function(i) pi_l[,1])
      log_pi_l <- sapply(1:max(1,L), function(i) log_pi_l[,1])
    }
    if (!pi_k_weighted) {
      pi_k <- pi_k[T:1,1:max(1,K),drop = FALSE]
      log_pi_k <- log_pi_k[T:1,1:max(1,K),drop = FALSE]
    } else {
      pi_k <- sapply(1:max(1,K), function(i) pi_k[,1])
      log_pi_k <- sapply(1:max(1,K), function(i) log_pi_k[,1])
    }
    if (!pi_j_weighted) {
      pi_j <- pi_j[T:1,1:max(1,J),drop = FALSE]
      log_pi_j <- log_pi_j[T:1,1:max(1,J),drop = FALSE]
    } else {
      pi_j <- sapply(1:max(1,J), function(i) pi_j[,1])
      log_pi_j <- sapply(1:max(1,J), function(i) log_pi_j[,1])
    }
    
    # reversing mean components
    L_seq = seq_len(L)
    pi_bar_l <- matrix(nrow = T, ncol = L)
    log_pi_bar_l <- matrix(nrow = T, ncol = L)
    b_bar_l <- matrix(nrow = T, ncol = L)
    omega_bar_l <- matrix(nrow = T, ncol = L)
    
    for (l in L_seq) {
      tau_l <- which.max(fit$mean_model$pi_bar[,l])
      mu_bar_l <- fit$mean_model$mu_bar[,l]
      r_tilde_l <- fit$residual + mu_bar_l - mean(mu_bar_l[tau_l:T])
      r_tilde_l <- r_tilde_l[T:1]
      
      fit_scp <- mean_scp(r_tilde_l, lambda_bar, omega_l, log_pi_l[,l])
      
      pi_bar_l[,l] <- fit_scp$pi_bar
      log_pi_bar_l[,l] <- fit_scp$log_pi_bar
      b_bar_l[,l] <- fit_scp$b_bar
      omega_bar_l[,l] <- fit_scp$omega_bar
    }
    
    # reversing variance components
    K_seq = seq_len(K)
    pi_bar_k <- matrix(nrow = T, ncol = K)
    log_pi_bar_k <- matrix(nrow = T, ncol = K)
    v_bar_k <- matrix(nrow = T, ncol = K)
    
    for (k in K_seq) {
      tau_k <- which.max(fit$var_model$pi_bar[,k])
      lambda_bar_k <- fit$lambda / (fit$var_model$lambda_bar[,k] / mean(fit$var_model$lambda_bar[tau_k:T, k]))
      lambda_bar_k <- lambda_bar_k[T:1]
      
      fit_scp <- var_scp(r_tilde, lambda_bar_k, u_bar_k, lgamma_u_bar_k, 
                         rep(v_k, T), log_pi_k[,k])
      
      pi_bar_k[,k] <- fit_scp$pi_bar
      log_pi_bar_k[,k] <- fit_scp$log_pi_bar
      v_bar_k[,k] <- fit_scp$v_bar
    }
    
    # reversing mean-variance components 
    J_seq = seq_len(J)
    pi_bar_j <- matrix(nrow = T, ncol = J)
    log_pi_bar_j <- matrix(nrow = T, ncol = J)
    b_bar_j <- matrix(nrow = T, ncol = J)
    omega_bar_j <- matrix(nrow = T, ncol = J)
    v_bar_j <- matrix(nrow = T, ncol = J)
    
    for (j in J_seq) {
      tau_j <- which.max(fit$meanvar_model$pi_bar[,j])
      mu_lambda_bar_j <- fit$meanvar_model$mu_lambda_bar[,j]
      mu2_lambda_bar_j <- fit$meanvar_model$mu2_lambda_bar[,j]
      lambda_bar_j <- fit$meanvar_model$lambda_bar[,j]
      r_tilde_j <- fit$residual + mu_lambda_bar_j / lambda_bar_j - mean(mu_lambda_bar_j[tau_j:T] / lambda_bar_j[tau_j:T])
      r_tilde_j <- r_tilde_j[T:1]
      lambda_bar_j <- fit$lambda / (fit$meanvar_model$lambda_bar[,j] / mean(fit$meanvar_model$lambda_bar[tau_j:T,j]))
      lambda_bar_j <- lambda_bar_j[T:1]
      
      fit_scp <- meanvar_scp(r_tilde_j, lambda_bar_j, omega_j, u_bar_j, 
                             lgamma_u_bar_j, rep(v_j, T), log_pi_j[,j])
      
      pi_bar_j[,j] <- fit_scp$pi_bar
      log_pi_bar_j[,j] <- fit_scp$log_pi_bar
      b_bar_j[,j] <- fit_scp$b_bar
      omega_bar_j[,j] <- fit_scp$omega_bar
      v_bar_j[,j] <- fit_scp$v_bar
    }
    
    #### fit model ####
    fit <- mich_cpp(y[T:1], J, L, K, mu_0, lambda_0, 
                    fit_intercept, fit_scale, refit = TRUE,
                    max_iter, (verbose & !any(J_auto, K_auto, L_auto)), tol, 
                    omega_j, u_j, v_j, log_pi_j, 
                    pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j,
                    lgamma_u_bar_j, digamma_u_bar_j,
                    omega_l, log_pi_l, 
                    pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                    u_k, v_k, log_pi_k, 
                    pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                    lgamma_u_bar_k, digamma_u_bar_k)
  }
  
  #### return model ####
  return(fit)
}