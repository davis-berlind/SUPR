mich_vector_fast <- function(y, fit_intercept, fit_scale,
                             #merge = FALSE, merge_tol = 1 / T,
                             J, L, K, J_auto, L_auto, K_auto,
                             tol, B_l, B_r, verbose, max_iter, n_restart, n_search,
                             tau_j, u_j, v_j, pi_j,
                             tau_l, pi_l,
                             u_k, v_k, pi_k) {
  
  # calculate dimensions of y
  T <- length(y) - B_r - B_l
  
  # extract y in (1-B_l):1
  y_0 <- numeric(0)
  if (B_l > 0) {
    y_0 <- y[1:B_l]
    y <- y[-c(1:B_l)]
  }
  
  #### initialize mu_0 ####
  if (fit_intercept) {
    if (B_l >= ceiling(log(T))) mu_0 <- mean(y_0)
    else mu_0 <- mean(c(y_0, y[1:(ceiling(log(T)) - B_l)]))
  } else mu_0 <- 0.0
  
  #### initialize lambda_0 ####
  if (fit_scale) {
    if (B_l >= ceiling(log(T))) lambda_0 <- 1 / var(y_0)
    else lambda_0 <- 1 / var(c(y_0, y[1:(ceiling(log(T)) - B_l)]))
  } else lambda_0 <- 1.0
  
  #### initializing posterior parameters ####
  # J components
  pi_bar_j <- matrix(1 / T, nrow = T, ncol = J)
  b_bar_j <- matrix(0.0, nrow = T, ncol = J)
  tau_bar_j <- matrix(1.0, nrow = T, ncol = J)
  u_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  v_bar_j <- matrix(v_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  
  # L components
  pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
  b_bar_l <- matrix(0.0, nrow = T, ncol = L)
  tau_bar_l <- matrix(1.0, nrow = T, ncol = L)
  
  # K components
  pi_bar_k <- matrix(1 / T, nrow = T, ncol = K)
  u_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  v_bar_k <-  matrix(v_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  
  #### fit model ####
  fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                  B_l, B_r, fit_intercept, fit_scale, 
                  max_iter, (verbose & !any(J_auto, K_auto, L_auto)), tol, 
                  tau_j, u_j, v_j, log(pi_j), 
                  pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j, 
                  tau_l, log(pi_l), 
                  pi_bar_l, b_bar_l, tau_bar_l, 
                  u_k, v_k, log(pi_k), 
                  pi_bar_k, u_bar_k, v_bar_k)
  
  #### merging redundant components ####
  # while(TRUE & J > 0) {
  #   merge_probs <- t(pi_bar_j) %*% pi_bar_j
  #   merge_probs <- merge_probs[lower.tri(merge_probs)]
  # }
  
  #### auto procedure with single component ####
  if (sum(J_auto, L_auto, K_auto) == 1) {
    
    done <- FALSE
    counter <- n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo
    
    if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
    
    while (TRUE) {
      
      if (J_auto) {
        # increment dimension of parameters
        J <- J + 1 
        if (J > 1) pi_j <- cbind(pi_j[,1], pi_j)
        pi_bar_j <- cbind(1 / T, pi_bar_j)
        b_bar_j <- cbind(0.0, b_bar_j)
        tau_bar_j <- cbind(1.0, tau_bar_j)
        u_bar_j <- cbind(u_j + (T + B_r - 1:T + 1) / 2, u_bar_j) 
        v_bar_j <- cbind(v_j + (T + B_r - 1:T + 1) / 2, v_bar_j)
      }
      
      if (L_auto) {
        # increment dimension of parameters
        L <- L + 1
        if (L > 1) pi_l <- cbind(pi_l[,1], pi_l)
        pi_bar_l <- cbind(1 / T, pi_bar_l)
        b_bar_l <- cbind(0.0, b_bar_l)
        tau_bar_l <- cbind(1.0, tau_bar_l)
      }
      
      if (K_auto) {
        # increment dimension of parameters
        K <- K + 1
        if (K > 1) pi_k <- cbind(pi_k[,1], pi_k)
        pi_bar_k <- cbind(1 / T, pi_bar_k)
        u_bar_k <- cbind(u_k + (T + B_r - 1:T + 1) / 2, u_bar_k) 
        v_bar_k <- cbind(v_k + (T + B_r - 1:T + 1) / 2, v_bar_k)
      }
      
      # fit new model 
      fit_new <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                          B_l, B_r, fit_intercept, fit_scale, 
                          max_iter, verbose = FALSE, tol, 
                          tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, 
                          tau_bar_j, u_bar_j, v_bar_j,
                          tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                          u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k)
      
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
      
      if (elbo_new > elbo) {
        elbo <- elbo_new
        fit <- fit_new
      } else {  
        # if elbo decreases begin grid search
        if (J_auto) {
          J <- max(J - n_restart, 2)
          pi_j <- sapply(1:J, pi_j[,1])
        } else if (L_auto) {
          L <- max(L - n_restart, 2)
          pi_l <- sapply(1:L, pi_l[,1])
        } else {
          K <- max(K - n_restart, 2)
          pi_k <- sapply(1:K, pi_k[,1])
        }
        
        fit_new <- mich_vector(c(y_0, y), fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               tol, B_l, B_r, verbose = FALSE, max_iter,
                               tau_j, u_j, v_j, pi_j,
                               tau_l, pi_l,
                               u_k, v_k, pi_k)
        
        elbo_new <- fit_new$elbo[length(fit_new$elbo)]
        if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
        
        while (TRUE) {
          if (elbo_new > elbo) {
            fit <- fit_new 
            elbo <- elbo_new
            counter <- n_search
          } else {
            counter <- counter - 1
          }
          
          if (counter == 0) {
            done <- TRUE
            break
          }
          
          if (J_auto) {
            J <- J + 1
            pi_j <- cbind(pi_j[,1], pi_j)
          } else if (L_auto) {
            L <- L + 1 
            pi_l <- cbind(pi_l[,1], pi_l)
          } else {
            K <- K + 1
            pi_k <- cbind(pi_k[,1], pi_k)
          }
          
          fit_new <- mich_vector(c(y_0, y), fit_intercept, fit_scale, J, L, K,
                                 J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                                 tol, B_l, B_r, verbose = FALSE, max_iter,
                                 tau_j, u_j, v_j, pi_j,
                                 tau_l, pi_l,
                                 u_k, v_k, pi_k)
          
          elbo_new <- fit_new$elbo[length(fit_new$elbo)]
          if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
        }
      }
      if (done) break
    }
  }
  
  #### auto procedure with multiple components ####
  
  if (sum(J_auto, L_auto, K_auto) > 1) {
    
    done <- FALSE
    counter <- (J_auto + L_auto + K_auto) * n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo))
    
    while (TRUE) {
      elbo_new <- rep(-Inf, 3)
      if (J_auto) {
        J <- J + 1
        if (J > 1) pi_j <- cbind(pi_j[,1], pi_j)
        
        # fit new model 
        fit_new <- mich_vector(c(y_0, y), fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               tol, B_l, B_r, verbose = FALSE, max_iter,
                               tau_j, u_j, v_j, pi_j,
                               tau_l, pi_l,
                               u_k, v_k, pi_k)
        
        elbo_new[1] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new[1]))
        if (elbo_new[1] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[1]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        J <- J - 1
        if (J > 0) pi_j <- pi_j[,-1, drop = FALSE]
      }
      
      if (L_auto) {
        L <- L + 1
        if (L > 1) pi_l <- cbind(pi_l[,1], pi_l)
        
        # fit new model 
        fit_new <- mich_vector(c(y_0, y), fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               tol, B_l, B_r, verbose = FALSE, max_iter,
                               tau_j, u_j, v_j, pi_j,
                               tau_l, pi_l,
                               u_k, v_k, pi_k)
        
        elbo_new[2] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new[2]))
        if (elbo_new[2] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[2]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        L <- L - 1
        if (L > 0) pi_l <- pi_l[,-1, drop = FALSE]
      }
      
      if (K_auto) {
        K <- K + 1
        if (K > 1) pi_k <- cbind(pi_k[,1], pi_k)
        
        # fit new model 
        fit_new <- mich_vector(c(y_0, y), fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               tol, B_l, B_r, verbose = FALSE, max_iter,
                               tau_j, u_j, v_j, pi_j,
                               tau_l, pi_l,
                               u_k, v_k, pi_k)
        
        elbo_new[3] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new[3]))
        if (elbo_new[3] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[3]
          counter <- (J_auto + L_auto + K_auto) * n_search
        } else {
          counter <- counter - 1
        }
        K <- K - 1
        if (K > 0) pi_k <- pi_k[,-1, drop = FALSE]
      }
      if (counter <= 0) {
        done <- TRUE
        break
      }
      
      if (which.max(elbo_new) == 1) {
        J <- J + 1 
        if (J > 1) pi_j <- cbind(pi_j[,1], pi_j)
      } else if (which.max(elbo_new) == 2) {
        L <- L + 1 
        if (L > 1) pi_l <- cbind(pi_l[,1], pi_l)
      } else {
        K <- K + 1
        if (K > 1) pi_k <- cbind(pi_k[,1], pi_k)
      }
      if (done) break
    }
  }
  
  #### return model ####
  # reassemble y
  fit$y <- c(y_0, y)
  
  # add zero probabilities in buffer 
  if (fit$J > 0) fit$J_model$pi_bar <- rbind(matrix(0, nrow = B_l, ncol = fit$J), 
                                             fit$J_model$pi_bar,
                                             matrix(0, nrow = B_r, ncol = fit$J))
  if (fit$L > 0) fit$L_model$pi_bar <- rbind(matrix(0, nrow = B_l, ncol = fit$L), 
                                             fit$L_model$pi_bar,
                                             matrix(0, nrow = B_r, ncol = fit$L))
  if (fit$K > 0) fit$K_model$pi_bar <- rbind(matrix(0, nrow = B_l, ncol = fit$K), 
                                             fit$K_model$pi_bar,
                                             matrix(0, nrow = B_r, ncol = fit$K))
  return(fit)
}