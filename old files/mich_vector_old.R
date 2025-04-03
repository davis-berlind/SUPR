mich_vector <- function(y, fit_intercept, fit_scale,
                        J, L, K, J_auto, L_auto, K_auto,
                        tol, B_l, B_r, verbose, max_iter, 
                        detect, merge_level, merge_prob, n_restart, n_search,
                        omega_j, u_j, v_j, pi_j,
                        omega_l, pi_l,
                        u_k, v_k, pi_k)

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
  # mean components
  pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
  b_bar_l <- matrix(0.0, nrow = T, ncol = L)
  omega_bar_l <- matrix(1.0, nrow = T, ncol = L)
  
  # variance components
  pi_bar_k <- matrix(1 / T, nrow = T, ncol = K)
  u_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  v_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  
  # mean-variance components
  pi_bar_j <- matrix(1 / T, nrow = T, ncol = J)
  b_bar_j <- matrix(0.0, nrow = T, ncol = J)
  omega_bar_j <- matrix(1.0, nrow = T, ncol = J)
  u_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  v_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  
  #### fit model ####
  fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                  B_l, B_r, fit_intercept, fit_scale, refit = FALSE,
                  max_iter, (verbose & !any(J_auto, K_auto, L_auto)), tol, 
                  omega_j, u_j, v_j, log(pi_j), 
                  pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j, 
                  omega_l, log(pi_l), 
                  pi_bar_l, b_bar_l, omega_bar_l, 
                  u_k, v_k, log(pi_k), 
                  pi_bar_k, u_bar_k, v_bar_k)
  
  #### merging redundant components ####
  cntr <- 0
  while (TRUE) {
    no_merges <- TRUE
    
    if (L > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      
      # compute adjacency matrix from pairwise merge probabilities and keep_mat
      merge_mat <- (t(pi_bar_l) %*% pi_bar_l > merge_prob * 1.1^cntr) & keep_mat
      
      if (any(merge_mat[lower.tri(merge_mat)])) {
        no_merges <- FALSE
        
        # force transitivity
        merge_mat <- force_adj(merge_mat)
        
        # combine and reset components
        for (i in 1:(L-1)) {
          for (j in (i+1):L) {
            if (merge_mat[i, j]) {
              if (max(pi_bar_l[,i]) < max(pi_bar_l[,j])) pi_bar_l[,i] <- pi_bar_l[,j]
              pi_bar_l[,j] <- 1 / T
              
              b_bar_l[,i] <- b_bar_l[,i] + b_bar_l[,j]
              b_bar_l[,j] <- 0.0
              
              omega_bar_l[,i] <- (omega_bar_l[,i] + omega_bar_l[,j]) / 2
              omega_bar_l[,j] <- 1.0         
            }
          }
        }
        order_l <- order(apply(pi_bar_l, 2, max))
        pi_bar_l <- pi_bar_l[,order_l]
        b_bar_l <- b_bar_l[,order_l]
        omega_bar_l <- omega_bar_l[,order_l]
      }
    }
    
    if (K > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = K, nrow = K)
      keep_mat[keep, keep] <- TRUE

      # compute adjacency matrix from pairwise merge probabilities 
      merge_mat <- (t(pi_bar_k) %*% pi_bar_k > merge_prob * 1.1^cntr) & keep_mat
      
      if (any(merge_mat[lower.tri(merge_mat)])) {
        no_merges <- FALSE
        
        # force transitivity
        merge_mat <- force_adj(merge_mat)
        
        # combine and reset components
        for (i in 1:(K-1)) {
          for (j in (i+1):K) {
            if (merge_mat[i, j]) {
              if (max(pi_bar_k[,i]) < max(pi_bar_k[,j])) pi_bar_k[,i] <- pi_bar_k[,j]
              pi_bar_k[,j] <- 1 / T
              
              v_bar_k[,i] <- v_bar_k[,i] * v_bar_k[,j] / u_bar_k[,j]
              v_bar_k[,j] <- u_bar_k[,j]
            }
          }
        }
        order_k <- order(apply(pi_bar_k, 2, max))
        pi_bar_k <- pi_bar_k[,order_k]
        v_bar_k <- v_bar_k[,order_k]
      }
    }
    
    if (J > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = J, nrow = J)
      keep_mat[keep, keep] <- TRUE

      # compute adjacency matrix from pairwise merge probabilities 
      merge_mat <- (t(pi_bar_j) %*% pi_bar_j > merge_prob * 1.1^cntr) & keep_mat
      
      if (any(merge_mat[lower.tri(merge_mat)])) {
        no_merges <- FALSE
        
        # force transitivity
        merge_mat <- force_adj(merge_mat)
        
        # combine and reset components
        for (i in 1:(J-1)) {
          for (j in (i+1):J) {
            if (merge_mat[i, j]) {

              if (max(pi_bar_j[,i]) < max(pi_bar_j[,j])) pi_bar_j[,i] <- pi_bar_j[,j]
              pi_bar_j[,j] <- 1 / T

              b_bar_j[,i] <- b_bar_j[,i] + b_bar_j[,j]
              b_bar_j[,j] <- 0.0

              omega_bar_j[,i] <- (omega_bar_j[,i] + omega_bar_j[,j]) / 2
              omega_bar_j[,j] <- 1.0    

              v_bar_j[,i] <- v_bar_j[,i] * v_bar_j[,j] / u_bar_j[,j]
              v_bar_j[,j] <- u_bar_j[,j]
            }
          }
        }
        order_j <- order(apply(pi_bar_j, 2, max))
        pi_bar_j <- pi_bar_j[,order_j]
        b_bar_j <- b_bar_j[,order_j]
        omega_bar_j <- omega_bar_j[,order_j]
        v_bar_j <- v_bar_j[,order_j]
      }
    }
    
    # check no more merges otherwise refit
    if (no_merges) break
    if (verbose & !any(J_auto, K_auto, L_auto)) print("Merging duplicate components")
    
    fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                    B_l, B_r, fit_intercept, fit_scale, refit = TRUE,
                    max_iter, (verbose & !any(J_auto, K_auto, L_auto)), tol, 
                    omega_j, u_j, v_j, log(pi_j), 
                    pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j, 
                    omega_l, log(pi_l), 
                    pi_bar_l, b_bar_l, omega_bar_l, 
                    u_k, v_k, log(pi_k), 
                    pi_bar_k, u_bar_k, v_bar_k)
    cntr <- cntr + 1
  }

  #### auto procedure with single component ####
  if (sum(J_auto, L_auto, K_auto) == 1) {
    
    counter <- n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo
    
    if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
    
    while (TRUE) {
      auto_done <- FALSE
      
      if (J_auto) {
        # increment dimension of parameters
        J <- J + 1 
        if (J > 1) pi_j <- cbind(pi_j[,1], pi_j)
        pi_bar_j <- cbind(1 / T, pi_bar_j)
        b_bar_j <- cbind(0.0, b_bar_j)
        omega_bar_j <- cbind(1.0, omega_bar_j)
        u_bar_j <- cbind(u_j + (T + B_r - 1:T + 1) / 2, u_bar_j) 
        v_bar_j <- cbind(v_j + (T + B_r - 1:T + 1) / 2, v_bar_j)
      }
      
      if (L_auto) {
        # increment dimension of parameters
        L <- L + 1
        if (L > 1) pi_l <- cbind(pi_l[,1], pi_l)
        pi_bar_l <- cbind(1 / T, pi_bar_l)
        b_bar_l <- cbind(0.0, b_bar_l)
        omega_bar_l <- cbind(1.0, omega_bar_l)
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
                          B_l, B_r, fit_intercept, fit_scale, refit = TRUE, 
                          max_iter, verbose = FALSE, tol, 
                          omega_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, 
                          omega_bar_j, u_bar_j, v_bar_j,
                          omega_l, log(pi_l), pi_bar_l, b_bar_l, omega_bar_l, 
                          u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k)
      
      
      # check for merges
      cntr <- 0
      while (TRUE) {
        no_merges <- TRUE
        
        if (L > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = L, nrow = L)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities and keep_mat
          merge_mat <- (t(pi_bar_l) %*% pi_bar_l > merge_prob * 1.1^cntr) & keep_mat
          
          if (any(merge_mat[lower.tri(merge_mat)])) {
            no_merges <- FALSE
            
            # force transitivity
            merge_mat <- force_adj(merge_mat)
            
            # combine and reset components
            for (i in 1:(L-1)) {
              for (j in (i+1):L) {
                if (merge_mat[i, j]) {
                  if (max(pi_bar_l[,i]) < max(pi_bar_l[,j])) pi_bar_l[,i] <- pi_bar_l[,j]
                  pi_bar_l[,j] <- 1 / T
                  
                  b_bar_l[,i] <- b_bar_l[,i] + b_bar_l[,j]
                  b_bar_l[,j] <- 0.0
                  
                  omega_bar_l[,i] <- (omega_bar_l[,i] + omega_bar_l[,j]) / 2
                  omega_bar_l[,j] <- 1.0         
                }
              }
            }
            order_l <- order(apply(pi_bar_l, 2, max))
            pi_bar_l <- pi_bar_l[,order_l]
            b_bar_l <- b_bar_l[,order_l]
            omega_bar_l <- omega_bar_l[,order_l]
          }
        }
        
        if (K > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = K, nrow = K)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities 
          merge_mat <- (t(pi_bar_k) %*% pi_bar_k > merge_prob * 1.1^cntr) & keep_mat
          
          if (any(merge_mat[lower.tri(merge_mat)])) {
            no_merges <- FALSE
            
            # force transitivity
            merge_mat <- force_adj(merge_mat)
            
            # combine and reset components
            for (i in 1:(K-1)) {
              for (j in (i+1):K) {
                if (merge_mat[i, j]) {
                  if (max(pi_bar_k[,i]) < max(pi_bar_k[,j])) pi_bar_k[,i] <- pi_bar_k[,j]
                  pi_bar_k[,j] <- 1 / T
                  
                  v_bar_k[,i] <- v_bar_k[,i] * v_bar_k[,j] / u_bar_k[,j]
                  v_bar_k[,j] <- u_bar_k[,j]
                }
              }
            }
            order_k <- order(apply(pi_bar_k, 2, max))
            pi_bar_k <- pi_bar_k[,order_k]
            v_bar_k <- v_bar_k[,order_k]
          }
        }
        
        if (J > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = J, nrow = J)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities 
          merge_mat <- (t(pi_bar_j) %*% pi_bar_j > merge_prob * 1.1^cntr) & keep_mat
          
          if (any(merge_mat[lower.tri(merge_mat)])) {
            no_merges <- FALSE
            
            # force transitivity
            merge_mat <- force_adj(merge_mat)
            
            # combine and reset components
            for (i in 1:(J-1)) {
              for (j in (i+1):J) {
                if (merge_mat[i, j]) {
                  
                  if (max(pi_bar_j[,i]) < max(pi_bar_j[,j])) pi_bar_j[,i] <- pi_bar_j[,j]
                  pi_bar_j[,j] <- 1 / T
                  
                  b_bar_j[,i] <- b_bar_j[,i] + b_bar_j[,j]
                  b_bar_j[,j] <- 0.0
                  
                  omega_bar_j[,i] <- (omega_bar_j[,i] + omega_bar_j[,j]) / 2
                  omega_bar_j[,j] <- 1.0    
                  
                  v_bar_j[,i] <- v_bar_j[,i] * v_bar_j[,j] / u_bar_j[,j]
                  v_bar_j[,j] <- u_bar_j[,j]
                }
              }
            }
            order_j <- order(apply(pi_bar_j, 2, max))
            pi_bar_j <- pi_bar_j[,order_j]
            b_bar_j <- b_bar_j[,order_j]
            omega_bar_j <- omega_bar_j[,order_j]
            v_bar_j <- v_bar_j[,order_j]
          }
        }
        
        # check no more merges otherwise refit
        if (no_merges) break

        fit_new <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                            B_l, B_r, fit_intercept, fit_scale, refit = TRUE,
                            max_iter, verbose = FALSE, tol, 
                            omega_j, u_j, v_j, log(pi_j), 
                            pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j, 
                            omega_l, log(pi_l), 
                            pi_bar_l, b_bar_l, omega_bar_l, 
                            u_k, v_k, log(pi_k), 
                            pi_bar_k, u_bar_k, v_bar_k)
        cntr <- cntr + 1
      }
      
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
      
      if (elbo_new > elbo) {
        elbo <- elbo_new
        fit <- fit_new
      } else {  
        # if elbo decreases begin grid search
        if (J_auto) {
          J <- max(J - n_restart, 2)
          pi_j <- sapply(1:J, function(i) pi_j[,1])
        } else if (L_auto) {
          L <- max(L - n_restart, 2)
          pi_l <- sapply(1:L, function(i) pi_l[,1])
        } else {
          K <- max(K - n_restart, 2)
          pi_k <- sapply(1:K, function(i) pi_k[,1])
        }
        
        fit_new <- mich_vector(c(y_0, y), fit_intercept, fit_scale, J, L, K,
                               J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                               tol, B_l, B_r, verbose = FALSE, max_iter, 
                               merge_prob, n_restart, n_search,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l,
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
            auto_done <- TRUE
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
                                 merge_prob, n_restart, n_search,
                                 omega_j, u_j, v_j, pi_j,
                                 omega_l, pi_l,
                                 u_k, v_k, pi_k)
          
          elbo_new <- fit_new$elbo[length(fit_new$elbo)]
          if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
        }
      }
      if (auto_done) break
    }
  }
  
  #### auto procedure with multiple components ####
  
  if (sum(J_auto, L_auto, K_auto) > 1) {
    
    auto_done <- FALSE
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
                               merge_prob, n_restart, n_search,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l,
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
                               merge_prob, n_restart, n_search,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l,
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
                               merge_prob, n_restart, n_search,
                               omega_j, u_j, v_j, pi_j,
                               omega_l, pi_l,
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
        auto_done <- TRUE
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
      if (auto_done) break
    }
  }
  
  #### return model ####
  # reassemble y
  fit$y <- c(y_0, y)
  
  # add zero probabilities in buffer 
  if (B_l > 0) {
    if (fit$J > 0) fit$meanvar_model$pi_bar <- rbind(matrix(0, nrow = B_l, ncol = fit$J), 
                                                     fit$J_model$pi_bar,
                                                     matrix(0, nrow = B_r, ncol = fit$J))
    if (fit$L > 0) fit$mean_model$pi_bar <- rbind(matrix(0, nrow = B_l, ncol = fit$L), 
                                                  fit$L_model$pi_bar,
                                                  matrix(0, nrow = B_r, ncol = fit$L))
    if (fit$K > 0) fit$var_model$pi_bar <- rbind(matrix(0, nrow = B_l, ncol = fit$K), 
                                                 fit$K_model$pi_bar,
                                                 matrix(0, nrow = B_r, ncol = fit$K))
  }

  return(fit)
}