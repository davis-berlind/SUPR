mich_vector_meanvar <- function(y, fit_intercept, fit_scale,
                        J, J_auto, pi_j_weighted,
                        tol, verbose, max_iter, reverse,
                        detect, merge_level, merge_prob,
                        restart, n_restart, n_search,
                        omega_j, u_j, v_j, pi_j) {

  #### set up ####
  # seeding sequence to increase merge probability
  merge_seq <- 2
  
  # calculate dimensions of y
  T <- length(y)
  
  #### initialize mu_0 ####
  if (fit_intercept) {
    mu_0 <- mean(y[1:ceiling(log(T))])
  } else mu_0 <- 0.0
  
  #### initialize lambda_0 ####
  if (fit_scale) {
    lambda_0 <- 1 / var(y[1:ceiling(log(T))])
  } else lambda_0 <- 1.0
  
  # initializing meanvar with mean model
  L <- J
  J <- 0
  
  #### initializing posterior parameters ####
  # mean components
  pi_l <- pi_j
  log_pi_l <- log(pi_l)
  pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
  log_pi_bar_l <- matrix(0, nrow = T, ncol = L)
  b_bar_l <- matrix(0.0, nrow = T, ncol = L)
  omega_bar_l <- matrix(1.0, nrow = T, ncol = L)

  # variance components
  K <- 0
  pi_k <- matrix(nrow = T, ncol = K)
  log_pi_k <- matrix(nrow = T, ncol = K)
  pi_bar_k <- matrix(nrow = T, ncol = K)
  log_pi_bar_k <- matrix(nrow = T, ncol = K)
  u_bar_k <- rep(0, T)
  lgamma_u_bar_k <- rep(0, T)
  digamma_u_bar_k <- rep(0, T)
  v_bar_k <- matrix(nrow = T, ncol = K)
  
  # mean-variance components
  log_pi_j <- log(pi_j)
  pi_bar_j <- matrix(nrow = T, ncol = J)
  log_pi_bar_j <- matrix(nrow = T, ncol = J)
  b_bar_j <- matrix(nrow = T, ncol = J)
  omega_bar_j <- matrix(nrow = T, ncol = J)
  u_bar_j <- u_j + (T - 1:T + 1) / 2
  lgamma_u_bar_j <- lgamma(u_bar_j)
  digamma_u_bar_j <- digamma(u_bar_j)
  v_bar_j <- matrix(nrow = T, ncol = J)
  
  #### fit model ####
  # initialize mean model
  cntr <- 0
  while (TRUE) {
    no_merges <- TRUE
    fit <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                    fit_intercept, fit_scale, refit = cntr > 0,
                    max_iter, verbose = verbose & !(J_auto), tol, 
                    omega_j, u_j, v_j, log_pi_j, 
                    pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j,
                    lgamma_u_bar_j, digamma_u_bar_j,
                    omega_l, log_pi_l, 
                    pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                    u_k, v_k, log_pi_k, 
                    pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                    lgamma_u_bar_k, digamma_u_bar_k)
    
    if (L > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      
      # compute adjacency matrix from pairwise merge probabilities and keep_mat
      merge_mat <- (t(pi_bar_l) %*% pi_bar_l > merge_prob * merge_seq^cntr) & keep_mat
      
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
              log_pi_bar_l[,j] <- 0
              
              b_bar_l[,i] <- b_bar_l[,i] + b_bar_l[,j]
              b_bar_l[,j] <- 0.0
              
              omega_bar_l[,i] <- (omega_bar_l[,i] + omega_bar_l[,j]) / 2
              omega_bar_l[,j] <- 1.0         
            }
          }
        }
        order_l <- order(apply(pi_bar_l, 2, max))
        pi_bar_l <- pi_bar_l[,order_l]
        log_pi_bar_l <- log_pi_bar_l[,order_l]
        b_bar_l <- b_bar_l[,order_l]
        omega_bar_l <- omega_bar_l[,order_l]
      }
    }
    if (no_merges) break
    cntr <- cntr + 1
    if (verbose & !(J_auto)) print(paste("Merging duplicate components iteration:", cntr))
  }
  
  # drop mean components
  J <- L
  L <- 0
  
  if (J > 0) {
    pi_l <- pi_l[, 1, drop = FALSE]
    log_pi_l <- log_pi_l[, 1, drop = FALSE]
    pi_bar_l <- matrix(nrow = T, ncol = L)
    log_pi_bar_l <- matrix(nrow = T, ncol = L)
    b_bar_l <- matrix(nrow = T, ncol = L)
    omega_bar_l <- matrix(nrow = T, ncol = L)
    
    # add mean-variance components
    pi_bar_j <- fit$mean_model$pi_bar
    log_pi_bar_j <- log(pi_bar_j)
    b_bar_j <- fit$mean_model$b_bar
    omega_bar_j <- fit$mean_model$omega_bar
    v_bar_j <- sapply(1:max(1,J), function(i) u_bar_j)
    
    # fit meanvar model with mean model initialization. 
    if (verbose & !J_auto) print("Refiting with mean model initialization.")
    cntr <- 0
    while (TRUE) {
      no_merges <- TRUE
      fit <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                      fit_intercept, fit_scale, refit = TRUE,
                      max_iter, verbose = verbose & !(J_auto), tol, 
                      omega_j, u_j, v_j, log_pi_j, 
                      pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j,
                      lgamma_u_bar_j, digamma_u_bar_j,
                      omega_l, log_pi_l, 
                      pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                      u_k, v_k, log_pi_k, 
                      pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                      lgamma_u_bar_k, digamma_u_bar_k)
      
      
      if (J > 1) {
        # test which columns actually contain change-points
        cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
        keep <- sapply(cred_sets, length) <= detect
        keep_mat <- matrix(FALSE, ncol = J, nrow = J)
        keep_mat[keep, keep] <- TRUE
        
        # compute adjacency matrix from pairwise merge probabilities 
        merge_mat <- (t(pi_bar_j) %*% pi_bar_j > merge_prob * merge_seq^cntr) & keep_mat
        
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
                log_pi_bar_j[,j] <- 0
                
                b_bar_j[,i] <- b_bar_j[,i] + b_bar_j[,j]
                b_bar_j[,j] <- 0.0
                
                omega_bar_j[,i] <- (omega_bar_j[,i] + omega_bar_j[,j]) / 2
                omega_bar_j[,j] <- 1.0    
                
                v_bar_j[,i] <- v_bar_j[,i] * v_bar_j[,j] / u_bar_j
                v_bar_j[,j] <- u_bar_j
              }
            }
          }
          order_j <- order(apply(pi_bar_j, 2, max))
          pi_bar_j <- pi_bar_j[,order_j]
          log_pi_bar_j <- log_pi_bar_j[,order_j]
          b_bar_j <- b_bar_j[,order_j]
          omega_bar_j <- omega_bar_j[,order_j]
          v_bar_j <- v_bar_j[,order_j]
        }
      }
      if (no_merges) break
      cntr <- cntr + 1
      if (verbose & !(J_auto)) print(paste("Merging duplicate components iteration:", cntr))
    }
  }
  
  #### auto procedure with single component ####
  last_restart <- ifelse(restart, 2, Inf)
  
  if (J_auto) {
    
    refit <- (J > 0)
    counter <- n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo
    
    if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new))
    while (TRUE) {
      if (J_auto) {
        # increment dimension of parameters
        L <- L + 1
        if (L > 1) {
          pi_l <- cbind(pi_l[,1], pi_l)
          log_pi_l <- cbind(log_pi_l[,1], log_pi_l)
        }
        pi_bar_l <- cbind(1 / T, pi_bar_l)
        log_pi_bar_l <- cbind(0.0, log_pi_bar_l)
        b_bar_l <- cbind(0.0, b_bar_l)
        omega_bar_l <- cbind(1.0, omega_bar_l)
      }
      
      # fit new model 
      cntr <- 0
      while (TRUE) {
        no_merges <- TRUE
        fit_new <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                            fit_intercept, fit_scale, refit,
                            max_iter, verbose = verbose & !(J_auto), tol, 
                            omega_j, u_j, v_j, log_pi_j, 
                            pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j,
                            lgamma_u_bar_j, digamma_u_bar_j,
                            omega_l, log_pi_l, 
                            pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                            u_k, v_k, log_pi_k, 
                            pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                            lgamma_u_bar_k, digamma_u_bar_k)
        
        if (L > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = L, nrow = L)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities and keep_mat
          merge_mat <- (t(pi_bar_l) %*% pi_bar_l > merge_prob * merge_seq^cntr) & keep_mat
          
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
                  log_pi_bar_l[,j] <- 0
                  
                  b_bar_l[,i] <- b_bar_l[,i] + b_bar_l[,j]
                  b_bar_l[,j] <- 0.0
                  
                  omega_bar_l[,i] <- (omega_bar_l[,i] + omega_bar_l[,j]) / 2
                  omega_bar_l[,j] <- 1.0         
                }
              }
            }
            order_l <- order(apply(pi_bar_l, 2, max))
            pi_bar_l <- pi_bar_l[,order_l]
            log_pi_bar_l <- log_pi_bar_l[,order_l]
            b_bar_l <- b_bar_l[,order_l]
            omega_bar_l <- omega_bar_l[,order_l]
          }
        }
        if (no_merges) break
        cntr <- cntr + 1
      }
      
      # drop mean component 
      if (L > 1) {
        J <- L
        pi_j <- pi_l
        log_pi_j <- log_pi_bar_l
        b_bar_j <- b_bar_l
        pi_bar_j <- pi_bar_l
        log_pi_bar_j <- log_pi_bar_l
        omega_bar_j <- omega_bar_l
        v_bar_j <- sapply(1:J, function(i) u_bar_j)
      } else {
        J <- J + 1
        if (J > 1) {
          pi_j <- cbind(pi_j[,1], pi_j)
          log_pi_j <- cbind(log_pi_j[,1], log_pi_j)
        }        
        b_bar_j <- cbind(b_bar_l, b_bar_j)
        pi_bar_j <- cbind(pi_bar_l, pi_bar_j)
        log_pi_bar_j <- cbind(log_pi_bar_l, log_pi_bar_j)
        omega_bar_j <- cbind(omega_bar_l, omega_bar_j)
        v_bar_j <- cbind(u_bar_j, v_bar_j)
      }
      
      L <- 0
      pi_l <- pi_l[, 1, drop = FALSE]
      log_pi_l <- log_pi_l[, 1, drop = FALSE]
      pi_bar_l <- matrix(nrow = T, ncol = L)
      log_pi_bar_l <- matrix(nrow = T, ncol = L)
      b_bar_l <- matrix(nrow = T, ncol = L)
      omega_bar_l <- matrix(nrow = T, ncol = L)
      
      # fit meanvar model with mean model initialization. 
      cntr <- 0
      while (TRUE) {
        no_merges <- TRUE
        fit_new <- mich_cpp(y, J, L, K, mu_0, lambda_0,
                            fit_intercept, fit_scale, refit = TRUE,
                            max_iter, verbose = TRUE, tol, 
                            omega_j, u_j, v_j, log_pi_j, 
                            pi_bar_j, log_pi_bar_j, b_bar_j, omega_bar_j, u_bar_j, v_bar_j,
                            lgamma_u_bar_j, digamma_u_bar_j,
                            omega_l, log_pi_l, 
                            pi_bar_l, log_pi_bar_l, b_bar_l, omega_bar_l, 
                            u_k, v_k, log_pi_k, 
                            pi_bar_k, log_pi_bar_k, u_bar_k, v_bar_k,
                            lgamma_u_bar_k, digamma_u_bar_k)
        
        if (J > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = J, nrow = J)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities 
          merge_mat <- (t(pi_bar_j) %*% pi_bar_j > merge_prob * merge_seq^cntr) & keep_mat
          
          if (any(merge_mat[lower.tri(merge_mat)])) {
            no_merges <- FALSE
            stop("here")
            # force transitivity
            merge_mat <- force_adj(merge_mat)
            
            # combine and reset components
            for (i in 1:(J-1)) {
              for (j in (i+1):J) {
                if (merge_mat[i, j]) {
                  
                  if (max(pi_bar_j[,i]) < max(pi_bar_j[,j])) pi_bar_j[,i] <- pi_bar_j[,j]
                  pi_bar_j[,j] <- 1 / T
                  log_pi_bar_j[,j] <- 0
                  
                  b_bar_j[,i] <- b_bar_j[,i] + b_bar_j[,j]
                  b_bar_j[,j] <- 0.0
                  
                  omega_bar_j[,i] <- (omega_bar_j[,i] + omega_bar_j[,j]) / 2
                  omega_bar_j[,j] <- 1.0    
                  
                  v_bar_j[,i] <- v_bar_j[,i] * v_bar_j[,j] / u_bar_j
                  v_bar_j[,j] <- u_bar_j
                }
              }
            }
            order_j <- order(apply(pi_bar_j, 2, max))
            pi_bar_j <- pi_bar_j[,order_j]
            log_pi_bar_j <- log_pi_bar_j[,order_j]
            b_bar_j <- b_bar_j[,order_j]
            omega_bar_j <- omega_bar_j[,order_j]
            v_bar_j <- v_bar_j[,order_j]
          }
        }
        if (no_merges) break
        cntr <- cntr + 1
      }
      
      if (!refit) refit <- TRUE
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("(L = ", L, ", K = ", K, ", J = ", J, "): ELBO = ", elbo_new))
      
      if (elbo_new > elbo) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- n_search
        if (last_restart < Inf) restart <- TRUE 
      } else if (restart & J_auto & J - n_restart >= last_restart) { 
        refit <- FALSE
        restart <- FALSE
        last_restart <- J
        counter <- counter + n_restart
        
        L <- J - n_restart - 1
        J <- 0
        pi_j <- pi_j[, 1, drop = FALSE]
        log_pi_j <- log_pi_j[, 1, drop = FALSE]
        pi_bar_j <- matrix(nrow = T, ncol = J)
        log_pi_bar_j <- matrix(nrow = T, ncol = J)
        b_bar_j <- matrix(nrow = T, ncol = J)
        omega_bar_j <- matrix(nrow = T, ncol = J)
        v_bar_j <- matrix(nrow = T, ncol = J)
        
        # mean components
        pi_l <- sapply(1:L, function(i) pi_j)
        log_pi_l <- sapply(1:L, function(i) log_pi_j)
        pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
        log_pi_bar_l <- matrix(0, nrow = T, ncol = L)
        b_bar_l <- matrix(0.0, nrow = T, ncol = L)
        omega_bar_l <- matrix(1.0, nrow = T, ncol = L)
      } else {
        counter <- counter - 1
      }
      if (counter == 0) break
    }
  }
  
  #### reversing model ####
  if (reverse){
    print("reversing model")
    
    L <- fit$L; K <- fit$K; J <- fit$J;
    
    # reversed residuals, variance, and intercepts
    r_tilde <- fit$residual[T:1] 
    lambda_bar <- fit$lambda[T:1]
    delta <- fit$delta[T:1]
    mu_0 <- fit$mu[T]
    lambda_0 <- lambda_bar[1]
    
    # don't reverse weighted priors
    if (!pi_j_weighted) {
      pi_j <- pi_j[1:T,1:max(1,J),drop = FALSE]
      log_pi_j <- log_pi_j[1:T,1:max(1,J),drop = FALSE]
    }
    
    # reversing mean-variance components 
    J_seq = seq_len(J)
    pi_bar_j <- pi_bar_j[, J_seq, drop = FALSE]
    log_pi_bar_j <- log_pi_bar_j[, J_seq, drop = FALSE]
    b_bar_j <- b_bar_j[, J_seq, drop = FALSE]
    omega_bar_j <- omega_bar_j[, J_seq, drop = FALSE]
    v_bar_j <- v_bar_j[, J_seq, drop = FALSE]
    
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