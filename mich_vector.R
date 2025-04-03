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
  # seeding sequence to increase merge probability
  # merge_seq <- exp((2/5) * (log(T) - log(log(T))))
  
  # Flag for new model 
  refit <- FALSE
  
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
                       restart, n_restart, n_search, increment,
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
      
      # compute adjacency matrix from pairwise merge probabilities and keep_mat
      merge_mat <- (t(pi_bar_l) %*% pi_bar_l > merge_prob) & keep_mat
      
      if (any(merge_mat[lower.tri(merge_mat)])) {
        merged <- FALSE
        
        # force transitivity
        merge_mat <- force_adj(merge_mat)
        diag(merge_mat) <- TRUE
        
        # combine and reset components
        keep <- which(apply(merge_mat, 2, sum) == 1)
        visited <- keep
        for (i in 1:L) {
          if (i %in% visited) next
          merge_cand <- which(merge_mat[,i])
          merger <- merge_cand[which.max(apply(pi_bar_l[,merge_cand], 2, max))]
          mergees <- merge_cand[merge_cand != merger]
          keep <- c(merger, keep)
          visited <- c(merge_cand , visited)
          
          b_bar_l[,merger] <- rowSums(b_bar_l[,merge_cand])
          omega_bar_l[,merger] <- rowMeans(omega_bar_l[,merge_cand])
        }
        
        L <- length(keep)
        if (L_auto) pi_l <- pi_l[,keep, drop = FALSE]
        if (L_auto) log_pi_l <- log_pi_l[,keep, drop = FALSE]
        pi_bar_l <- pi_bar_l[,keep, drop = FALSE]
        log_pi_bar_l <- log_pi_bar_l[,keep, drop = FALSE]
        b_bar_l <- b_bar_l[,keep, drop = FALSE]
        omega_bar_l <- omega_bar_l[,keep, drop = FALSE]
      }
    }
    
    if (K > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = K, nrow = K)
      keep_mat[keep, keep] <- TRUE
      
      # compute adjacency matrix from pairwise merge probabilities 
      merge_mat <- (t(pi_bar_k) %*% pi_bar_k > merge_prob) & keep_mat
      
      if (any(merge_mat[lower.tri(merge_mat)])) {
        merged <- FALSE
        
        # force transitivity
        merge_mat <- force_adj(merge_mat)
        diag(merge_mat) <- TRUE
        
        # combine and reset components
        keep <- which(apply(merge_mat, 2, sum) == 1)
        visited <- keep
        for (i in 1:K) {
          if (i %in% visited) next
          merge_cand <- which(merge_mat[,i])
          merger <- merge_cand[which.max(apply(pi_bar_k[,merge_cand], 2, max))]
          mergees <- merge_cand[merge_cand != merger]
          keep <- c(merger, keep)
          visited <- c(merge_cand , visited)
          
          v_bar_k[,merger] <- apply(v_bar_k[,merge_cand], 1, prod) / u_bar_k
          v_bar_k[,mergees] <- u_bar_k
        }
        
        K <- length(keep)
        if (K_auto) pi_k <- pi_k[,keep, drop = FALSE]
        if (K_auto) log_pi_k <- log_pi_k[,keep, drop = FALSE]
        pi_bar_k <- pi_bar_k[,keep, drop = FALSE]
        log_pi_bar_k <- log_pi_bar_k[,keep, drop = FALSE]
        v_bar_k <- v_bar_k[,keep, drop = FALSE]
      }
    }
    
    if (J > 1) {
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = J, nrow = J)
      keep_mat[keep, keep] <- TRUE
      
      # compute adjacency matrix from pairwise merge probabilities 
      merge_mat <- (t(pi_bar_j) %*% pi_bar_j > merge_prob) & keep_mat
      
      if (any(merge_mat[lower.tri(merge_mat)])) {
        merged <- FALSE
        
        # force transitivity
        merge_mat <- force_adj(merge_mat)
        diag(merge_mat) <- TRUE
        
        # combine and reset components
        keep <- which(apply(merge_mat, 2, sum) == 1)
        visited <- keep            
        for (i in 1:J) {
          if (i %in% visited) next
          merge_cand <- which(merge_mat[,i])
          merger <- merge_cand[which.max(apply(pi_bar_j[,merge_cand], 2, max))]
          mergees <- merge_cand[merge_cand != merger]
          keep <- c(merger, keep)
          visited <- c(merge_cand , visited)
          
          b_bar_j[,merger] <- rowSums(b_bar_j[,merge_cand])
          omega_bar_j[,merger] <- rowMeans(omega_bar_j[,merge_cand])
          v_bar_j[,merger] <- apply(v_bar_j[,merge_cand], 1, prod) / u_bar_j
        }
        
        J <- length(keep)
        if (J_auto) pi_j <- pi_j[,keep, drop = FALSE]
        if (J_auto) log_pi_j <- log_pi_j[,keep, drop = FALSE]
        pi_bar_j <- pi_bar_j[,keep, drop = FALSE]
        log_pi_bar_j <- log_pi_bar_j[,keep, drop = FALSE]
        b_bar_j <- b_bar_j[,keep, drop = FALSE]
        omega_bar_j <- omega_bar_j[,keep, drop = FALSE]
        v_bar_j <- v_bar_j[,keep, drop = FALSE]
      }
    }
    if (verbose & !merged) print("Merging components.")
  }
  
  # if components were merged out use auto procedure to increase to desired JLK
  merge_flag <- (J < J_max & !J_auto) | (L < L_max & !L_auto) | (K < K_max & !K_auto)
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components
  
  #### auto procedures ####
  #### auto procedure with single component ####
  last_restart <- ifelse(restart, 2, Inf)
  
  if (sum(J_auto, L_auto, K_auto) == 1 | merge_flag) {
    
    refit <- (L > 0 | K > 0 | J > 0)
    counter <- n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo
    merged_at_l <- c(); merged_at_k <- c(); merged_at_j <- c();
    
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
      
      # fit new model and merge ####
      merged <- FALSE
      while (!merged) {
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
        
        merged <- TRUE
        refit <- TRUE
        
        if (L > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = L, nrow = L)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities and keep_mat
          merge_mat <- (t(pi_bar_l) %*% pi_bar_l > merge_prob) & keep_mat
          
          if (any(merge_mat[lower.tri(merge_mat)]) & !(L %in% merged_at_l)) {
            merged <- FALSE
            
            # keep track of where merges occur
            merged_at_l <- c(merged_at_l, L)
            
            # force transitivity
            merge_mat <- force_adj(merge_mat)
            diag(merge_mat) <- TRUE
            
            # combine and reset components
            keep <- which(apply(merge_mat, 2, sum) == 1)
            visited <- keep
            for (i in 1:L) {
              if (i %in% visited) next
              merge_cand <- which(merge_mat[,i])
              merger <- merge_cand[which.max(apply(pi_bar_l[,merge_cand], 2, max))]
              mergees <- merge_cand[merge_cand != merger]
              keep <- c(merger, keep)
              visited <- c(merge_cand , visited)

              b_bar_l[,merger] <- rowSums(b_bar_l[,merge_cand])
              omega_bar_l[,merger] <- rowMeans(omega_bar_l[,merge_cand])
            }
            
            counter <- counter + L - length(keep)
            
            L <- length(keep)
            if (L_auto) pi_l <- pi_l[,keep, drop = FALSE]
            if (L_auto) log_pi_l <- log_pi_l[,keep, drop = FALSE]
            pi_bar_l <- pi_bar_l[,keep, drop = FALSE]
            log_pi_bar_l <- log_pi_bar_l[,keep, drop = FALSE]
            b_bar_l <- b_bar_l[,keep, drop = FALSE]
            omega_bar_l <- omega_bar_l[,keep, drop = FALSE]
          }
        }
        
        if (K > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_k, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = K, nrow = K)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities 
          merge_mat <- (t(pi_bar_k) %*% pi_bar_k > merge_prob) & keep_mat
          
          if (any(merge_mat[lower.tri(merge_mat)]) & !(K %in% merged_at_k)) {
            merged <- FALSE
            
            # keep track of where merges occur
            merged_at_k <- c(merged_at_k, K)
            
            # force transitivity
            merge_mat <- force_adj(merge_mat)
            diag(merge_mat) <- TRUE
            
            # combine and reset components
            keep <- which(apply(merge_mat, 2, sum) == 1)
            visited <- keep
            for (i in 1:K) {
              if (i %in% visited) next
              merge_cand <- which(merge_mat[,i])
              merger <- merge_cand[which.max(apply(pi_bar_k[,merge_cand], 2, max))]
              mergees <- merge_cand[merge_cand != merger]
              keep <- c(merger, keep)
              visited <- c(merge_cand , visited)
              
              v_bar_k[,merger] <- apply(v_bar_k[,merge_cand], 1, prod) / u_bar_k
              v_bar_k[,mergees] <- u_bar_k
            }
            
            counter <- counter + K - length(keep)
            K <- length(keep)
            if (K_auto) pi_k <- pi_k[,keep, drop = FALSE]
            if (K_auto) log_pi_k <- log_pi_k[,keep, drop = FALSE]
            pi_bar_k <- pi_bar_k[,keep, drop = FALSE]
            log_pi_bar_k <- log_pi_bar_k[,keep, drop = FALSE]
            v_bar_k <- v_bar_k[,keep, drop = FALSE]
          }
        }
        
        if (J > 1) {
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_j, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = J, nrow = J)
          keep_mat[keep, keep] <- TRUE
          
          # compute adjacency matrix from pairwise merge probabilities 
          merge_mat <- (t(pi_bar_j) %*% pi_bar_j > merge_prob) & keep_mat
          
          if (any(merge_mat[lower.tri(merge_mat)]) & !(J %in% merged_at_j)) {
            merged <- FALSE
            
            # keep track of where merges occur
            merged_at_j <- c(merged_at_j, J)
            
            # force transitivity
            merge_mat <- force_adj(merge_mat)
            diag(merge_mat) <- TRUE
            
            # combine and reset components
            keep <- which(apply(merge_mat, 2, sum) == 1)
            visited <- keep            
            for (i in 1:J) {
              if (i %in% visited) next
              merge_cand <- which(merge_mat[,i])
              merger <- merge_cand[which.max(apply(pi_bar_j[,merge_cand], 2, max))]
              mergees <- merge_cand[merge_cand != merger]
              keep <- c(merger, keep)
              visited <- c(merge_cand , visited)
              
              b_bar_j[,merger] <- rowSums(b_bar_j[,merge_cand])
              omega_bar_j[,merger] <- rowMeans(omega_bar_j[,merge_cand])
              v_bar_j[,merger] <- apply(v_bar_j[,merge_cand], 1, prod) / u_bar_j
            }
            
            counter <- counter + J - length(keep)
            J <- length(keep)
            if (J_auto) pi_j <- pi_j[,keep, drop = FALSE]
            if (J_auto) log_pi_j <- log_pi_j[,keep, drop = FALSE]
            pi_bar_j <- pi_bar_j[,keep, drop = FALSE]
            log_pi_bar_j <- log_pi_bar_j[,keep, drop = FALSE]
            b_bar_j <- b_bar_j[,keep, drop = FALSE]
            omega_bar_j <- omega_bar_j[,keep, drop = FALSE]
            v_bar_j <- v_bar_j[,keep, drop = FALSE]
          }
        }
      }
      
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
        pi_bar_j <- matrix(1/T, nrow = T, ncol = J)
        log_pi_bar_j <- matrix(0, nrow = T, ncol = J)
        b_bar_j <- matrix(0.0, nrow = T, ncol = J)
        omega_bar_j <- matrix(1.0, nrow = T, ncol = J)
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
    }
    if (!pi_k_weighted) {
      pi_k <- pi_k[T:1,1:max(1,K),drop = FALSE]
      log_pi_k <- log_pi_k[T:1,1:max(1,K),drop = FALSE]
    }
    if (!pi_j_weighted) {
      pi_j <- pi_j[T:1,1:max(1,J),drop = FALSE]
      log_pi_j <- log_pi_j[T:1,1:max(1,J),drop = FALSE]
    }
    
    # reversing mean components
    L_seq = seq_len(L)
    pi_bar_l <- pi_bar_l[, L_seq, drop = FALSE]
    log_pi_bar_l <- log_pi_bar_l[, L_seq, drop = FALSE]
    b_bar_l <- b_bar_l[, L_seq, drop = FALSE]
    omega_bar_l <- omega_bar_l[, L_seq, drop = FALSE]
    
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
    pi_bar_k <- pi_bar_k[, K_seq, drop = FALSE]
    log_pi_bar_k <- log_pi_bar_k[, K_seq, drop = FALSE]
    v_bar_k <- v_bar_k[, K_seq, drop = FALSE]
    
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