supr_sim_single <- function(y, sigma2, sigma2_0, u_0, v_0) {
  
  # model parameters
  T <- length(y)
  
  # uniform prior on change point locations
  pi <- 1 / T    
  
  # posterior parameters
  
  # updating q(b_i | s_i^2, gamma_i)
  sigma2_ij <- 1 / ((T - 1:T + 1) / sigma2 + 1 / sigma2_0)
  mu_ij <- (sigma2_ij / sigma2) * revcumsum(y)
    
  # updating q(s^2_i | gamma_i)
  u_ij <- u_0 + (T - 1:T + 1) / 2
  v_ij <- v_0 - mu_ij^2 / (2 * sigma2_ij) + revcumsum(y^2) / (2 * sigma2)
  
  # updating q(gamma_i)
  log_C_ij <- log(pi) + 0.5 * log(sigma2_ij) + lgamma(u_ij) - u_ij * log(v_ij) - c(0, cumsum(y^2))[-(T+1)] / (2 * sigma2)
  pi_ij <- prop.table(exp(log_C_ij - max(log_C_ij)))
  
  return(list(pi = pi_ij, mu = mu_ij, sigma = sigma2_ij, 
              v = v_ij, u = u_ij,
              beta = mu_ij * pi_ij, 
              lambda2 = lambda_update(u_ij, v_ij, pi_ij)))
}
  
supr_sim <- function(y, sigma2 = 1, sigma2_0, L, tol = 1e-5, u_0 = 1e-3, v_0 = 1e-3) {
  
  # return exact posterior when L = 1
  if (L == 1) return(supr_sim_single(y, sigma2, sigma2_0, u_0, v_0))
  
  # model parameters
  T <- length(y)
  
  # uniform prior on change point locations
  pi <- 1 / T    
  
  # initializing posterior parameters
  pi_ij <- matrix(1 / T, nrow = T, ncol = L)
  mu_ij <- matrix(0, nrow = T, ncol = L)
  sigma2_ij <- matrix(1, nrow = T, ncol = L)
  u_ij <- u_0 + (T - 1:T + 1) / 2
  v_ij <- matrix(u_ij, nrow = T, ncol = L)
  
  # intermediate calculations
  lambda2_l <- matrix(0, nrow = T, ncol = L)
  
  for (l in 1:L) {
    lambda2_l[,l] <- lambda_update(u_ij, v_ij[,l], pi_ij[,l])
  }
  
  lambda2 <- apply(lambda2_l, 1, prod)
  
  beta_lambda2 <- mu_ij * (u_ij / v_ij) * pi_ij
  beta_lambda2_cs <- apply(beta_lambda2, 2, cumsum)
  beta2_lambda2 <- (mu_ij^2 * (u_ij / v_ij) + sigma2_ij) * pi_ij
  beta2_lambda2_cs <- apply(beta2_lambda2, 2, cumsum)
  
  # store current parameter values
  params <- c(sigma2_ij, mu_ij, pi_ij, v_ij)
  
  n_iter <- 1
  
  # vb routine
  while (TRUE) {
    
    for (i in 1:L) {
      
      lambda2 <- lambda2 / lambda2_l[,i]           # factor out ith precision
      lambda2_minus_ij <- lambda2 / lambda2_l[,-i] # factor out jth precision, j =/= i 
      
      # calculate \sum_{k \neq i} x'_t E[\beta_k \lambda^2_{k,t}] \prod_{\ell \neq k,i} E[\lambda_{\ell,t}^2]
      if (L == 2) xbl <- beta_lambda2_cs[,-i] * lambda2_minus_ij
      else xbl <- rowSums(beta_lambda2_cs[,-i] * lambda2_minus_ij)
      
      # updating q(b_i | s_i^2, gamma_i)
      sigma2_ij[,i] <- 1 / (revcumsum(lambda2) / sigma2 + 1 / sigma2_0)
      mu_ij[,i] <- (sigma2_ij[,i] / sigma2) * revcumsum(y * lambda2 - xbl)
      
      # updating q(s^2_i | gamma_i)
      # calculate sum_{k \neq i} sum_{k' \neq i} E[\beta'_{k'} x_t x'_t \beta_k prod_{\ell \neq i} \lambda_{\ell,t}^2]
      if (L == 2) bxxb <- beta2_lambda2_cs[,-i] * lambda2_minus_ij    # when L = 2, taking out ith column turns everything to a vector
      else bxxb <- rowSums(beta2_lambda2_cs[,-i] * lambda2_minus_ij)
      
      # note there are only cross-product terms when  L > 2
      if (L > 2) {
        for (l in 1:(L-2)) {
          for (j in (l+1):(L-1)) {
            bxxb <- bxxb + 2 * beta_lambda2_cs[,-i][,l] * beta_lambda2_cs[,-i][,j] * lambda2_minus_ij[,l] / lambda2_l[,-i][,j]
          }
        }
      }
      
      v_ij[,i] <- v_0 - mu_ij[,i]^2 / (2 * sigma2_ij[,i]) + revcumsum(lambda2 * y^2 + bxxb - 2 * y * xbl) / (2 * sigma2)
      
      # updating q(gamma_i)
      log_C_ij <- log(pi) + 0.5 * log(sigma2_ij[,i]) + lgamma(u_ij) - u_ij * log(v_ij[,i]) - c(0, cumsum(lambda2 * y^2 + bxxb - 2 * y * xbl))[-(T+1)] / (2 * sigma2)
      pi_ij[,i] <- prop.table(exp(log_C_ij - max(log_C_ij)))
      
      # intermediate updates
      lambda2_l[,i] <- lambda_update(u_ij, v_ij[,i], pi_ij[,i])
      lambda2 <- lambda2 * lambda2_l[,i]
      
      beta_lambda2[,i] <- mu_ij[,i] * (u_ij / v_ij[,i]) * pi_ij[,i]
      beta_lambda2_cs[,i] <- cumsum(beta_lambda2[,i])
      
      beta2_lambda2[,i] <- (mu_ij[,i]^2 * (u_ij / v_ij[,i]) + sigma2_ij[,i]) * pi_ij[,i]
      beta2_lambda2_cs[,i] <- cumsum(beta2_lambda2[,i])
      
    }
    
    new_params <- c(sigma2_ij, mu_ij, pi_ij, v_ij)
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 1000 == 0) print(paste("Iteration", n_iter))
    n_iter <- n_iter + 1
    params <- new_params
  }
  
  return(list(pi = pi_ij, mu = mu_ij, sigma2 = sigma2_ij, 
              v = v_ij, u = u_ij,
              beta = mu_ij * pi_ij, lambda2 = lambda2))
}
