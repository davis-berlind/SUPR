supr_sim_single <- function(y, sigma2, sigma2_0, u_0, v_0) {
  
  # model parameters
  T <- length(y)
  
  # uniform prior on change point locations
  pi <- 1 / T    
  
  # posterior parameters
  
  # location component
  sigma2_ij <- 1 / ((T - 1:T + 1) / sigma2 + 1 / sigma2_0)
  mu_ij <- (sigma2_ij / sigma2) * revcumsum(y)
    
  # scale component
  u_ij <- u_0 + (T - 1:T + 1) / 2
  v_ij <- v_0 - mu_ij^2 / (2 * sigma2_ij) + revcumsum(y^2) / (2 * sigma2)
  
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
  
  # gamma_post <- matrix(1 / T, nrow = T, ncol = L)
  pi_ij <- matrix(0, nrow = T, ncol = L)
  pi_ij[201,1] <- 1
  pi_ij[401,2] <- 1
  pi_ij[,-c(1,2)] <- 1/T
  
  mu_ij <- matrix(0, nrow = T, ncol = L)
  sigma2_ij <- matrix(1, nrow = T, ncol = L)
  # beta <- mu_post * gamma_post
  u_ij <- u_0 + (T - 1:T + 1) / 2
  v_ij <- matrix(u_ij, nrow = T, ncol = L)
  
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
  
  while (TRUE) {
    
    new_params <- c()
    
    for (i in 1:L) {
      
      lambda2 <- lambda2 / lambda2_l[,i]
      lambda2_minus_ij <- lambda2 / lambda2_l[,-i]
      
      if (L == 2) xbl <- beta_lambda2_cs[,-i] * lambda2_minus_ij
      else xbl <- rowSums(beta_lambda2_cs[,-i] * lambda2_minus_ij)
      
      # updating q(b_i | gamma_i)
      
      sigma2_ij[,i] <- 1 / (revcumsum(lambda2) / sigma2 + 1 / sigma2_0)
      mu_ij[,i] <- (sigma2_ij[,i] / sigma2) * revcumsum(y * lambda2 - xbl)
      
      # updating q(s^2_i | gamma_i)
      
      # the problem is the next couple of lines
      
      if (L == 2) bxxb <- beta2_lambda2_cs[,-i] * lambda2_minus_ij
      else bxxb <- rowSums(beta2_lambda2_cs[,-i] * lambda2_minus_ij)
      
      if (L > 2) {
        for (l in 1:(L-2)) {
          for (j in (l+1):(L-1)) {
            bxxb <- bxxb + 2 * beta_lambda2_cs[,-i][,l] * beta_lambda2_cs[,-i][,j] * lambda2_minus_ij[,l] / lambda2_l[,-i][,j]
          }
        }
      }
      
      v_ij[,i] <- v_0 - mu_ij[,i]^2 / (2 * sigma2_ij[,i]) + revcumsum(lambda2 * y^2 + bxxb - 2 * y * xbl) / (2 * sigma2)
      
      # updating q(gamma_i)
      
      # gamma_post[,l] <- log(pi) + lgamma(u_post) - u_post * log(v_post[,l]) + 0.5 * log(sigma_post[,l]) - cumsum(c(0, v_f))[-(T+1)] / (2 * sigma2)
      # gamma_post[,l] <- prop.table(exp(gamma_post[,l] - max(gamma_post[,l])))
      
      # intermediate updates
      
      lambda2_l[,i] <- lambda_update(u_ij, v_ij[,i], pi_ij[,i])
      lambda2 <- lambda2 * lambda2_l[,i]
      
      beta_lambda2[,i] <- mu_ij[,i] * (u_ij / v_ij[,i]) * pi_ij[,i]
      beta_lambda2_cs[,i] <- cumsum(beta_lambda2[,i])
      beta2_lambda2[,i] <- (mu_ij[,i]^2 * (u_ij / v_ij[,i]) + sigma2_ij[,i]) * pi_ij[,i]
      beta2_lambda2_cs[,i] <- cumsum(beta2_lambda2[,i])
    }
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 100 == 0) print(paste("Iteration", n_iter))
    n_iter <- n_iter + 1
    params <- new_params
  }
  
  return(list(gamma = pi_ij, mu = mu_ij, sigma = sigma2_ij, 
              v = v_ij, u = u_ij,
              beta = mu_ij * pi_ij, lambda2 = lambda2))
}
