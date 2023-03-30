supr_sim <- function(y, sigma2 = 1, sigma2_0, L, tol = 1e-5, u_0 = 1e-3, v_0 = 1e-3) {
  
  # model parameters
  T <- length(y)
  
  # uniform prior on change point locations
  pi <- 1 / T    
  
  # initializing posterior parameters
  
  # gamma_post <- matrix(1 / T, nrow = T, ncol = L)
  gamma_ij <- matrix(0, nrow = T, ncol = L)
  gamma_ij[201,1] <- 1
  gamma_ij[401,2] <- 1
  gamma_ij[,-c(1,2)] <- 1/T
  
  mu_ij <- matrix(0, nrow = T, ncol = L)
  sigma_ij <- matrix(1, nrow = T, ncol = L)
  # beta <- mu_post * gamma_post
  u_ij <- u_0 + (T - 1:T + 1) / 2
  v_ij <- matrix(u_ij, nrow = T, ncol = L)
  
  lambda2_l <- matrix(0, nrow = T, ncol = L)
  
  for (l in 1:L) {
    lambda2_l[,l] <- lambda_update(u_ij, v_ij[,l], gamma_ij[,l])
  }
  
  lambda2 <- apply(lambda2_l, 1, prod)
  
  beta_lambda <- mu_ij * (u_ij / v_ij) * gamma_ij
  beta_lambda_cs <- apply(beta_lambda, 2, cumsum)
  beta2_lambda <- (mu_ij^2 * (u_ij / v_ij) + sigma_ij) * gamma_ij
  beta2_lambda_cs <- apply(beta2_lambda, 2, cumsum)
  
  # store current parameter values
  params <- c(sigma_ij, mu_ij, gamma_ij, v_ij)
  
  n_iter <- 1
  
  while (TRUE) {
    
    new_params <- c()
    
    for (i in 1:L) {
      
      lambda2 <- lambda2 / lambda2_l[,i]
      lambda_minus_ij <- lambda2 / lambda2_l[,-i]
      
      xbl <- rowSums(beta_lambda_cs[,-i] * lambda_minus_ij)
      
      # updating q(b_i | gamma_i)
      
      sigma_ij[,i] <- 1 / (revcumsum(lambda2) / sigma2 + 1 / sigma2_0)
      mu_ij[,i] <- (sigma_ij[,i] / sigma2) * revcumsum(y * lambda2 - xbl)
      
      # updating q(s^2_i | gamma_i)
      
      bxxb <- rowSums(beta2_lambda_cs[,-i] * lambda_minus_ij)
      
      for (l in 1:(L-2)) {
        for (j in (l+1):(L-1)) {
          bxxb <- bxxb + 2 * beta_lambda_cs[,-i][,l] * beta_lambda_cs[,-i][,j] * lambda_minus_ij[,l] / lambda2_l[,-i][,j]
        }
      }
      
      v_ij[,i] <-  v_0 - mu_ij[,i]^2 / (2 * sigma_ij[,i]) + revcumsum(lambda2 * y^2 + bxxb - 2 * y * xbl) / (2 * sigma2)
      
      # updating q(gamma_i)
      
      # gamma_post[,l] <- log(pi) + lgamma(u_post) - u_post * log(v_post[,l]) + 0.5 * log(sigma_post[,l]) - cumsum(c(0, v_f))[-(T+1)] / (2 * sigma2)
      # gamma_post[,l] <- prop.table(exp(gamma_post[,l] - max(gamma_post[,l])))
      
      # intermediate updates
      
      lambda2_l[,i] <- lambda_update(u_ij, v_ij[,i], gamma_ij[,i])
      lambda2 <- lambda2 * lambda2_l[,i]
      
      beta_lambda[,i] <- mu_ij[,i] * (u_ij / v_ij[,i]) * gamma_ij[,i]
      beta_lambda_cs[,i] <- cumsum(beta_lambda[,i])
      beta2_lambda[,i] <- (mu_ij[,i]^2 + sigma_ij[,i]) * (u_ij / v_ij[,i]) * gamma_ij[,i]
      beta2_lambda_cs[,i] <-  cumsum(beta2_lambda[,i])
    }
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 100 == 0) print(paste("Iteration", n_iter))
    n_iter <- n_iter + 1
    params <- new_params
  }
  
  return(list(gamma = gamma_ij, mu = mu_ij, sigma = sigma_ij, 
              v = v_ij, u = u_ij,
              beta = mu_ij * gamma_ij, lambda2 = lambda2))
}