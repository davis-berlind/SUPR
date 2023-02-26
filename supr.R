# helper functions

revcumsum <- function(x){
  rev(cumsum(rev(x)))
}

lambda_update <- function(u, v, prob) {
  cumsum((u / v) * prob) + c(revcumsum(prob[-1]), 0)
}

# main function

supr <- function(y, sigma2 = 1, sigma2_0, K, L, tol = 1e-5, u_0 = 1e-3, v_0 = 1e-3) {
  
  # model parameters
  T <- length(y)
  
  # uniform prior on change point locations
  pi <- 1 / T    
  omega <- 1 / T
  
  # initializing posterior parameters
  
  gamma_post <- matrix(1 / T, nrow = T, ncol = L)
  alpha_post <- matrix(1 / T, nrow = T, ncol = K)
  
  mu_post <- matrix(0, nrow = T, ncol = L)
  sigma_post <- rep(1, length = T)
  beta <- mu_post * gamma_post
  u_post <- u_0 + (T - 1:T + 1) / 2
  v_post <- matrix(u_post, nrow = T, ncol = K)
  
  lambda2_k <- matrix(0, nrow = T, ncol = K)
  
  for (k in 1:K) {
    lambda2_k[,k] <- lambda_update(u_post, v_post[,k], alpha_post[,k])
  }
  
  lambda2 <- apply(lambda2_k, 1, prod)
  
  # store current parameter values
  params <- c(sigma_post, mu_post, gamma_post, v_post, alpha_post)
  
  # initialize residual
  r_bar <- y - rowSums(apply(beta,2,cumsum))
  
  n_iter <- 1
  
  while (TRUE) {
    
    new_params <- c()
    sigma_post <- 1 / (revcumsum(lambda2) / sigma2 + 1 / sigma2_0)
    
    new_params <- c(new_params, sigma_post)
    
    # mean update step
    for (l in 1:L) {
      
      r_bar <- r_bar + cumsum(beta[,l])
      
      mu_post[,l] <- (sigma_post / sigma2) * revcumsum(lambda2 * r_bar)
      
      gamma_post[,l] <- log(pi) + 0.5 * log(sigma_post) + mu_post[,l]^2 / (2 * sigma_post)
      gamma_post[,l] <- prop.table(exp(gamma_post[,l] - max(gamma_post[,l])))
      
      beta[,l] <- mu_post[,l] * gamma_post[,l]
      
      r_bar <- r_bar - cumsum(beta[,l])
    }
    
    new_params <- c(new_params, mu_post, gamma_post)
    
    # precision update step
    Xb2 <- rowSums(apply(beta,2,cumsum)^2)
    
    for (k in 1:K) {
      lambda2 <- lambda2 / lambda2_k[,k]
      
      z2_tilde <- lambda2 * (r_bar^2 - Xb2 + cumsum(rowSums((mu_post^2 + sigma_post) * gamma_post)))
      
      v_post[,k] <- v_0 + revcumsum(z2_tilde) / (2 * sigma2)
      
      alpha_post[,k] <- log(omega) + lgamma(u_post) - u_post * log(v_post[,k]) - cumsum(c(0,z2_tilde))[-(T+1)] / (2 * sigma2)
      alpha_post[,k] <- prop.table(exp(alpha_post[,k] - max(alpha_post[,k])))
      
      lambda2_k[,k] <- lambda_update(u_post, v_post[,k], alpha_post[,k])
      
      lambda2 <- lambda2 * lambda2_k[,k]
    }
    
    new_params <- c(new_params, v_post, alpha_post)
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 100 == 0) print(paste("Iteration", n_iter))
    n_iter <- n_iter + 1
    params <- new_params
  }
  return(list(gamma = gamma_post, mu = mu_post, sigma = sigma_post, 
              alpha = alpha_post, v = v_post, u = u_post,
              beta = beta, lambda2 = lambda2))
}

cred_set <- function(probs, level = 0.9) {
  L <- ncol(probs) 
  cs <- list()
  for(l in 1:L) {
    set <- order(probs[,l],decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,l], decreasing = TRUE), l]) > level)]
    cs[[l]] <- set
  }
  return(cs)
}

cred_set(alpha_post)
cred_set(gamma_post)
