#### helper functions #### 

prod_l <- function(l, mat) {
  apply(mat[,-l], 1, prod)
}

revcumsum <- function(x){
  rev(cumsum(rev(x)))
}

lambda_update <- function(u, v, prob) {
  cumsum((u / v) * prob) + c(revcumsum(prob[-1]), 0)
}

#### Main #### 

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
  r_bar <- y - rowSums(apply(beta, 2, cumsum))
  
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

supr_sim <- function(y, sigma2 = 1, sigma2_0, L, tol = 1e-5, u_0 = 1e-3, v_0 = 1e-3) {
  
  # model parameters
  T <- length(y)
  
  # uniform prior on change point locations
  pi <- 1 / T    

  # initializing posterior parameters
  
  # gamma_post <- matrix(1 / T, nrow = T, ncol = L)
  gamma_post <- matrix(0, nrow = T, ncol = L)
  gamma_post[201,1] <- 1
  gamma_post[401,2] <- 1
  gamma_post[,-c(1,2)] <- 1/T

  mu_post <- matrix(0, nrow = T, ncol = L)
  sigma_post <- matrix(1, nrow = T, ncol = L)
  # beta <- mu_post * gamma_post
  u_post <- u_0 + (T - 1:T + 1) / 2
  v_post <- matrix(u_post, nrow = T, ncol = L)
  
  lambda2_l <- matrix(0, nrow = T, ncol = L)

  for (l in 1:L) {
    lambda2_l[,l] <- lambda_update(u_post, v_post[,l], gamma_post[,l])
  }
  
  lambda2 <- apply(lambda2_l, 1, prod)

  beta_lambda <- mu_post * (u_post / v_post) * gamma_post
  beta_lambda_cs <- apply(beta_lambda, 2, cumsum)
  beta2_lambda <- (mu_post^2 + sigma_post) * (u_post / v_post) * gamma_post

  # store current parameter values
  params <- c(sigma_post, mu_post, gamma_post, v_post)
  
  n_iter <- 1
  
  while (TRUE) {
    new_params <- c()
    
    for (l in 1:L) {
      
      lambda2 <- lambda2 / lambda2_l[,l]
      lambda_minus_lj <- lambda2 / lambda2_l[,-l]
      
      xbl <- rowSums(beta_lambda_cs[,-l] * lambda_minus_lj)
      
      # updating q(b_i | gamma_i)
      
      s2_lj <- u_post / v_post[,l]

      sigma_post[,l] <- 1 / (s2_lj * revcumsum(lambda2) / sigma2 + 1 / sigma2_0)
      mu_post[,l] <- (s2_lj * sigma_post[,l] / sigma2) * revcumsum(y * lambda2 - xbl)
      
      # updating q(s^2_i | gamma_i)

      v <- revcumsum(lambda2 * y^2)
      v_f <- cumsum(lambda2 * y^2)

      v <- v + (mu_post[,l]^2 + sigma_post[,l]) * revcumsum(lambda2)
      v_f <- v_f + (mu_post[,l]^2 + sigma_post[,l]) * cumsum(lambda2)

      v <- v - 2 * mu_post[,l] * revcumsum(lambda2 * y)
      v_f <- v_f - 2 * mu_post[,l] * cumsum(lambda2 * y)

      v <- v + 2 * (mu_post[,l] * revcumsum(xbl) - revcumsum(y * xbl))
      v_f <- v_f + 2 * (mu_post[,l] * cumsum(xbl) - cumsum(y * xbl))

      v <- v + revcumsum(rowSums(apply(beta2_lambda[,-l], 2, cumsum) * lambda_minus_lj))
      v_f <- v_f + cumsum(rowSums(apply(beta2_lambda[,-l], 2, cumsum) * lambda_minus_lj))

      bxxb <- rep(0, length(v))

      for (i in 1:(L-2)) {
        for (j in (i+1):(L-1)) {
          bxxb <- bxxb + 2 * beta_lambda_cs[,-l][,i] * beta_lambda_cs[,-l][,j] * lambda_minus_lj[,i] / lambda2_l[,-l][,j]
        }
      }

      v <- v + revcumsum(bxxb)
      v_f <- v_f + cumsum(bxxb)

      v_post[,l] <- v / (2 * sigma2) + v_0

      # updating q(gamma_i)

      # gamma_post[,l] <- log(pi) + lgamma(u_post) - u_post * log(v_post[,l]) + 0.5 * log(sigma_post[,l]) - cumsum(c(0, v_f))[-(T+1)] / (2 * sigma2)
      # gamma_post[,l] <- prop.table(exp(gamma_post[,l] - max(gamma_post[,l])))

      # intermediate updates
      
      lambda2_l[,l] <- lambda_update(u_post, v_post[,l], gamma_post[,l])
      lambda2 <- lambda2 * lambda2_l[,l]
      
      beta_lambda[,l] <- mu_post[,l] * (u_post / v_post[,l]) * gamma_post[,l]

      beta_lambda_cs[,l] <- cumsum(beta_lambda[,l])
      beta2_lambda[,l] <- (mu_post[,l]^2 + sigma_post[,l]) * (u_post / v_post[,l]) * gamma_post[,l]

    }
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 100 == 0) print(paste("Iteration", n_iter))
    n_iter <- n_iter + 1
    params <- new_params
  }
  
  return(list(gamma = gamma_post, mu = mu_post, sigma = sigma_post, 
              v = v_post, u = u_post,
              beta = mu_post * gamma_post, lambda2 = lambda2))
}

v_update <- function(y, lambda2, lambda2_l, beta_lambda, beta2_lambda, mu_l, sigma2_l, xbl){
  L <- ncol(lambda2_l)
  v <- revcumsum(lambda2 * y^2) + (mu_l^2 + sigma2_l) * revcumsum(lambda2) - 2 * mu_l * revcumsum(lambda2 * y)
  v <- v + 2 * (mu_l * revcumsum(xbl) - revcumsum(y * xbl))
  v <- v + revcumsum(rowSums(apply(beta2_lambda[,-l], 2, cumsum) * sapply(1:(L-1), prod_l, lambda2_l[,-l])))
  beta_lambda_cs <- apply(beta_lambda, 2, cumsum)
  bxxb <- rep(0, length(v))
  for (i in 1:L) {
    for (j in 1:L) {
      if (i == j) {
        next
      } else {
        bxxb <- bxxb + beta_lambda_cs[,i] * beta_lambda_cs[,j] * prod_l(c(i,j), lambda2_l)
      }
    }
  }
  return(v + revcumsum(bxxb))
}


#### Post Processing #### 

cred_set_susie <- function(probs, level = 0.95) {
  L <- ncol(probs)
  T <- nrow(probs)
  probs <- probs[-T,]
  cs <- list()
  for(l in 1:L) {
    set <- order(probs[,l], decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,l], decreasing = TRUE), l]) > level)]
    if(length(set) == 1) {
      cs <- c(cs, list(set))
      next
    } 
    cmbs <- combn(set, 2)
    i <- cmbs[1,]
    j <- cmbs[2,]
    purity <- min(exp(log(abs(T * min(i,j) - j * i)) - 0.5 * (log(j) + log(T - j) + log(i) + log(T - i))))
    if (purity > 0.5) {
      cs <- c(cs,list(set))
    }
  }
  
  if (length(cs) == 0) return(cs)
  
  cs <- unique(cs)
  nset <- length(cs)
  subset <- c()
  for (i in 1:nset) {
    for(j in 1:nset) {
       if (i == j) next
       if (length(setdiff(cs[[i]], cs[[j]])) == 0) subset <- c(subset, i)
    }
  }
  
  if (length(subset) > 0) cs[[subset]] <- NULL
  
  return(cs)
}

cred_set_prisca <- function(probs, level = 0.95) {
  L <- ncol(probs)
  T <- nrow(probs)
  cs <- list()
  for(l in 1:L) {
    set <- order(probs[,l], decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,l], decreasing = TRUE), l]) > level)]
    if (length(set) <= T / 2) {
      cs <- c(cs,list(set))
    }
  }
  
  if (length(cs) == 0) return(cs)
  
  cs <- unique(cs)
  nset <- length(cs)
  subset <- c()
  for (i in 1:nset) {
    for(j in 1:nset) {
      if (i == j) next
      if (length(setdiff(cs[[i]], cs[[j]])) == 0) subset <- c(subset, i)
    }
  }
  
  if (length(subset) > 0) cs[[subset]] <- NULL
  
  return(cs)
}

cred_set_prisca(supr_mod$alpha)
cred_set_susie(supr_mod$gamma)
cred_set_prisca(supr_mod$gamma)

