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
  pi_ij <- matrix(1 / T, nrow = T, ncol = L)
  mu_ij <- matrix(0, nrow = T, ncol = L)
  sigma2_ij <- rep(1, length = T)
  beta <- mu_ij * pi_ij
  
  omega_ij <- matrix(1 / T, nrow = T, ncol = K)
  u_ij <- u_0 + (T - 1:T + 1) / 2
  v_ij <- matrix(u_ij, nrow = T, ncol = K)
  lambda2_k <- matrix(0, nrow = T, ncol = K)
  
  for (k in 1:K) {
    lambda2_k[,k] <- lambda_update(u_ij, v_ij[,k], omega_ij[,k])
  }
  
  lambda2 <- apply(lambda2_k, 1, prod)
  
  # store current parameter values
  params <- c(sigma2_ij, mu_ij, pi_ij, v_ij, omega_ij)
  
  # initialize residual
  r_bar <- y - rowSums(apply(beta, 2, cumsum))
  
  n_iter <- 1
  
  while (TRUE) {
    
    new_params <- c()
    sigma2_ij <- 1 / (revcumsum(lambda2) / sigma2 + 1 / sigma2_0)
    
    new_params <- c(new_params, sigma2_ij)
    
    # mean update step
    for (l in 1:L) {
      
      r_bar <- r_bar + cumsum(beta[,l])
      
      mu_ij[,l] <- (sigma2_ij / sigma2) * revcumsum(lambda2 * r_bar)
      
      pi_ij[,l] <- log(pi) + 0.5 * log(sigma2_ij) + mu_ij[,l]^2 / (2 * sigma2_ij)
      pi_ij[,l] <- prop.table(exp(pi_ij[,l] - max(pi_ij[,l])))
      
      beta[,l] <- mu_ij[,l] * pi_ij[,l]
      
      r_bar <- r_bar - cumsum(beta[,l])
    }
    
    new_params <- c(new_params, mu_ij, pi_ij)
    
    # precision update step
    Xb2 <- rowSums(apply(beta,2,cumsum)^2)
    
    for (k in 1:K) {
      lambda2 <- lambda2 / lambda2_k[,k]
      
      z2_tilde <- lambda2 * (r_bar^2 - Xb2 + cumsum(rowSums((mu_ij^2 + sigma2_ij) * pi_ij)))
      
      v_ij[,k] <- v_0 + revcumsum(z2_tilde) / (2 * sigma2)
      
      omega_ij[,k] <- log(omega) + lgamma(u_ij) - u_ij * log(v_ij[,k]) - cumsum(c(0,z2_tilde))[-(T+1)] / (2 * sigma2)
      omega_ij[,k] <- prop.table(exp(omega_ij[,k] - max(omega_ij[,k])))
      
      lambda2_k[,k] <- lambda_update(u_ij, v_ij[,k], omega_ij[,k])
      
      lambda2 <- lambda2 * lambda2_k[,k]
    }
    
    new_params <- c(new_params, v_ij, omega_ij)
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 100 == 0) print(paste("Iteration", n_iter))
    n_iter <- n_iter + 1
    params <- new_params
  }
  
  return(list(gamma = pi_ij, mu = mu_ij, sigma = sigma2_ij, 
              alpha = omega_ij, v = v_ij, u = u_ij,
              beta = beta, lambda2 = lambda2))
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

