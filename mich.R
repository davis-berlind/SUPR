#### MICH I ####

mich_i <- function(y, L, K, tol = 1e-5, fit.intercept = TRUE,
                  tau = 0.1, u = 1e-3, v = 1e-3, tau_0 = 0.1, u_0 = 1e-3, v_0 = 1e-3,
                  pi = NULL, omega = NULL) {
  
  if (!is.numeric(y) | !is.vector(y)) stop("y must be a numeric vector.")
  
  # model parameters
  if (fit.intercept) {
    y_0 <- y[1]
    y <- y[-1]
  }
  
  T <- length(y)
  
  # checking that priors are proper
  tau <- prior_check(tau, L)
  if (is.null(tau)) stop("tau must either be a positive number or length L vector of positive number.")
  u <- prior_check(u, K)
  if (is.null(u)) stop("u must either be a positive number or length K vector of positive number.")
  v <- prior_check(v, K)
  if (is.null(v)) stop("v must either be a positive number or length K vector of positive numbers.")
  
  if (fit.intercept) {
    tau_0 <- prior_check(tau_0, 1)
    if (is.null(tau)) stop("tau_0 must be a positive number")
    u_0 <- prior_check(u_0, 1)
    if (is.null(u)) stop("u_0 must either be a positive number.")
    v_0 <- prior_check(v_0, 1)
    if (is.null(v_0)) stop("v_0 must be a positive number.")
  }
  
  # uniform prior on change point locations if not specified
  pi <- prob_check(pi, fit.intercept, L, T)
  if(is.null(pi)) stop("If given, pi must be a T x L matrix with columns that sum to one.")
  omega <- prob_check(omega, fit.intercept, K, T)
  if(is.null(omega)) stop("If given, omega must be a T x K matrix with columns that sum to one.")
  
  # initializing posterior parameters
  
  # mean components
  pi_bar_ij <- matrix(1 / T, nrow = T, ncol = L)
  mu_bar_ij <- matrix(0, nrow = T, ncol = L)
  tau_bar_ij <- matrix(1, nrow = T, ncol = L)
  beta_bar_l <- apply(pi_bar_ij * mu_bar_ij, 2, cumsum)
  var_beta <- rowSums(apply(pi_bar_ij * (mu_bar_ij^2 + 1 / tau_bar_ij), 2, cumsum) - beta_bar_l^2)

  # scale components
  omega_bar_ij <- matrix(1 / T, nrow = T, ncol = K)
  u_bar_ij <- matrix(u, nrow = T, ncol = K, byrow = TRUE) + (T - 1:T + 1) / 2
  v_bar_ij <- matrix(revcumsum(y^2) + var_beta, nrow = T, ncol = K)
  lambda_bar_k <- matrix(0, nrow = T, ncol = K)
  
  for (k in 1:K) {
    lambda_bar_k[,k] <- lambda_bar_fn(u_bar_ij[,k], v_bar_ij[,k], omega_bar_ij[,k])
  }
  
  lambda_bar <- apply(lambda_bar_k, 1, prod)
  
  # store current parameter values
  params <- c(tau_bar_ij, mu_bar_ij, pi_bar_ij, v_bar_ij, omega_bar_ij)
  
  # intercept parameters
  if (fit.intercept) {
    ng_fit <- ng(c(y_0, y), c(1, lambda_bar), tau_0, v_0 + sum(lambda_bar * var_beta) / 2)
    
    mu_bar_0 <- ng_fit$mu
    tau_bar_0 <- ng_fit$tau
    v_bar_0 <- ng_fit$v
    u_bar_0 <- (T + 1) / 2 + u_0
    
    s_bar_0 <- u_bar_0 / v_bar_0 
    lambda_bar <- s_bar_0 * lambda_bar
    
    params <- c(mu_bar_0, v_bar_0, params)
  }

  # initialize residual
  r_bar <- y - ifelse(fit.intercept, mu_bar_0, 0)
  
  n_iter <- 1
  
  while (TRUE) {
    
    new_params <- c()
    
    # updating q(b_0, s_0)
    if (fit.intercept) {
      lambda_bar <- lambda_bar / s_bar_0 # dividing out intercept precision
      
      r_bar <- r_bar + mu_bar_0 # intercept partial residual
      
      ng_fit <- ng(c(y_0, r_bar), c(1, lambda_bar), tau_0, v_0 + sum(lambda_bar * var_beta) / 2) # fit ng model on intercept residual
      
      mu_bar_0 <- ng_fit$mu
      tau_bar_0 <- ng_fit$tau
      v_bar_0 <- ng_fit$v
      
      r_bar <- r_bar - mu_bar_0 # full residual
      s_bar_0 <- u_bar_0 / v_bar_0 # expected posterior intercept precision 
      lambda_bar <- s_bar_0 * lambda_bar # multiplying by intercept precision
      
      new_params <- c(mu_bar_0, v_bar_0)
    }
    
    # updating q(b_i, gamma_i)
    for (l in 1:L) {
      r_bar <- r_bar + beta_bar_l[,l] # partial mean residual 
      
      smcp_fit <- smcp(r_bar,  lambda_bar, tau[l], pi[,l]) # fit single smcp model on partial mean residual
      
      # store posterior parameters
      mu_bar_ij[,l] <- smcp_fit$mu
      tau_bar_ij[,l] <- smcp_fit$tau
      pi_bar_ij[,l] <- smcp_fit$pi
      
      beta_bar_l[,l] <- cumsum(pi_bar_ij[,l] * mu_bar_ij[,l])
      r_bar <- r_bar - beta_bar_l[,l] # subtracting l^{th} mean change back from residual
    }
    
    var_beta <- rowSums(apply(pi_bar_ij * (mu_bar_ij^2 + 1 / tau_bar_ij), 2, cumsum) - beta_bar_l^2)
    
    z2_bar <- lambda_bar * (r_bar^2 + fit.intercept * (v_bar_0 / (tau_bar_0 * (u_bar_0 - 1))) + var_beta) # squared scale residual
    
    # updating q(s_i, alpha_i)

    for (k in 1:K) {
      z2_bar <- z2_bar / lambda_bar_k[,k] # squared partial scale residual
      lambda_bar <- lambda_bar / lambda_bar_k[,k] # dividing out lambda_k
        
      sscp_fit <- sscp(sqrt(z2_bar), u[k], v[k], omega[,k]) # fit single sscp model on partial scale residual

      # store posterior parameters
      v_bar_ij[,k] <- sscp_fit$v
      omega_bar_ij[,k] <- sscp_fit$omega
      
      lambda_bar_k[,k] <- lambda_bar_fn(u_bar_ij[,k], v_bar_ij[,k], omega_bar_ij[,k])
      z2_bar <- z2_bar * lambda_bar_k[,k] # multiplying k^{th} mean change back to residual
      lambda_bar <- lambda_bar * lambda_bar_k[,k] # multiplying back lambda_k
    }
    
    new_params <- c(new_params, c(tau_bar_ij, mu_bar_ij, pi_bar_ij, v_bar_ij, omega_bar_ij))
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 1000 == 0) print(paste0("Iteration ", n_iter,", Error: ", sqrt(sum((params - new_params)^2))))
    n_iter <- n_iter + 1
    params <- new_params
  }
  
  ret <- list(y = y, 
              beta = rowSums(beta_bar_l),
              lambda = lambda_bar,
              pi = pi_bar_ij, mu = mu_bar_ij, sigma2 = 1 / tau_bar_ij, 
              omega = omega_bar_ij, v = v_bar_ij, u = u_bar_ij,
              fit.intercept = fit.intercept)
  
  if (fit.intercept) {
    ret <- c(ret, list(mu_0 = mu_bar_0, u_0 = u_bar_0, v_0 = v_bar_0))
    ret$y <- c(y_0, y)
    ret$beta <- c(mu_bar_0, mu_bar_0 + ret$beta)
    ret$lambda <- c(s_bar_0, ret$lambda)
  }
  
  return(ret)
}

#### MICH II ####

mich_ii <- function(y, L, tol = 1e-5, fit.intercept = TRUE,
                    tau = 0.1, u = 1e-3, v = 1e-3, tau_0 = 0.1, u_0 = 1e-3, v_0 = 1e-3,
                    pi = NULL) {
  
  if (!is.numeric(y) | !is.vector(y)) stop("y must be a numeric vector.")
  
  # model parameters
  if (fit.intercept) {
    y_0 <- y[1]
    y <- y[-1]
  }
  
  T <- length(y)
  
  # checking that priors are proper
  tau <- prior_check(tau, L)
  if (is.null(tau)) stop("tau must either be a positive number or length L vector of positive number.")
  u <- prior_check(u, L)
  if (is.null(u)) stop("u must either be a positive number or length L vector of positive number.")
  v <- prior_check(v, L)
  if (is.null(v)) stop("v must either be a positive number or length L vector of positive numbers.")
  
  if (fit.intercept) {
    tau_0 <- prior_check(tau_0, 1)
    if (is.null(tau)) stop("tau_0 must be a positive number")
    u_0 <- prior_check(u_0, 1)
    if (is.null(u)) stop("u_0 must either be a positive number.")
    v_0 <- prior_check(v_0, 1)
    if (is.null(v_0)) stop("v_0 must be a positive number.")
  }
  
  # uniform prior on change point locations if not specified
  pi <- prob_check(pi, fit.intercept, L, T)
  if(is.null(pi)) stop("If given, pi must be a T x L matrix with columns that sum to one.")
  
  # mean components
  mu_bar_ij <- matrix(0, nrow = T, ncol = L)
  tau_bar_ij <- matrix(1, nrow = T, ncol = L)
  pi_bar_ij <- matrix(1 / T, nrow = T, ncol = L)
  
  # scale components
  u_bar_ij <- matrix(u, nrow = T, ncol = L, byrow = TRUE) + (T - 1:T + 1) / 2
  v_bar_ij <- u_bar_ij
  lambda_bar_l <- matrix(0, nrow = T, ncol = L)
  
  for (l in 1:L) {
    lambda_bar_l[,l] <- lambda_bar_fn(u_bar_ij[,l], v_bar_ij[,l], pi_bar_ij[,l])
  }
  
  lambda_bar <- apply(lambda_bar_l, 1, prod)
  
  beta_lambda <- apply(pi_bar_ij * mu_bar_ij * (u_bar_ij / v_bar_ij), 2, cumsum) / lambda_bar_l
  beta2_lambda <- apply((mu_bar_ij^2 * (u_bar_ij / v_bar_ij) + 1 / tau_bar_ij) * pi_bar_ij, 2, cumsum) / lambda_bar_l
  var_beta <- rowSums(beta2_lambda - beta_lambda^2)
  
  # store current parameter values
  params <- c(tau_bar_ij, mu_bar_ij, pi_bar_ij, v_bar_ij)
  
  # intercept parameters
  if (fit.intercept) {
    delta_0 <- sum(lambda_bar * var_beta) / 2 # variance correction term
    ng_fit <- ng(c(y_0, y), c(1, lambda_bar), tau_0, v_0 + delta_0)
    
    mu_bar_0 <- ng_fit$mu
    tau_bar_0 <- ng_fit$tau
    v_bar_0 <- ng_fit$v
    u_bar_0 <- (T + 1) / 2 + u_0
    
    s_bar_0 <- u_bar_0 / v_bar_0 
    lambda_bar <- s_bar_0 * lambda_bar
    
    params <- c(mu_bar_0, v_bar_0, params)
  }
  
  # initialize residual
  r_tilde <- y - ifelse(fit.intercept, mu_bar_0, 0)
  
  n_iter <- 1
  
  while (TRUE) {
    
    new_params <- c()
    
    # updating q(b_0, s_0)
    if (fit.intercept) {
      lambda_bar <- lambda_bar / s_bar_0 # dividing out intercept precision
      
      r_tilde <- r_tilde + mu_bar_0 # intercept partial residual
      delta_0 <- sum(lambda_bar * var_beta) / 2 # variance correction term
      
      ng_fit <- ng(c(y_0, r_tilde), c(1, lambda_bar), tau_0, v_0 + delta_0) # fit ng model on intercept residual
      
      mu_bar_0 <- ng_fit$mu
      tau_bar_0 <- ng_fit$tau
      v_bar_0 <- ng_fit$v
      
      r_tilde <- r_tilde - mu_bar_0 # full residual
      s_bar_0 <- u_bar_0 / v_bar_0 # expected posterior intercept precision 
      lambda_bar <- s_bar_0 * lambda_bar # multiplying by intercept precision
      
      new_params <- c(mu_bar_0, v_bar_0)
    }
    
    # updating q(b_i, s_i, gamma_i)
    for (l in 1:L) {
      lambda_bar <- lambda_bar / lambda_bar_l[,l] # divide out l^{th} precision component
      r_tilde <- r_tilde + beta_lambda[,l] # partial mean residual 
      var_beta <- var_beta - beta2_lambda[,l] + beta_lambda[,l]^2
      delta_l <- revcumsum(lambda_bar * (var_beta + fit.intercept / (s_bar_0 * tau_bar_0))) / 2 # variance correction term
      
      smscp_fit <- smscp(r_tilde, lambda_bar, tau[l], u[l], v[l] + delta_l, pi[,l]) # fit single smscp model on modified partial residual
      
      # store posterior parameters
      mu_bar_ij[,l] <- smscp_fit$mu
      tau_bar_ij[,l] <- smscp_fit$tau
      v_bar_ij[,l] <- smscp_fit$v
      pi_bar_ij[,l] <- smscp_fit$pi
      
      lambda_bar_l[,l] <- lambda_bar_fn(u_bar_ij[,l], v_bar_ij[,l], pi_bar_ij[,l])
      lambda_bar <- lambda_bar * lambda_bar_l[,l] # multiplying back lambda_l
      beta_lambda[,l] <- cumsum(pi_bar_ij[,l] * mu_bar_ij[,l] * (u_bar_ij[,l] / v_bar_ij[,l])) / lambda_bar_l[,l]
      beta2_lambda[,l] <- cumsum((mu_bar_ij[,l]^2 * (u_bar_ij[,l] / v_bar_ij[,l]) + 1 / tau_bar_ij[,l]) * pi_bar_ij[,l]) / lambda_bar_l[,l]
      var_beta <- var_beta + beta2_lambda[,l] - beta_lambda[,l]^2
      
      r_tilde <- r_tilde - beta_lambda[,l] # subtracting l^{th} mean change back from residual
    }
    
    new_params <- c(new_params, c(tau_bar_ij, mu_bar_ij, pi_bar_ij, v_bar_ij))
    
    # l^2 convergence check
    if (sqrt(sum((params - new_params)^2)) < tol) break
    if (n_iter %% 1000 == 0) print(paste0("Iteration ", n_iter,", Error: ", sqrt(sum((params - new_params)^2))))
    n_iter <- n_iter + 1
    params <- new_params
    
  }
  
  ret <- list(y = y, 
              beta = cumsum(rowSums(mu_bar_ij * pi_bar_ij)),
              lambda = lambda_bar,
              pi = pi_bar_ij, mu = mu_bar_ij, sigma2 = 1 / tau_bar_ij, 
              v = v_bar_ij, u = u_bar_ij,
              fit.intercept = fit.intercept)
  
  if (fit.intercept) {
    ret <- c(ret, list(mu_0 = mu_bar_0, u_0 = u_bar_0, v_0 = v_bar_0))
    ret$y <- c(y_0, y)
    ret$beta <- c(mu_bar_0, mu_bar_0 + ret$beta)
    ret$lambda <- c(s_bar_0, ret$lambda)
  }
  
  return(ret)
}




