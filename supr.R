#### helper functions #### 

revcumsum <- function(x){
  rev(cumsum(rev(x)))
}

lambda_bar_fn <- function(u, v, prob) {
  cumsum((u / v) * prob) + c(revcumsum(prob[-1]), 0)
}

prob_check <- function(probs, fit.intercept, n, T) {
  test <- TRUE
  if (is.null(probs)) {
    probs <- matrix(1/T, ncol = n, nrow = T)
  } else {
    if (fit.intercept & nrow(probs) == (T + 1)) probs <- probs[-1,] / colSums(probs[-1,])  # delete prob that first point is a change if fitting intercept and too many pi
    if (!is.numeric(probs)) test <- FALSE
    if (is.array(probs) & (nrow(probs) != T | ncol(probs) != n | any(round(colSums(probs), 10) != 1))) test <- FALSE
  }
  if (test) return(probs)
  else return(NULL)
}

prior_check <- function(prior, n) {
  test <- TRUE
  if (!is.numeric(prior)) test <- FALSE
  else if (any(prior < 0) | (length(prior) > 1 & length(prior) != n)) test <- FALSE
  if (!test) return(NULL)
  else if (length(prior) == 1) return(rep(prior, n))
  else return (prior)
}

ng <- function(y, lambda, tau, v) {
  tau <- sum(lambda) + tau
  mu <- sum(y * lambda) / tau
  v <- v - (tau * mu^2) / 2 + sum(lambda * y^2) / 2
  return(list(tau = tau, mu = mu, v = v))
} 

smcp <- function(y, pi, tau, lambda) {
  tau <- tau + revcumsum(lambda) 
  mu <- revcumsum(lambda * y) / tau
  log_pi <- log(pi) - 0.5 * log(tau) + (tau * mu^2) / 2
  pi <- prop.table(exp(log_pi - max(log_pi)))
  return(list(tau = tau, mu = mu, pi = pi))
}

sscp <- function(y, omega, u, v) {
  T <- length(y)
  v <- v + revcumsum(y^2) / 2
  u <- u + (T - 1:T + 1) / 2
  log_omega <- log(omega) + lgamma(u) - u * log(v) - cumsum(c(0,y[-T]^2)) / 2 
  omega <- prop.table(exp(log_omega - max(log_omega)))
  return(list(v = v, omega = omega))
}

#### Main #### 

supr <- function(y, L, K, tol = 1e-5, fit.intercept = TRUE,
                 tau = 0.1, u = 1e-3, v = 1e-3, tau_0 = 0.1, u_0 = 1e-3, v_0 = 1e-3,
                 pi = NULL, omega = NULL) {
  
  if (!is.numeric(y)) stop("y must be a numeric vector.")
  
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
    if (is.null(tau)) stop("tau must be a positive number")
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
      
      smcp_fit <- smcp(r_bar, pi[,l], tau[l], lambda_bar) # fit single smcp model on partial mean residual
      
      # store posterior parameters
      mu_bar_ij[,l] <- smcp_fit$mu
      tau_bar_ij[,l] <- smcp_fit$tau
      pi_bar_ij[,l] <- smcp_fit$pi
      
      beta_bar_l[,l] <- cumsum(pi_bar_ij[,l] * mu_bar_ij[,l])
      r_bar <- r_bar - beta_bar_l[,l] # adding l^{th} mean change back to residual
    }
    
    var_beta <- rowSums(apply(pi_bar_ij * (mu_bar_ij^2 + 1 / tau_bar_ij), 2, cumsum) - beta_bar_l^2)
    
    z2_bar <- lambda_bar * (r_bar^2 + fit.intercept * (v_bar_0 / (tau_bar_0 * (u_bar_0 - 1))) + var_beta) # squared scale residual
    
    # updating q(s_i, alpha_i)

    for (k in 1:K) {
      z2_bar <- z2_bar / lambda_bar_k[,k] # squared partial scale residual
      lambda_bar <- lambda_bar / lambda_bar_k[,k] # dividing out lambda_k
        
      sscp_fit <- sscp(sqrt(z2_bar), omega[,k], u[k], v[k]) # fit single sscp model on partial scale residual

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
    if (n_iter %% 1000 == 0) print(paste("Iteration", n_iter))
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


#### Post Processing #### 

plot_supr_fit <- function(supr_fit, level = 0.95) {
  alph <- (1 - level) / 2
  q <- -qnorm(alph)
  plot(supr_fit$y, ylim = c(min(supr_fit$y - q / sqrt(supr_fit$lambda)), max(supr_fit$y + q / sqrt(supr_fit$lambda))),
       ylab = "")
  lines(supr_fit$beta, col = "red", lwd = 3)
  lines(supr_fit$beta - q / sqrt(supr_fit$lambda), lty = 3, col = "blue", lwd = 3)
  lines(supr_fit$beta + q / sqrt(supr_fit$lambda), lty = 3, col = "blue", lwd = 3)
}

cred_set_susie <- function(probs, fit.intercept = TRUE, level = 0.95) {
  L <- ncol(probs)
  T <- nrow(probs) + fit.intercept
  cs <- list()
  for(l in 1:L) {
    set <- order(probs[,l], decreasing = TRUE)[1:which.max(cumsum(probs[order(probs[,l], decreasing = TRUE), l]) > level)]
    set <- set[order(set)]
    if(length(set) == 1) {
      cs <- c(cs, list(set))
      next
    } 
    cmbs <- combn(set, 2)
    i <- cmbs[1,] + fit.intercept
    j <- cmbs[2,] + fit.intercept
    E_i <- 1 - (i - 1) / T  
    E_j <- 1 - (j - 1) / T
    purity <- min(sqrt(E_j) * (1 - E_i) / sqrt(E_i * (1 - E_i) * (1 - E_j)))
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

