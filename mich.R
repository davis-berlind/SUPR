mich <- function(y, J = 0, L = 0, K = 0, fit.intercept = TRUE, fit.scale = TRUE,
                 tol = 1e-5, B_l = 1, B_r = B_l, verbose = FALSE, max_iter = 10000,
                 conv_crit = "ELBO",
                 tau_j = 1e-3, u_j = 1e-3, v_j = u_j, theta = NULL,
                 tau_l = tau_j, pi = NULL,
                 u_k = u_j, v_k = v_j, omega = NULL) {
  
  if (!is.numeric(y) | !is.vector(y)) stop("y must be a numeric vector.")
  
  # extract y in (1-B_l):1
  y_0 <- c()
  if (B_l > 0) {
    y_0 <- y[1:B_l]
    y <- y[-c(1:B_l)]
  }
  
  T <- length(y) - B_r
  tail <- (T+1):(T+B_r)
  
  # checking that priors are proper
  if (J > 0) {
    tau_j <- prior_check(tau_j, J)
    if (is.null(tau_j)) stop("tau_j must either be a positive number or length J vector of positive number.")
    u_j <- prior_check(u_j, J)
    if (is.null(u_j)) stop("u_j must either be a positive number or length J vector of positive number.")
    v_j <- prior_check(v_j, J)
    if (is.null(v_j)) stop("v_j must either be a positive number or length J vector of positive numbers.") 
  }
  if (L > 0) {
    tau_l <- prior_check(tau_l, L)
    if (is.null(tau_l)) stop("tau_l must either be a positive number or length L vector of positive number.")
  }
  if (K > 0) { 
    u_k <- prior_check(u_k, K)
    if (is.null(u_k)) stop("u_k must either be a positive number or length K vector of positive number.")
    v_k <- prior_check(v_k, K)
    if (is.null(v_k)) stop("v_K must either be a positive number or length K vector of positive numbers.") 
  }
  
  # uniform prior on change point locations if not specified
  if(J > 0) {
    theta <- prob_check(theta, J, T)
    if(is.null(theta)) stop("If given, theta must be a T x J matrix with columns that sum to one.")
  }
  if(L > 0) {
    pi <- prob_check(pi, L, T)
    if(is.null(pi)) stop("If given, pi must be a T x L matrix with columns that sum to one.")
  }
  if(K > 0) {
    omega <- prob_check(omega, K, T)
    if(is.null(omega)) stop("If given, omega must be a T x K matrix with columns that sum to one.")
  }
  
  # initialize parameter vector
  if (conv_crit == "l2") params <- c()
  
  # initializing posterior parameters
  if (J > 0) {
    theta_bar_j <- matrix(1 / T, nrow = T, ncol = J)
    b_bar_j <- matrix(0, nrow = T, ncol = J)
    tau_bar_j <- matrix(1, nrow = T, ncol = J)
    mu_lambda_j <- matrix(0, nrow = T + B_r, ncol = J)
    mu2_lambda_j <- matrix(0, nrow = T + B_r, ncol = J)
    u_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
    v_bar_j <- u_bar_j
    lambda_bar_j <- matrix(1, nrow = T + B_r, ncol = J)
    if (conv_crit == "l2") params <- c(params, theta_bar_j, b_bar_j, tau_bar_j, v_bar_j)
  }
  if (L > 0) {
    pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
    b_bar_l <- matrix(0, nrow = T, ncol = L)
    tau_bar_l <- matrix(1, nrow = T, ncol = L)
    mu_bar_l <- matrix(0, nrow = T + B_r, ncol = L)
    mu2_bar_l <- matrix(0, nrow = T + B_r, ncol = L)
    if (conv_crit == "l2") params <- c(params, pi_bar_l, b_bar_l, tau_bar_l)
  }
  if (K > 0) {
    omega_bar_k <- matrix(1 / T, nrow = T, ncol = K)
    u_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
    v_bar_k <- u_bar_k
    lambda_bar_k <- matrix(1, nrow = T + B_r, ncol = K)
    if (conv_crit == "l2") params <- c(params, omega_bar_k, v_bar_k)
  }
  
  # quick fit to initialize mu_0 and lambda_0
  # if (fit.intercept | fit.scale) {
  #   init_fit <- mich(c(y_0,y), 
  #                    ifelse(J > 0, J+1, 0), 
  #                    ifelse(L > 0, L+1, 0), 
  #                    ifelse(K > 0, K+1, 0),
  #                    fit.intercept = FALSE, fit.scale = FALSE,
  #                    B_l = 0, B_r = B_r, max_iter = 5)
  # }
 
  # initialize mu_0
  if (fit.intercept) {
    #mu_0 <- init_fit$mu[1]
    if (B_l > 0) mu_0 <- mean(y_0)
    else mu_0 <- y[1]
    if (conv_crit == "l2") params <- c(params, mu_0)
  } else mu_0 <- 0
  
  # initialize lambda_0. 
  if (fit.scale) {
    #lambda_0 <- init_fit$lambda[1]
    if (B_l == 1) lambda_0 <- 1
    else if (B_l > 1) lambda_0 <- 1 / var(y_0)
    else lambda_0 <- 1 / var(y[1:2])
    if (conv_crit == "l2") params <- c(params, lambda_0)
  } else lambda_0 <- 1

  # initialize residual and scale vector
  r_tilde <- y - mu_0
  lambda_bar <- rep(lambda_0, T + B_r)
  
  # initialize correction term
  delta <- rep(0, T + B_r)
  
  # initialize ELBO
  elbo <- -Inf
  elbo_track <- rep(elbo, max_iter + 1)
  
  # VB algorithm 
  for (iter in 1:max_iter) {
    if (conv_crit == "l2") new_params <- c()
    
    # updating q(b_j, s_j, beta_j)
    if (J > 0) {
      for (j in 1:J) {
        # deleting j^{th} component from residual terms
        r_tilde <- r_tilde + mu_lambda_j[,j] / lambda_bar_j[,j]
        lambda_bar <- lambda_bar / lambda_bar_j[,j]
        delta <- delta - mu2_lambda_j[,j] / lambda_bar_j[,j] + (mu_lambda_j[,j] / lambda_bar_j[,j])^2
        
        # modified priors
        delta_j <- c(0, cumsum(lambda_bar[1:(T-1)] * delta[1:(T-1)]))
        v_tilde_j <- v_j[j] + 0.5 * ((lambda_bar[T] * delta[T] + delta_j[T]) - delta_j)
        if (B_r > 0) v_tilde_j <- v_tilde_j + 0.5 * sum(lambda_bar[tail] * delta[tail])
        log_theta_tilde_j <- log(theta[,j]) - 0.5 * delta_j
        theta_tilde_j <- prop.table(exp(log_theta_tilde_j - max(log_theta_tilde_j)))
        
        # fit single smscp model on modified partial residual
        smscp_fit <- smscp(r_tilde, lambda_bar,  tau_j[j], u_j[j], v_tilde_j, theta_tilde_j, B_r) 
        
        # store posterior parameters
        b_bar_j[,j] <- smscp_fit$mu
        tau_bar_j[,j] <- smscp_fit$tau
        v_bar_j[,j] <- smscp_fit$v
        theta_bar_j[,j] <- smscp_fit$theta
        
        # update j^{th} mean and scale parameters
        mu_lambda_j[1:T,j] <- cumsum(b_bar_j[,j] * (u_bar_j[,j] / v_bar_j[,j]) * theta_bar_j[,j])
        mu2_lambda_j[1:T,j] <- cumsum((b_bar_j[,j]^2 * (u_bar_j[,j] / v_bar_j[,j]) + 1 / tau_bar_j[,j]) * theta_bar_j[,j])
        lambda_bar_j[1:T,j] <- lambda_bar_fn(u_bar_j[,j], v_bar_j[,j], theta_bar_j[,j])
        if (B_r > 0) {
          mu_lambda_j[tail,j] <- mu_lambda_j[T,j]
          mu2_lambda_j[tail,j] <- mu2_lambda_j[T,j]
          lambda_bar_j[tail,j] <- lambda_bar_j[T,j]
        }

        # update residual terms
        r_tilde <- r_tilde - mu_lambda_j[,j] / lambda_bar_j[,j]
        lambda_bar <- lambda_bar * lambda_bar_j[,j]
        delta <- delta + mu2_lambda_j[,j] / lambda_bar_j[,j] - (mu_lambda_j[,j] / lambda_bar_j[,j])^2
      }
      if (conv_crit == "l2") new_params <- c(new_params, theta_bar_j, b_bar_j, tau_bar_j, v_bar_j)
    }

    # updating q(b_l, gamma_l)
    if (L > 0) {
      for (l in 1:L) {
        # deleting l^{th} component from residual terms
        r_tilde <- r_tilde + mu_bar_l[,l]
        delta <- delta - (mu2_bar_l[,l] - mu_bar_l[,l]^2)
        
        # fit single smcp model on modified partial residual
        smcp_fit <- smcp(r_tilde, lambda_bar, tau_l[l], pi[,l], B_r) 
        
        # store posterior parameters
        b_bar_l[,l] <- smcp_fit$mu
        tau_bar_l[,l] <- smcp_fit$tau
        pi_bar_l[,l] <- smcp_fit$pi
        
        # update l^{th} mean parameters
        mu_bar_l[1:T,l] <- cumsum(b_bar_l[,l] * pi_bar_l[,l])
        mu2_bar_l[1:T,l] <- cumsum((b_bar_l[,l]^2 + 1 / tau_bar_l[,l]) * pi_bar_l[,l])
        if (B_r > 0) {
          mu_bar_l[tail,l] <- mu_bar_l[T,l]
          mu2_bar_l[tail,l] <- mu2_bar_l[T,l]
        }
        
        # update residual terms
        r_tilde <- r_tilde - mu_bar_l[,l]
        delta <- delta + (mu2_bar_l[,l] - mu_bar_l[,l]^2)
      }
      if (conv_crit == "l2") new_params <- c(new_params, pi_bar_l, b_bar_l, tau_bar_l)
    }

    # updating q(s_k, alpha_k)
    if (K > 0) {
      for (k in 1:K) {
        # deleting k^{th} component from residual terms
        lambda_bar <- lambda_bar / lambda_bar_k[,k]
        
        # modified priors
        delta_k <- c(0, cumsum(lambda_bar[1:(T-1)] * delta[1:(T-1)]))
        v_tilde_k <- v_k[k] + 0.5 * ((lambda_bar[T] * delta[T] + delta_k[T]) - delta_k)
        if (B_r > 0) v_tilde_k <- v_tilde_k + 0.5 * sum(lambda_bar[tail] * delta[tail])
        log_omega_tilde_k <- log(omega[,k]) - 0.5 * delta_k
        omega_tilde_k <- prop.table(exp(log_omega_tilde_k - max(log_omega_tilde_k)))
        
        # fit single sscp model on modified partial residual
        sscp_fit <- sscp(r_tilde, lambda_bar, u_k[k], v_tilde_k, omega_tilde_k, B_r) 
        
        # store posterior parameters
        v_bar_k[,k] <- sscp_fit$v
        omega_bar_k[,k] <- sscp_fit$omega
        
        # update j^{th} scale parameters
        lambda_bar_k[1:T,k] <- lambda_bar_fn(u_bar_k[,k], v_bar_k[,k], omega_bar_k[,k])
        if (B_r > 0) lambda_bar_k[tail,k] <- lambda_bar_k[T,k]
        
        # update residual terms
        lambda_bar <- lambda_bar * lambda_bar_k[,k] # multiplying back lambda_k
      }
      if (conv_crit == "l2") new_params <- c(new_params, omega_bar_k, v_bar_k)
    }
    
    # updating mu_0 and lambda_0
    r_tilde <- r_tilde + mu_0 
    lambda_bar <- lambda_bar / lambda_0
    
    if (fit.intercept) {
      mu_0 <- (sum(y_0) + sum(lambda_bar * r_tilde)) / (B_l + sum(lambda_bar))
      if (conv_crit == "l2") new_params <- c(new_params, mu_0)
    }
    if (fit.scale){
      lambda_0 <- (T + B_l + B_r) / (sum((y_0 - mu_0)^2) + sum(lambda_bar * (r_tilde - mu_0)^2 + delta))
      if (conv_crit == "l2") new_params <- c(new_params, lambda_0)
    }
    
    r_tilde <- r_tilde - mu_0 
    lambda_bar <- lambda_bar * lambda_0
    
    # calculate ELBO (up to constant)
    elbo <- ((T + B_l + B_r) * log(lambda_0) - sum(lambda_0 * (y_0 - mu_0)^2) - sum(lambda_bar * (r_tilde^2 + delta))) / 2
    if (J > 0) {
      elbo <- elbo + sum((T+B_r):(1+B_r) * rowSums((digamma(u_bar_j) - log(v_bar_j)) * theta_bar_j)) / 2
      
      # need special care for log odds when prob is near zero
      log_odds <- log(theta) - log(theta_bar_j)
      log_odds[round(theta_bar_j, 10) == 0] <- 0
      
      # E[log p - log q]
      # variance component
      log_pq_var <- (u_j - u_bar_j) * digamma(u_bar_j)  + u_j * log(v_bar_j) + lgamma(u_bar_j) + u_bar_j * (1 - v_j / v_bar_j)
      
      # mean component
      log_pq_mean <- tau_j * (b_bar_j^2 * u_bar_j / v_bar_j + 1 / tau_bar_j) + log(tau_bar_j)
      
      log_pq <- log_odds - log_pq_mean / 2 + log_pq_var
      elbo <- elbo + sum(theta_bar_j * log_pq)
    }
    if (L > 0) {
      # need special care for log odds when prob is near zero
      log_odds <- log(pi) - log(pi_bar_l)
      log_odds[round(pi_bar_l, 10) == 0] <- 0

      # E[log p - log q]
      log_pq <- log_odds - (tau_l * (b_bar_l^2 + 1 / tau_bar_l) + log(tau_bar_l)) / 2
      elbo <- elbo + sum(pi_bar_l * log_pq)
    } 
    if (K > 0) {
      elbo <- elbo + sum((T+B_r):(1+B_r) * rowSums((digamma(u_bar_k) - log(v_bar_k) ) * omega_bar_k)) / 2

      # need special care for log odds when prob is near zero
      log_odds <- log(omega) - log(omega_bar_k)
      log_odds[round(omega_bar_k, 10) == 0] <- 0

      # E[log p - log q]
      log_pq <- log_odds + (u_k - u_bar_k) * digamma(u_bar_k)  + u_k * log(v_bar_k) + lgamma(u_bar_k) + u_bar_k * (1 - v_k / v_bar_k)
      elbo <- elbo + sum(omega_bar_k * log_pq)
    }
    
    # l2 convergence check
    if (conv_crit == "l2") {
      error <- sqrt(sum((params - new_params)^2))
      params <- new_params
    } 
    
    # elbo convergence check
    elbo_track[iter+1] <- elbo
    if (conv_crit == "ELBO") {
      error <- elbo - elbo_track[iter]
    }
    
    if (error < tol) break
    if (verbose & iter %% 1000 == 0) print(paste0("Iteration ", iter,", Error: ", error))      
  }
  
  # reassemble y 
  if (B_l > 0) {
    y <- c(y_0, y)
  }
  
  # construct mean signal 
  mu <- rep(mu_0, T)
  if (J > 0) mu <- mu + rowSums(apply(b_bar_j * theta_bar_j, 2, cumsum)) 
  if (L > 0) mu <- mu + rowSums(apply(b_bar_l * pi_bar_l, 2, cumsum)) 
  mu <- c(rep(mu_0, B_l), mu, rep(mu[T], B_r))
  
  # construct scale signal 
  lambda_bar <- c(rep(lambda_0, B_l), lambda_bar)
  
  ret <- list(y = y, mu = mu, lambda = lambda_bar, 
              J = J, K = K, L = L, 
              elbo = elbo_track[1:iter+1], converged = (max_iter > iter))
  
  if (fit.intercept) ret$mu_0 <- mu_0
  if (fit.scale) ret$lambda_0 <- lambda_0
  if (J > 0) {
    ret$mean.scale.model <- list(b = b_bar_j, tau = tau_bar_j, 
                                 u = u_bar_j, v = v_bar_j,
                                 probs = rbind(matrix(0, ncol = J, nrow = B_l), 
                                               theta_bar_j,
                                               matrix(0, ncol = J, nrow = B_r)))
  } 
  if (L > 0) {
    ret$mean.model <- list(b = b_bar_l, tau = tau_bar_l, 
                           probs = rbind(matrix(0, ncol = L, nrow = B_l), 
                                         pi_bar_l,
                                         matrix(0, ncol = L, nrow = B_r)))
  }
  if (K > 0) {
    ret$scale.model <- list(u = u_bar_k, v = v_bar_k, 
                            probs = rbind(matrix(0, ncol = K, nrow = B_l), 
                                          omega_bar_k,
                                          matrix(0, ncol = K, nrow = B_r)))
  }
  return(ret)
}

