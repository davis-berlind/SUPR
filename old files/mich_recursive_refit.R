mich <- function(y, 
                 J = 0, J_auto = FALSE, 
                 L = 0, L_auto = FALSE, 
                 K = 0, K_auto = FALSE,
                 fit_intercept = FALSE, fit_scale = FALSE,
                 tol = 1e-5, B_l = 1, B_r = 1, max_iter = 10000, verbose = FALSE,
                 tau_j = 1e-3, u_j = 1e-3, v_j = 1e-3, pi_j = NULL,
                 tau_l = 1e-3, pi_l = NULL,
                 u_k = 1e-3, v_k = 1e-3, pi_k = NULL) {
  
  #### checking parameters ####
  if (!is.numeric(y) | !is.vector(y)) stop("y must be a numeric vector.")
  
  B_l <- component_check(B_l)
  B_r <- component_check(B_r)
  if (is.null(B_l) | is.null(B_r)) stop("B_l and B_r must be integers >= 0")
  
  # calculate length of sequence
  T = length(y) - B_r - B_l
  
  # extract y in (1-B_l):1
  y_0 <- c()
  if (B_l > 0) {
    y_0 <- y[1:B_l]
    y <- y[-c(1:B_l)]
  }
  
  max_iter <- component_check(max_iter)
  if (is.null(max_iter)) stop("max_iter must be a positive integer or Inf.")
  
  if (!is.logical(fit_intercept) | !is.logical(fit_scale)) { 
    stop("fit_intercept and fit_scale must be either TRUE or FALSE.")
  }
  
  if (!is.logical(verbose)) stop("verbose must be either TRUE or FALSE.")
  
  J <- component_check(J)  
  if (is.null(J)) stop("J must be an integer >= 0 or set to 'auto'.")
  # J_auto <- (J == "auto")
  # if (J_auto) J = 0 
  
  L <- component_check(L)  
  if (is.null(L)) stop("L must be an integer >= 0 or set to 'auto'.")
  L_auto <- (L == "auto")
  if (L_auto) L = 0 
  
  K <- component_check(K)  
  if (is.null(K)) stop("K must be an integer >= 0 or set to 'auto'.")
  K_auto <- (K == "auto")
  if (K_auto) K = 0 
  
  ##### checking that priors are proper ####
  # J components
  tau_j <- prior_check(tau_j, J)
  if (is.null(tau_j)) stop("tau_j must either be a positive number or length J vector of positive number.")
  u_j <- prior_check(u_j, J)
  if (is.null(u_j)) stop("u_j must either be a positive number or length J vector of positive number.")
  v_j <- prior_check(v_j, J)
  if (is.null(v_j)) stop("v_j must either be a positive number or length J vector of positive numbers.") 
  
  # L components
  tau_l <- prior_check(tau_l, L)
  if (is.null(tau_l)) stop("tau_l must either be a positive number or length L vector of positive number.")
  
  # K components
  u_k <- prior_check(u_k, K)
  if (is.null(u_k)) stop("u_k must either be a positive number or length K vector of positive number.")
  v_k <- prior_check(v_k, K)
  if (is.null(v_k)) stop("v_K must either be a positive number or length K vector of positive numbers.") 
  
  # uniform prior on change point locations if not specified
  pi_j <- prob_check(pi_j, J, T)
  if(is.null(pi_j)) stop("If given, pi_j must be a T x J matrix with columns that sum to one.")
  pi_l <- prob_check(pi_l, L, T)
  if(is.null(pi_l)) stop("If given, pi_l must be a T x L matrix with columns that sum to one.")
  pi_k <- prob_check(pi_k, K, T)
  if(is.null(pi_k)) stop("If given, pi_k must be a T x K matrix with columns that sum to one.")
  
  #### initializing posterior parameters ####
  # J components
  pi_bar_j <- matrix(1 / T, nrow = T, ncol = J)
  b_bar_j <- matrix(0.0, nrow = T, ncol = J)
  tau_bar_j <- matrix(1.0, nrow = T, ncol = J)
  mu_lambda_j <- matrix(0.0, nrow = T + B_r, ncol = J)
  mu2_lambda_j <- matrix(0.0, nrow = T + B_r, ncol = J)
  u_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  v_bar_j <- matrix(v_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  lambda_bar_j <- matrix(1.0, nrow = T + B_r, ncol = J)
  
  # L components
  pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
  b_bar_l <- matrix(0.0, nrow = T, ncol = L)
  tau_bar_l <- matrix(1.0, nrow = T, ncol = L)
  mu_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
  mu2_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
  
  # K components
  pi_bar_k <- matrix(1 / T, nrow = T, ncol = K)
  u_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  v_bar_k <-  matrix(v_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
  lambda_bar_k <- matrix(1.0, nrow = T + B_r, ncol = K)
  
  #### initialize mu_0 ####
  if (fit_intercept) {
    if (B_l >= 5) mu_0 <- mean(y_0)
    else mu_0 <- mean(c(y_0, y[1:(5-B_l)]))
  } else mu_0 <- 0.0
  
  #### initialize lambda_0 ####
  if (fit_scale) {
    if (B_l >= 5) lambda_0 <- 1 / var(y_0)
    else lambda_0 <- 1 / var(c(y_0, y[1:(5-B_l)]))
  } else lambda_0 <- 1.0
  
  #### fit model ####
  fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                  B_l, B_r, fit_intercept, fit_scale, 
                  max_iter, (verbose & !any(J_auto, K_auto, L_auto)), tol, 
                  tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                  mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                  tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                  mu_bar_l, mu2_bar_l,
                  u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                  lambda_bar_k)
  
  #### auto procedure with single component ####
  if (sum(J_auto, L_auto, K_auto) == 1) {
    while (TRUE) {
      elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
      if (J_auto) {
        if(verbose) print(paste0("J = ", J, ": ELBO = ", elbo))
        
        # increment dimension of parameters
        J <- J + 1 
        tau_j <- c(tau_j, 0.001); u_j <- c(u_j, 0.001); v_j <- c(v_j, 0.001);
        pi_j <- cbind(1 / T, pi_j)
        pi_bar_j <- cbind(1 / T, pi_bar_j)
        b_bar_j <- cbind(0.0, b_bar_j)
        tau_bar_j <- cbind(1.0, tau_bar_j)
        mu_lambda_j <- cbind(0.0, mu_lambda_j)
        mu2_lambda_j <- cbind(0.0, mu2_lambda_j)
        u_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, u_bar_j) 
        v_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, v_bar_j)
        lambda_bar_j <- cbind(1.0, lambda_bar_j)
        
        # fit new model 
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
        
        # if elbo decreases restart fitting
        if (fit$elbo[length(fit$elbo)] < elbo) {
          new_fit <- mich(y, J = J-1, J_auto = FALSE, 
                          L, L_auto = FALSE, 
                          K, K_auto = FALSE, 
                          fit_intercept, fit_scale,
                          tol, B_l, B_r, max_iter)
          if (new_fit$elbo[length(new_fit$elbo)] < elbo) break
          else {
            return(mich(y, J = J-1, J_auto = TRUE, 
                        L, L_auto = FALSE, 
                        K, K_auto = FALSE, 
                        fit_intercept, fit_scale,
                        tol, B_l, B_r, max_iter, verbose))
          }
        }
      }
      
      if (L_auto) {
        if(verbose) print(paste0("L = ", L, ": ELBO = ", elbo))
        
        # increment dimension of parameters
        L <- L + 1
        tau_l <- c(tau_l, 0.001)
        pi_l <- cbind(1 / T, pi_l)
        pi_bar_l <- cbind(1 / T, pi_bar_l)
        b_bar_l <- cbind(0.0, b_bar_l)
        tau_bar_l <- cbind(1.0, tau_bar_l)
        mu_bar_l <- cbind(0.0, mu_bar_l)
        mu2_bar_l <- cbind(0.0, mu2_bar_l)
        
        # fit new model 
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
        
        # if elbo decreases restart fitting
        if (fit$elbo[length(fit$elbo)] < elbo) {
          L <- L - 1
          tau_l <- tau_l[-1]
          pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
          b_bar_l <- matrix(0.0, nrow = T, ncol = L)
          tau_bar_l <- matrix(1.0, nrow = T, ncol = L)
          mu_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
          mu2_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
          
          fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                          B_l, B_r, fit_intercept, fit_scale, 
                          max_iter, verbose = FALSE, tol, 
                          tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                          mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                          tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                          mu_bar_l, mu2_bar_l,
                          u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                          lambda_bar_k)
        }
      }
      
      if (K_auto) {
        if(verbose) print(paste0("K = ", K, ": ELBO = ", elbo))
        
        # increment dimension of parameters
        K <- K + 1
        u_k <- c(u_k, 0.001); v_k <- c(v_k, 0.001);
        pi_k <- cbind(1 / T, pi_k)
        pi_bar_k <- matrix(1 / T, nrow = T, ncol = K)
        u_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
        v_bar_k <-  matrix(v_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
        lambda_bar_k <- matrix(1.0, nrow = T + B_r, ncol = K)
        
        # fit new model
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
        
        # if elbo decreases restart fitting
        if (fit$elbo[length(fit$elbo)] < elbo) {
          K <- K - 1
          u_k <- u_k[-1]; v_k <- v_k[-1];
          pi_k <- pi_k[,-1, drop = FALSE]
          pi_bar_k <- pi_bar_k[,-1, drop = FALSE]
          u_bar_k <- u_bar_k[,-1, drop = FALSE]
          v_bar_k <- v_bar_k[,-1, drop = FALSE]
          lambda_bar_k <- lambda_bar_k[,-1, drop = FALSE]
          
          fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                          B_l, B_r, fit_intercept, fit_scale, 
                          max_iter, verbose = FALSE, tol, 
                          tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                          mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                          tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                          mu_bar_l, mu2_bar_l,
                          u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                          lambda_bar_k)
        }
      }
      
      # return model if elbo stops increasing
      if (fit$elbo[length(fit$elbo)] < elbo) break
    }
  }
  
  #### auto procedure with single component ####
  
  if (sum(J_auto, L_auto, K_auto) > 1) {
    while (TRUE) {
      elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
      if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo))
      elbo_new <- rep(-Inf, 3)
      
      if (J_auto) {
        # increment dimension of parameters
        J <- J + 1 
        tau_j <- c(tau_j, 0.001); u_j <- c(u_j, 0.001); v_j <- c(v_j, 0.001);
        pi_j <- cbind(1 / T, pi_j)
        pi_bar_j <- cbind(1 / T, pi_bar_j)
        b_bar_j <- cbind(0.0, b_bar_j)
        tau_bar_j <- cbind(1.0, tau_bar_j)
        mu_lambda_j <- cbind(0.0, mu_lambda_j)
        mu2_lambda_j <- cbind(0.0, mu2_lambda_j)
        u_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, u_bar_j) 
        v_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, v_bar_j)
        lambda_bar_j <- cbind(1.0, lambda_bar_j)
        
        # fit new model 
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
        
        # extract elbo
        elbo_new[1] = fit$elbo[length(fit$elbo)]
        
        # decrement dimension of parameters
        J <- J - 1
        tau_j <- tau_j[-1]; u_j <- u_j[-1]; v_j <- v_j[-1];
        pi_j <- pi_j[,-1, drop = FALSE]
        pi_bar_j <- pi_bar_j[,-1, drop = FALSE]
        b_bar_j <- b_bar_j[,-1, drop = FALSE]
        tau_bar_j <- tau_bar_j[,-1, drop = FALSE]
        mu_lambda_j <- mu_lambda_j[,-1, drop = FALSE]
        mu2_lambda_j <- mu2_lambda_j[,-1, drop = FALSE]
        u_bar_j <- u_bar_j[,-1, drop = FALSE]
        v_bar_j <- v_bar_j[,-1, drop = FALSE]
        lambda_bar_j <- lambda_bar_j[,-1, drop = FALSE]
        
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
      }
      
      if (L_auto) {
        # increment dimension of parameters
        L <- L + 1
        tau_l <- c(tau_l, 0.001)
        pi_l <- cbind(1 / T, pi_l)
        pi_bar_l <- cbind(1 / T, pi_bar_l)
        b_bar_l <- cbind(0.0, b_bar_l)
        tau_bar_l <- cbind(1.0, tau_bar_l)
        mu_bar_l <- cbind(0.0, mu_bar_l)
        mu2_bar_l <- cbind(0.0, mu2_bar_l)
        
        # fit new model 
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
        
        # extract elbo
        elbo_new[2] = fit$elbo[length(fit$elbo)]
        
        # decrement dimension of parameters
        L <- L - 1
        tau_l <- tau_l[-1]
        pi_l <- pi_l[,-1, drop = FALSE]
        pi_bar_l <- pi_bar_l[,-1, drop = FALSE]
        b_bar_l <- b_bar_l[,-1, drop = FALSE]
        tau_bar_l <- tau_bar_l[,-1, drop = FALSE]
        mu_bar_l <- mu_bar_l[,-1, drop = FALSE]
        mu2_bar_l <- mu2_bar_l[,-1, drop = FALSE]
        
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
      }
      
      if (K_auto) {
        # increment dimension of parameters
        K <- K + 1
        u_k <- c(u_k, 0.001); v_k <- c(v_k, 0.001);
        pi_k <- cbind(1 / T, pi_k)
        pi_bar_k <- cbind(1 / T, pi_bar_k)
        u_bar_k <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, u_bar_k) 
        v_bar_k <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, v_bar_k)
        lambda_bar_k <- cbind(1.0, lambda_bar_k)
        
        # fit new model
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
        
        # extract elbo
        elbo_new[3] = fit$elbo[length(fit$elbo)]
        
        # decrement dimension of parameters
        K <- K - 1
        u_k <- u_k[-1]; v_k <- v_k[-1];
        pi_k <- pi_k[,-1, drop = FALSE]
        pi_bar_k <- pi_bar_k[,-1, drop = FALSE]
        u_bar_k <- u_bar_k[,-1, drop = FALSE]
        v_bar_k <- v_bar_k[,-1, drop = FALSE]
        lambda_bar_k <- lambda_bar_k[,-1, drop = FALSE]
        
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
      }
      
      # if elbo decreases restart fitting
      if (max(elbo_new) < elbo) {
        # J components
        pi_bar_j <- matrix(1 / T, nrow = T, ncol = J)
        b_bar_j <- matrix(0.0, nrow = T, ncol = J)
        tau_bar_j <- matrix(1.0, nrow = T, ncol = J)
        mu_lambda_j <- matrix(0.0, nrow = T + B_r, ncol = J)
        mu2_lambda_j <- matrix(0.0, nrow = T + B_r, ncol = J)
        u_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
        v_bar_j <- matrix(v_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
        lambda_bar_j <- matrix(1.0, nrow = T + B_r, ncol = J)
        
        # L components
        pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
        b_bar_l <- matrix(0.0, nrow = T, ncol = L)
        tau_bar_l <- matrix(1.0, nrow = T, ncol = L)
        mu_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
        mu2_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
        
        # K components
        pi_bar_k <- matrix(1 / T, nrow = T, ncol = K)
        u_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
        v_bar_k <-  matrix(v_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
        lambda_bar_k <- matrix(1.0, nrow = T + B_r, ncol = K)
        
        # fit new model 
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
        
        if (fit$elbo[length(fit$elbo)] < elbo) break
      } else if (which.max(elbo_new) == 1) {
        # increment dimension of parameters
        J <- J + 1 
        tau_j <- c(tau_j, 0.001); u_j <- c(u_j, 0.001); v_j <- c(v_j, 0.001);
        pi_j <- cbind(1 / T, pi_j)
        pi_bar_j <- cbind(1 / T, pi_bar_j)
        b_bar_j <- cbind(0.0, b_bar_j)
        tau_bar_j <- cbind(1.0, tau_bar_j)
        mu_lambda_j <- cbind(0.0, mu_lambda_j)
        mu2_lambda_j <- cbind(0.0, mu2_lambda_j)
        u_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, u_bar_j) 
        v_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, v_bar_j)
        lambda_bar_j <- cbind(1.0, lambda_bar_j)
        
        # fit new model 
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
      } else if (which.max(elbo_new) == 2) {
        # increment dimension of parameters
        L <- L + 1
        tau_l <- c(tau_l, 0.001)
        pi_l <- cbind(1 / T, pi_l)
        pi_bar_l <- cbind(1 / T, pi_bar_l)
        b_bar_l <- cbind(0.0, b_bar_l)
        tau_bar_l <- cbind(1.0, tau_bar_l)
        mu_bar_l <- cbind(0.0, mu_bar_l)
        mu2_bar_l <- cbind(0.0, mu2_bar_l)
        
        # fit new model 
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
      } else {
        # increment dimension of parameters
        K <- K + 1
        u_k <- c(u_k, 0.001); v_k <- c(v_k, 0.001);
        pi_k <- cbind(1 / T, pi_k)
        pi_bar_k <- cbind(1 / T, pi_bar_k)
        u_bar_k <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, u_bar_k) 
        v_bar_k <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, v_bar_k)
        lambda_bar_k <- cbind(1.0, lambda_bar_k)
        
        # fit new model
        fit <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                        B_l, B_r, fit_intercept, fit_scale, 
                        max_iter, verbose = FALSE, tol, 
                        tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, tau_bar_j, u_bar_j, v_bar_j,
                        mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                        tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                        mu_bar_l, mu2_bar_l,
                        u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                        lambda_bar_k)
      }
      # extract elbo
      elbo <- fit$elbo[length(fit$elbo)]
    }
  }
  
  #### return model ####
  # reassemble y
  fit$y <- c(y_0, y)
  
  # add zero probabilities in buffer 
  if (J > 0) fit$J_model$pi <- rbind(matrix(0, nrow = B_l, ncol = J), 
                                     fit$J_model$pi,
                                     matrix(0, nrow = B_r, ncol = J))
  if (L > 0) fit$L_model$pi <- rbind(matrix(0, nrow = B_l, ncol = L), 
                                     fit$L_model$pi,
                                     matrix(0, nrow = B_r, ncol = L))
  if (K > 0) fit$L_model$pi <- rbind(matrix(0, nrow = B_l, ncol = K), 
                                     fit$K_model$pi,
                                     matrix(0, nrow = B_r, ncol = K))
  return(fit)
}

