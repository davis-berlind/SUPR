mich4 <- function(y, fit_intercept = FALSE, fit_scale = FALSE,
                 J = 0, L = 0, K = 0,
                 J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                 tol = 1e-3, B_l = 0, B_r = 0, verbose = FALSE, max_iter = 1e4,
                 tau_j = 1e-3, u_j = 1e-3, v_j = 1e-3, pi_j = NULL,
                 tau_l = 1e-3, pi_l = NULL,
                 u_k = 1e-3, v_k = 1e-3, pi_k = NULL) {
  
  #### checking parameters ####
  if (!is.numeric(y) | !is.vector(y)) stop("y must be a numeric vector.")
  
  B_l <- integer_check(B_l)
  B_r <- integer_check(B_r)
  if (is.null(B_l) | is.null(B_r)) stop("B_l and B_r must be integers >= 0")
  
  # calculate length of sequence
  T = length(y) - B_r - B_l
  
  # extract y in (1-B_l):1
  y_0 <- numeric(0)
  if (B_l > 0) {
    y_0 <- y[1:B_l]
    y <- y[-c(1:B_l)]
  }
  
  max_iter <- integer_check(max_iter)
  if (is.null(max_iter)) stop("max_iter must be a positive integer or Inf.")
  
  if (!is.logical(fit_intercept) | !is.logical(fit_scale)) { 
    stop("fit_intercept and fit_scale must be either TRUE or FALSE.")
  }
  
  if (!is.logical(J_auto)) stop("J_auto must be either TRUE or FALSE.")
  if (!is.logical(L_auto)) stop("L_auto must be either TRUE or FALSE.")
  if (!is.logical(K_auto)) stop("K_auto must be either TRUE or FALSE.")
  
  if (!is.logical(verbose)) stop("verbose must be either TRUE or FALSE.")
  
  J <- integer_check(J)  
  if (is.null(J)) stop("J must be an integer >= 0.")
  
  L <- integer_check(L)  
  if (is.null(L)) stop("L must be an integer >= 0.")
  
  K <- integer_check(K)  
  if (is.null(K)) stop("K must be an integer >= 0.")
  
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
  
  # checking prior on change-point locations.
  if (J_auto & !is.vector(pi_j) & !is.null(pi_j)) stop("when J_auto = TRUE, pi_j must be a length T vector or NULL.")
  if ((J_auto | J > 0) & is.null(pi_j)) pi_j <- prior_j(T, B_r, 1e5, 1e-3, 1e-3, 1e-3)
  pi_j <- prob_check(pi_j, J, T)
  if(is.null(pi_j)) stop("If given, pi_j must be a length T vector or a T x J matrix with columns that sum to one.")
  
  if (L_auto & !is.vector(pi_l) & !is.null(pi_l)) stop("when L_auto = TRUE, pi_l must be a length T vector or NULL.")
  if ((L_auto | L > 0) & is.null(pi_l)) pi_l <- prior_l(T, B_r, 1e5, 1e-3)
  pi_l <- prob_check(pi_l, L, T)
  if(is.null(pi_l)) stop("If given, pi_l must be a length T vector or a T x L matrix with columns that sum to one.")
  
  if (K_auto & !is.vector(pi_k) & !is.null(pi_k)) stop("when K_auto = TRUE, pi_k must be a length T vector or NULL.")
  if ((K_auto | K > 0) & is.null(pi_k)) pi_k <- prior_k(T, B_r, 1e5, 1e-3, 1e-3)
  pi_k <- prob_check(pi_k, K, T)
  if(is.null(pi_k)) stop("If given, pi_k must be a length T vector or a T x K matrix with columns that sum to one.")
  
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
    counter <- ceiling(log(T))  
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo
    restart <- TRUE
    
    if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
    
    while (TRUE) {
      
      if (J_auto) {
        # increment dimension of parameters
        J <- J + 1 
        tau_j <- c(tau_j, 0.001); u_j <- c(u_j, 0.001); v_j <- c(v_j, 0.001);
        if (J > 1) pi_j <- cbind(pi_j[,1], pi_j)
        pi_bar_j <- cbind(1 / T, pi_bar_j)
        b_bar_j <- cbind(0.0, b_bar_j)
        tau_bar_j <- cbind(1.0, tau_bar_j)
        mu_lambda_j <- cbind(0.0, mu_lambda_j)
        mu2_lambda_j <- cbind(0.0, mu2_lambda_j)
        u_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, u_bar_j) 
        v_bar_j <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, v_bar_j)
        lambda_bar_j <- cbind(1.0, lambda_bar_j)
      }
      
      if (L_auto) {
        # increment dimension of parameters
        L <- L + 1
        tau_l <- c(tau_l, 0.001)
        if (L > 1) pi_l <- cbind(pi_l[,1], pi_l)
        pi_bar_l <- cbind(1 / T, pi_bar_l)
        b_bar_l <- cbind(0.0, b_bar_l)
        tau_bar_l <- cbind(1.0, tau_bar_l)
        mu_bar_l <- cbind(0.0, mu_bar_l)
        mu2_bar_l <- cbind(0.0, mu2_bar_l)
      }
      
      if (K_auto) {
        # increment dimension of parameters
        K <- K + 1
        u_k <- c(u_k, 0.001); v_k <- c(v_k, 0.001);
        if (K > 1) pi_k <- cbind(pi_k[,1], pi_k)
        pi_bar_k <- cbind(1 / T, pi_bar_k)
        u_bar_k <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, u_bar_k) 
        v_bar_k <- cbind(0.001 + (T + B_r - 1:T + 1) / 2, v_bar_k)
        lambda_bar_k <- cbind(1.0, lambda_bar_k)
      }
      
      # fit new model 
      fit_new <- mich_cpp(y_0, y, J, L, K, mu_0, lambda_0,
                          B_l, B_r, fit_intercept, fit_scale, 
                          max_iter, verbose = FALSE, tol, 
                          tau_j, u_j, v_j, log(pi_j), pi_bar_j, b_bar_j, 
                          tau_bar_j, u_bar_j, v_bar_j,
                          mu_lambda_j, mu2_lambda_j, lambda_bar_j, 
                          tau_l, log(pi_l), pi_bar_l, b_bar_l, tau_bar_l, 
                          mu_bar_l, mu2_bar_l,
                          u_k, v_k, log(pi_k), pi_bar_k, u_bar_k, v_bar_k,
                          lambda_bar_k)
      
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new))
      
      if (elbo_new > elbo) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- ceiling(log(T))  
        restart <- TRUE
      } else {
        if (J_auto) {
          if (restart & J > 1) {
            J <- J - 1
            tau_j <- tau_j[-1]; u_j <- u_j[-1]; v_j <- v_j[-1];
            pi_j <- pi_j[,-1, drop = FALSE]
          }
          pi_bar_j <- matrix(1 / T, nrow = T, ncol = J)
          b_bar_j <- matrix(0.0, nrow = T, ncol = J)
          tau_bar_j <- matrix(1.0, nrow = T, ncol = J)
          mu_lambda_j <- matrix(0.0, nrow = T + B_r, ncol = J)
          mu2_lambda_j <- matrix(0.0, nrow = T + B_r, ncol = J)
          u_bar_j <- matrix(u_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
          v_bar_j <- matrix(v_j, nrow = T, ncol = J, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
          lambda_bar_j <- matrix(1.0, nrow = T + B_r, ncol = J)
        }
        
        if (L_auto) {
          if (restart & L > 1) {
            L <- L - 1
            tau_l <- tau_l[-1];
            pi_l <- pi_l[,-1, drop = FALSE]
          }
          pi_bar_l <- matrix(1 / T, nrow = T, ncol = L)
          b_bar_l <- matrix(0.0, nrow = T, ncol = L)
          tau_bar_l <- matrix(1.0, nrow = T, ncol = L)
          mu_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
          mu2_bar_l <- matrix(0.0, nrow = T + B_r, ncol = L)
        }
        
        if (K_auto) {
          if (restart & K > 1) {
            K <- K - 1
            u_k <- u_k[-c(1,2)]; v_k <- v_k[-c(1,2)];
            pi_k <- pi_k[,-c(1,2), drop = FALSE]
          }
          pi_bar_k <- matrix(1 / T, nrow = T, ncol = K)
          u_bar_k <- matrix(u_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
          v_bar_k <-  matrix(v_k, nrow = T, ncol = K, byrow = TRUE) + (T + B_r - 1:T + 1) / 2
          lambda_bar_k <- matrix(1.0, nrow = T + B_r, ncol = K)
        }
        
        restart <- FALSE
        counter <- counter - 1
      } 
      
      if (counter == 0) break
    }
  }
  
  #### auto procedure with multiple components ####
  
  if (sum(J_auto, L_auto, K_auto) > 1) {
    done <- FALSE
    counter <- (J_auto + L_auto + K_auto) * ceiling(log(T)) 
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo))
    
    while (TRUE) {
      elbo_new <- rep(-Inf, 3)
      if (J_auto) {
        J <- J + 1
        if (J > 1) pi_j <- cbind(pi_j[,1], pi_j)
        
        # fit new model 
        fit_new <- mich(c(y_0, y), fit_intercept, fit_scale, 
                        J = J, L = L, K = K,
                        pi_j = pi_j, pi_l = pi_l, pi_k = pi_k,
                        B_l = B_l, B_r = B_r, tol = tol, max_iter = max_iter)
        
        elbo_new[1] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new[1]))
        if (elbo_new[1] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[1]
          counter <- (J_auto + L_auto + K_auto) * ceiling(log(T)) 
        } else {
          counter <- counter - 1
        }
        J <- J - 1
        if (J > 0) pi_j <- pi_j[,-1, drop = FALSE]
      }
      
      if (L_auto) {
        L <- L + 1
        if (L > 1) pi_l <- cbind(pi_l[,1], pi_l)
        
        # fit new model 
        fit_new <- mich(c(y_0, y), fit_intercept, fit_scale, 
                        J = J, L = L, K = K,
                        pi_j = pi_j, pi_l = pi_l, pi_k = pi_k,
                        B_l = B_l, B_r = B_r, tol = tol, max_iter = max_iter)
        
        elbo_new[2] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new[2]))
        if (elbo_new[2] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[2]
          counter <- (J_auto + L_auto + K_auto) * ceiling(log(T)) 
        } else {
          counter <- counter - 1
        }
        L <- L - 1
        if (L > 0) pi_l <- pi_l[,-1, drop = FALSE]
      }
      
      if (K_auto) {
        K <- K + 1
        if (K > 1) pi_k <- cbind(pi_k[,1], pi_k)
        
        # fit new model 
        fit_new <- mich(c(y_0, y), fit_intercept, fit_scale, 
                        J = J, L = L, K = K,
                        pi_j = pi_j, pi_l = pi_l, pi_k = pi_k,
                        B_l = B_l, B_r = B_r, tol = tol, max_iter = max_iter)
        
        elbo_new[3] <- fit_new$elbo[length(fit_new$elbo)]
        
        if (verbose) print(paste0("(J = ", J, ", L = ", L, ", K = ", K, "): ELBO = ", elbo_new[3]))
        if (elbo_new[3] > elbo) {
          fit <- fit_new 
          elbo <- elbo_new[3]
          counter <- (J_auto + L_auto + K_auto) * ceiling(log(T)) 
        } else {
          counter <- counter - 1
        }
        K <- K - 1
        if (K > 0) pi_k <- pi_k[,-1, drop = FALSE]
      }
      if (counter <= 0) {
        done <- TRUE
        break
      }
      
      if (which.max(elbo_new) == 1) {
        J <- J + 1 
        if (J > 1) pi_j <- cbind(pi_j[,1], pi_j)
      } else if (which.max(elbo_new) == 2) {
        L <- L + 1 
        if (L > 1) pi_l <- cbind(pi_l[,1], pi_l)
      } else {
        K <- K + 1
        if (K > 1) pi_k <- cbind(pi_k[,1], pi_k)
      }
      if (done) break
    }
  }
  
  #### return model ####
  # reassemble y
  fit$y <- c(y_0, y)
  
  # add zero probabilities in buffer 
  if (fit$J > 0) fit$J_model$pi <- rbind(matrix(0, nrow = B_l, ncol = fit$J), 
                                         fit$J_model$pi,
                                         matrix(0, nrow = B_r, ncol = fit$J))
  if (fit$L > 0) fit$L_model$pi <- rbind(matrix(0, nrow = B_l, ncol = fit$L), 
                                         fit$L_model$pi,
                                         matrix(0, nrow = B_r, ncol = fit$L))
  if (fit$K > 0) fit$K_model$pi <- rbind(matrix(0, nrow = B_l, ncol = fit$K), 
                                         fit$K_model$pi,
                                         matrix(0, nrow = B_r, ncol = fit$K))
  return(fit)
}