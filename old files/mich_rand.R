mich_rand <- function(y, fit_intercept = TRUE, fit_scale = TRUE,
                 #merge = FALSE, merge_tol = 1 / T,
                 J = 0, L = 0, K = 0,
                 J_auto = FALSE, L_auto = FALSE, K_auto = FALSE,
                 tol = 1e-3, B_l = 0, B_r = 0, verbose = FALSE, max_iter = 1e4,
                 tau_j = 1e-3, u_j = 1e-3, v_j = 1e-3, pi_j = "weighted",
                 tau_l = 1e-3, pi_l = "weighted",
                 u_k = 1e-3, v_k = 1e-3, pi_k = "weighted") {
  
  
  # checking that y is a numeric vector/matrix
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.numeric(y)) stop("y must contain only numeric data.")
  if (!is.vector(y) & !is.matrix(y)) stop("y must be a numeric vector or matrix.")
  
  # calculate dimensions of y and handle NAs
  if (all(is.na(y))) stop("y does not contain any data.")
  if (is.matrix(y)) {
    if (any(is.na(y))) {
      warning("NAs detected in y. Removing missing rows and columns and replacing missing values with preceding observation.") 
      y <- y[!apply(y, 1, function(x) all(is.na(x))), !apply(y, 2, function(x) all(is.na(x))), drop == FALSE]
      for (t in 1:nrow(y)) {
        for (i in 1:ncol(y)) {
          if (is.na(y[t, i])){
            if (t == 1) y[t, i] <- y[which.min(is.na(y[,i])),i]
            else y[t, i] <- y[t-1, i]
          }
        }
      }
    }
    d <- ncol(y)
    T <- nrow(y) - B_r - B_l
    if (d == 1) y <- as.vector(y)
  } else {
    if (any(is.na(y))) {
      warning("NAs detected in y. Removing missing values.")
      y <- y[!is.na(y)]
    }
    d <- 1
    T <- length(y) - B_r - B_l
  }
  
  if (T < 2) stop("Number of rows/elements of y smaller than B_l + B_r + 1.")
  
  #### checking other model parameters ####
  
  B_l <- integer_check(B_l)
  B_r <- integer_check(B_r)
  max_iter <- integer_check(max_iter)
  tol <- scalar_check(tol)
  fit_intercept <- logical_check(fit_intercept)
  fit_scale <- logical_check(fit_scale)
  verbose <- logical_check(verbose)
  
  J_auto <- logical_check(J_auto)
  L_auto <- logical_check(L_auto)
  K_auto <- logical_check(K_auto)
  
  J <- integer_check(J)  
  L <- integer_check(L)  
  K <- integer_check(K) 
  
  if (J > 0 & J_auto) warning(paste0("J_auto = TRUE and J > 0. Beginning search from J = ", J))
  if (L > 0 & L_auto) warning(paste0("L_auto = TRUE and L > 0. Beginning search from L = ", L))
  if (K > 0 & K_auto) warning(paste0("K_auto = TRUE and K > 0. Beginning search from K = ", K))

  ##### checking that priors are proper ####
  
  # J components
  tau_j <- scalar_check(tau_j)
  u_j <- scalar_check(u_j)
  v_j <- scalar_check(v_j)
  
  pi_j <- prob_check(pi_j, J, T)
  if (J_auto & length(pi_j) > 1) stop("When J_auto = TRUE, pi_j must one of 'uniform' or 'weighted'.")
  if (J == 0 & !J_auto) {
    pi_j <- rep(1, T) 
  } else if (pi_j == "weighted" & (J_auto | J > 0)) {
    pi_j <- meanvar_prior(T, B_r)
  } else if (pi_j == "uniform" & (J_auto | J > 0)) {
    pi_j <- rep(1/T, T) 
  }
  pi_j <- sapply(1:max(1, J), function(i) pi_j)
  
  # L components
  tau_l <- scalar_check(tau_l)
  
  pi_l <- prob_check(pi_l, L, T)
  if (L_auto & length(pi_l) > 1) stop("When L_auto = TRUE, pi_l must one of 'uniform' or 'weighted'.")
  if (L == 0 & !L_auto) {
    pi_l <- rep(1, T) 
  } else if (pi_l == "weighted" & (L_auto | L > 0)) {
    pi_l <- mean_prior(T, B_r, d)
  } else if (pi_l == "uniform" & (L_auto | L > 0)) {
    pi_l <- rep(1/T, T) 
  }
  pi_l <- sapply(1:max(1, L), function(i) pi_l)

  # K components
  u_k <- scalar_check(u_k)
  v_k <- scalar_check(v_k)
  
  pi_k <- prob_check(pi_k, K, T)
  if (K_auto & length(pi_k) > 1) stop("When K_auto = TRUE, pi_k must one of 'uniform' or 'weighted'.")
  if (K == 0 & !K_auto) {
    pi_k <- rep(1, T) 
  } else if (pi_k == "weighted" & (K_auto | K > 0)) {
    pi_k <- meanvar_prior(T, B_r)
  } else if (pi_k == "uniform" & (K_auto | K > 0)) {
    pi_k <- rep(1/T, T) 
  }
  pi_k <- sapply(1:max(1, K), function(i) pi_k)

  # call multivariate or univariate MICH
  if (is.matrix(y)) {
    if (J > 0 | K > 0) warning("MICH currently only suports mean changepoint detection for multivariate y. Setting J = K = 0.")
    if (T + B_l + B_r < d & fit_scale) {
      warning("y has more columns than rows. MICH does not currently support high-dimensional variance estimation. Assuming Var(y) = I.")
      fit_scale <- FALSE
    }
    
    mich_matrix(y, fit_intercept, fit_scale,
                L, L_auto, tol, B_l, B_r, verbose, max_iter,
                tau_l, pi_l)
    
  } else {
    mich_vector_rand(y, fit_intercept, fit_scale,
                J, L, K, J_auto, L_auto, K_auto,
                tol, B_l, B_r, verbose, max_iter,
                tau_j, u_j, v_j, pi_j,
                tau_l, pi_l,
                u_k, v_k, pi_k)
  }
}

