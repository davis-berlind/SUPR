mich <- function(y, fit_intercept = TRUE, fit_scale = TRUE,
                 J = 0, L = 0, K = 0,
                 J_auto = FALSE, L_auto = FALSE, K_auto = FALSE, 
                 J_max = Inf, L_max = Inf, K_max = Inf, 
                 tol = 1e-5, merge_prob = NULL, merge_level = 0.95, 
                 max_iter = 1e4, verbose = FALSE, reverse = FALSE, 
                 restart = TRUE, n_restart = 0, increment = 1,
                 omega_j = 1e-3, u_j = 1e-3, v_j = 1e-3, pi_j = "weighted",
                 omega_l = 1e-3, pi_l = "weighted",
                 u_k = 1e-3, v_k = 1e-3, pi_k = "weighted") {
  
  # checking that y is a numeric vector/matrix
  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.numeric(y)) stop("y must contain only numeric data.")
  if (!is.vector(y) & !is.matrix(y)) stop("y must be a numeric vector or matrix.")
  
  # calculate dimensions of y and handle NAs
  if (all(is.na(y)) | length(y) == 0) stop("y does not contain any data.")
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
    T <- nrow(y)
    if (d == 1) y <- as.vector(y)
  } else {
    if (any(is.na(y))) {
      warning("NAs detected in y. Removing missing values.")
      y <- y[!is.na(y)]
    }
    d <- 1
    T <- length(y)
  }
  
  if (T < 2) stop("y must have at least 2 observations.")
  
  #### merge and detect defaults ####
  delta <- 0.5
  merge_level <- scalar_check(merge_level)
  if (merge_level < 0 | merge_level > 1) stop("merge_level must be in [0,1].")
  detect <- ceiling(log(T)^(1 + delta))
  if (is.null(merge_prob)) merge_prob <- detect / T^2
  merge_prob <- scalar_check(merge_prob)
  n_search <- max(4, ceiling(log(T) / ((1 + restart) * increment)))

  #### checking other model parameters ####
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
  
  J_max <- integer_check(J_max)  
  L_max <- integer_check(L_max)  
  K_max <- integer_check(K_max) 
  
  if (J_max < J) "J_max must be >= J."
  if (L_max < L) "L_max must be >= L."
  if (K_max < K) "K_max must be >= K."
  
  if (!J_auto) J_max <- J
  if (!L_auto) L_max <- L
  if (!K_auto) K_max <- K
  
  if (J > 0 & J_auto) warning(paste0("J_auto = TRUE and J > 0. Beginning search from J = ", J))
  if (L > 0 & L_auto) warning(paste0("L_auto = TRUE and L > 0. Beginning search from L = ", L))
  if (K > 0 & K_auto) warning(paste0("K_auto = TRUE and K > 0. Beginning search from K = ", K))

  ##### checking that priors are proper ####
  
  # J components
  omega_j <- scalar_check(omega_j)
  u_j <- scalar_check(u_j)
  v_j <- scalar_check(v_j)
  
  pi_j <- prob_check(pi_j, J, T)
  pi_j_weighted <- (pi_j == "weighted")
  if (J_auto & length(pi_j) > 1) stop("When J_auto = TRUE, pi_j must one of 'uniform' or 'weighted'.")
  if (is.character(pi_j)) {
    if (pi_j == "weighted") log_pi_j <- log_meanvar_prior(T)
    else if (pi_j == "uniform") log_pi_j <- rep(0, T) 
    log_pi_j <- sapply(1:max(1, J), function(i) log_pi_j)
  } else log_pi_j <- log(pi_j)
  
  # L components
  omega_l <- scalar_check(omega_l)
  
  pi_l <- prob_check(pi_l, L, T)
  pi_l_weighted <- (pi_l == "weighted")
  if (L_auto & length(pi_l) > 1) stop("When L_auto = TRUE, pi_l must one of 'uniform' or 'weighted'.")
  if (is.character(pi_l)) {
    if (pi_l == "weighted") log_pi_l <- log_mean_prior(T, floor(1.5 * d))
    else if (pi_l == "uniform") log_pi_l <- rep(0, T) 
    log_pi_l <- sapply(1:max(1, L), function(i) log_pi_l)
  } else log_pi_l <- log(pi_l)

  # K components
  u_k <- scalar_check(u_k)
  v_k <- scalar_check(v_k)
  
  pi_k <- prob_check(pi_k, K, T)
  pi_k_weighted <- (pi_k == "weighted")
  if (K_auto & length(pi_k) > 1) stop("When K_auto = TRUE, pi_k must one of 'uniform' or 'weighted'.")
  if (is.character(pi_k)) {
    if (pi_k == "weighted") log_pi_k <- log_var_prior(T)
    else if (pi_k == "uniform") log_pi_k <- rep(0, T) 
    log_pi_k <- sapply(1:max(1, K), function(i) log_pi_k)
  } else log_pi_k <- log(pi_k)

  # call multivariate or univariate MICH
  if (is.matrix(y)) {
    if (J > 0 | K > 0) stop("MICH currently only suports mean changepoint detection for multivariate y. Try setting J = K = 0.")
    if (T < d & fit_scale) {
      warning("y has more columns than rows. MICH does not currently support high-dimensional variance estimation. Assuming Var(y) = I.")
      fit_scale <- FALSE
    }

    if (reverse) {
      if (!pi_l_weighted) pi_l <- pi_l[T:1,,drop = FALSE]
      
      mich_matrix(y[T:1,], fit_intercept = TRUE, fit_scale,
                  L, L_auto, L_max, pi_l_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob, 
                  restart, n_restart, n_search, increment,
                  omega_l, log_pi_l)
    } else {
      mich_matrix(y, fit_intercept, fit_scale,
                  L, L_auto, L_max, pi_l_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob, 
                  restart, n_restart, n_search, increment,
                  omega_l, log_pi_l)
    }
  } else {
    if (reverse) {
      if (!pi_l_weighted) log_pi_l <- log_pi_l[T:1,,drop = FALSE]
      if (!pi_k_weighted) log_pi_k <- log_pi_k[T:1,,drop = FALSE]
      if (!pi_j_weighted) log_pi_j <- log_pi_j[T:1,,drop = FALSE]
      
      mich_vector(y[T:1], fit_intercept = TRUE, fit_scale = TRUE,
                  J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
                  pi_j_weighted, pi_l_weighted, pi_k_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob, 
                  restart, n_restart, n_search, increment,
                  omega_j, u_j, v_j, log_pi_j, 
                  omega_l, log_pi_l, 
                  u_k, v_k, log_pi_k)
    } else {
      mich_vector(y, fit_intercept, fit_scale,
                  J, L, K, J_auto, L_auto, K_auto, J_max, L_max, K_max,
                  pi_j_weighted, pi_l_weighted, pi_k_weighted,
                  tol, verbose, max_iter, reverse,
                  detect, merge_level, merge_prob, 
                  restart, n_restart, n_search, increment,
                  omega_j, u_j, v_j, log_pi_j, 
                  omega_l, log_pi_l, 
                  u_k, v_k, log_pi_k)
    }
  }
}

