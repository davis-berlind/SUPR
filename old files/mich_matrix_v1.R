mich_matrix <- function(y, fit_intercept, fit_scale,
                        L, L_auto, L_max, pi_l_weighted,
                        tol, verbose, max_iter, reverse,
                        detect, merge_level, merge_prob, 
                        restart, n_restart, n_search, increment,
                        omega_l, pi_l) {
  
  #### set up ####
  # Flag for new model 
  refit <- FALSE
  
  # store original y
  y_raw <- y
  
  # calculate dimensions of y
  T <- nrow(y)
  d <- ncol(y)
  
  # estimate precision matrix 
  if (fit_scale) {
    y_diff <- diff(y)
    y_diff_norm <- sqrt(rowSums(y_diff^2))
    
    # remove outliers due to big mean changes
    y_diff <- y_diff[y_diff_norm <= quantile(y_diff_norm, p = 0.75) +  1.5 * IQR(y_diff_norm), ] 
    
    # estimate variance
    Sigma <- var(y_diff) / 2 
    Sigma_eigen <- eigen(Sigma)
    if (any(Sigma_eigen$values == 0)) {
      warning("Var(y) is singular. Consider removing collinear columns.")
      Sigma_eigen$values <- Sigma_eigen$values + 1e-5
    }
    Lambda_sqrt <- Sigma_eigen$vectors %*% diag(1 / sqrt(Sigma_eigen$values)) %*% t(Sigma_eigen$vectors)
    
    # rescale data
    y <- y %*% Lambda_sqrt
  }
  
  # initialize mu_0 
  if (fit_intercept) {
    mu_0 <- colMeans(y[1:ceiling(log(T)),,drop=FALSE])
  } else mu_0 <- rep(0.0, d)
  
  # initializing posterior parameters 
  log_pi_l <- log(pi_l)
  omega_bar_l <- omega_l + T - 1:T + 1
  log_omega_bar_l <- log(omega_bar_l)
  post_params <- list()
  if (L > 0) {
    for (l in 1:L) {
      post_params[[l]] <- list(pi_bar = rep(1 / T, T),
                               log_pi_bar = rep(0.0, T),
                               b_bar = matrix(0.0, nrow = T, ncol = d))
    }
  }
  
  #### fit model and merge ####
  merged <- FALSE
  while (!merged) {
    fit <- multi_mich_cpp(y, mu_0, fit_intercept, refit,
                          max_iter, tol, verbose = verbose & !L_auto,  
                          omega_l, log_pi_l, omega_bar_l, log_omega_bar_l,
                          post_params)
    
    merged <- TRUE
    refit <- TRUE
    
    if (L > 1) {
      # extract prob matrix
      pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])
      
      # test which columns actually contain change-points
      cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
      keep <- sapply(cred_sets, length) <= detect
      keep_mat <- matrix(FALSE, ncol = L, nrow = L)
      keep_mat[keep, keep] <- TRUE
      diag(keep_mat) <- FALSE
      
      # compute pairwise merge probabilities 
      merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
      diag(merge_prob_mat) <- 0
      
      merge_residual <- fit$residual
      
      while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
        merged <- FALSE
        L <- L - 1
        merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
        
        mu_bar_merge <- Reduce("+", lapply(merge_dex, function(l) multi_mu_bar_fn(post_params[[l]][["b_bar"]], post_params[[l]][["pi_bar"]])))
        merge_residual <- merge_residual + mu_bar_merge 
        
        merge_fit <- multi_mean_scp(merge_residual, omega_bar_l, log_omega_bar_l, log_pi_l[,merge_dex[1]])
        
        post_params[[merge_dex[1]]][["pi_bar"]] <- merge_fit$pi_bar
        pi_bar_l[,merge_dex[1]] <- merge_fit$pi_bar
        post_params[[merge_dex[1]]][["log_pi_bar"]] <- merge_fit$log_pi_bar
        post_params[[merge_dex[1]]][["b_bar"]] <- merge_fit$b_bar
        
        merge_residual <- merge_residual - multi_mu_bar_fn(post_params[[merge_dex[1]]][["b_bar"]], post_params[[merge_dex[1]]][["pi_bar"]]) 
        
        if (length(cred_set(pi_bar_l[,merge_dex[1]], level = merge_level)) > detect) {
          keep_mat[merge_dex[1], ] <- FALSE
          keep_mat[, merge_dex[1]] <- FALSE
        }
        
        merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_l[,merge_dex[1]] %*% pi_bar_l[,-merge_dex])
        merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
        merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
        keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
        
        post_params[merge_dex[2]] <- NULL
        pi_bar_l <- pi_bar_l[,-merge_dex[2]]
        if (L_auto) {
          pi_l <- pi_l[,-merge_dex[2], drop=FALSE]
          log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
        }     
      }
    }
    if (verbose & !merged) print("Merging components.")
  }
  
  # if components were merged out use auto procedure to increase to desired JLK
  merge_flag <- (L < L_max & !L_auto)
  if (merge_flag) increment <- 1 # ensures we don't overshoot number of components
  
  #### auto procedure ####
  last_restart <- ifelse(restart, 2, Inf)
  
  if (L_auto | merge_flag) {
    refit <- (L > 0)
    counter <- n_search
    elbo <- fit$elbo[length(fit$elbo)] # current value of elbo
    elbo_new <- elbo
    merged_at_l <- c()
    
    if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new))
    
    while (L < L_max & counter > 0) {
      # increment dimension of parameters
      for (i in 1:increment) {
        post_params[[L+i]] <- list(pi_bar = rep(1 / T, T),
                                   log_pi_bar = rep(0.0, T),
                                   b_bar = matrix(0.0, nrow = T, ncol = d))
        if (L > 1) post_params <- post_params[c(L+i, 1:(L+i-1))]
      }
      
      L <- L + increment
      if (L > 1 & L_auto) {
        pi_l <- cbind(matrix(pi_l[,1], nrow = T, ncol = increment), pi_l)
        log_pi_l <- cbind(matrix(log_pi_l[,1], nrow = T, ncol = increment), log_pi_l)
      }
      
      #### fit new model and merge ####
      merged <- FALSE
      while (!merged) {
        fit_new <- multi_mich_cpp(y, mu_0, fit_intercept, refit,
                                  max_iter, tol, verbose = FALSE,  
                                  omega_l, log_pi_l, omega_bar_l, log_omega_bar_l,
                                  post_params)
        
        merged <- TRUE
        refit <- TRUE
        
        if (L > 1 & !(L %in% merged_at_l)) {
          # keep track of where merges occur
          merged_at_l <- c(merged_at_l, L)
          
          # extract prob matrix
          pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])
          
          # test which columns actually contain change-points
          cred_sets <- apply(pi_bar_l, 2, cred_set, level = merge_level, simplify = FALSE)
          keep <- sapply(cred_sets, length) <= detect
          keep_mat <- matrix(FALSE, ncol = L, nrow = L)
          keep_mat[keep, keep] <- TRUE
          diag(keep_mat) <- FALSE
          
          # compute pairwise merge probabilities 
          merge_prob_mat <- t(pi_bar_l) %*% pi_bar_l
          diag(merge_prob_mat) <- 0
          
          merge_residual <- fit_new$residual
          
          while (L > 1 & any(merge_prob_mat[keep_mat] > merge_prob)) {
            merged <- FALSE
            L <- L - 1
            counter <- counter + 1
            merge_dex <- which(merge_prob_mat == max(merge_prob_mat[keep_mat]), arr.ind=TRUE)[1,]
            
            mu_bar_merge <- Reduce("+", lapply(merge_dex, function(l) multi_mu_bar_fn(post_params[[l]][["b_bar"]], post_params[[l]][["pi_bar"]])))
            merge_residual <- merge_residual + mu_bar_merge 
            
            merge_fit <- multi_mean_scp(merge_residual, omega_bar_l, log_omega_bar_l, log_pi_l[,merge_dex[1]])
            
            post_params[[merge_dex[1]]][["pi_bar"]] <- merge_fit$pi_bar
            pi_bar_l[,merge_dex[1]] <- merge_fit$pi_bar
            post_params[[merge_dex[1]]][["log_pi_bar"]] <- merge_fit$log_pi_bar
            post_params[[merge_dex[1]]][["b_bar"]] <- merge_fit$b_bar
            
            merge_residual <- merge_residual - multi_mu_bar_fn(post_params[[merge_dex[1]]][["b_bar"]], post_params[[merge_dex[1]]][["pi_bar"]]) 
            
            if (length(cred_set(pi_bar_l[,merge_dex[1]], level = merge_level)) > detect) {
              keep_mat[merge_dex[1], ] <- FALSE
              keep_mat[, merge_dex[1]] <- FALSE
            }
            
            if (length(merge_prob_mat[merge_dex[1],-merge_dex]) != length(c(pi_bar_l[,merge_dex[1]] %*% pi_bar_l[,-merge_dex]))) browser()
            
            merge_prob_mat[merge_dex[1],-merge_dex] <- c(pi_bar_l[,merge_dex[1]] %*% pi_bar_l[,-merge_dex])
            merge_prob_mat[-merge_dex, merge_dex[1]] <- merge_prob_mat[merge_dex[1],-merge_dex] 
            merge_prob_mat <- merge_prob_mat[-merge_dex[1], -merge_dex[1]]
            keep_mat <- keep_mat[-merge_dex[1], -merge_dex[1]]
            
            post_params[merge_dex[2]] <- NULL
            pi_bar_l <- pi_bar_l[,-merge_dex[2]]
            if (L_auto) {
              pi_l <- pi_l[,-merge_dex[2], drop=FALSE]
              log_pi_l <- log_pi_l[,-merge_dex[2], drop=FALSE]
            }     
          }
        }
      }
      
      # test if model improved or restart ####
      if (!refit) refit <- TRUE
      elbo_new <- fit_new$elbo[length(fit_new$elbo)]
      if (verbose) print(paste0("L = ", L, ": ELBO = ", elbo_new, 
                                "; Counter: ", counter))
      
      if (elbo_new > elbo | merge_flag) {
        elbo <- elbo_new
        fit <- fit_new
        counter <- n_search
        if (last_restart < Inf) restart <- TRUE 
      } else if (restart & L_auto & L - n_restart > last_restart) {
        refit <- FALSE
        restart <- FALSE
        last_restart <- L
        n_restart <- n_restart + 1
        
        L <- L - increment
        pi_l <- sapply(1:L, function(i) pi_l[,1])
        log_pi_l <- sapply(1:L, function(i) log_pi_l[,1])
        post_params <- list()
        for (l in 1:L) {
          post_params[[l]] <- list(pi_bar = rep(1 / T, T),
                                   log_pi_bar = rep(00.0, T),
                                   b_bar = matrix(0.0, nrow = T, ncol = d))
        }
      } else {
        counter <- counter - 1
      }
    }
  }
  
  #### reversing model ####
  if (reverse){
    if(verbose) print("reversing model")
    y_raw <- y_raw[T:1,]
    L <- fit$L
    
    # reversed residuals, variance, and intercepts
    r_bar <- fit$residual[T:1,] 
    mu_0 <- fit$mu[T,]

    # don't reverse weighted priors
    if (!pi_l_weighted) {
      pi_l <- pi_l[T:1,1:max(1,L),drop = FALSE]
      log_pi_l <- log_pi_l[T:1,1:max(1,L),drop = FALSE]
    } else {
      pi_l <- sapply(1:max(1,L), function(i) pi_l[,1])
      log_pi_l <- sapply(1:max(1,L), function(i) log_pi_l[,1])
    }
    
    # reversing mean components
    post_params <- fit$post_params
    pi_bar_l <- sapply(1:L, function(l) post_params[[l]][["pi_bar"]])

    for (l in seq_len(L)) {
      tau_l <- which.max(pi_bar_l[,l])
      mu_bar_l <- multi_mu_bar_fn(post_params[[l]][["b_bar"]], post_params[[l]][["pi_bar"]])
      r_bar_l <- fit$residual + mu_bar_l - matrix(colMeans(mu_bar_l[tau_l:T,,drop = FALSE]), nrow = T, ncol = d, byrow = TRUE)
      r_bar_l <- r_bar_l[T:1,]
      
      fit_scp <- multi_mean_scp(r_bar_l, omega_bar_l, log_omega_bar_l, log_pi_l[,l])
      
      post_params[[l]][["pi_bar"]] <- fit_scp$pi_bar
      post_params[[l]][["log_pi_bar"]] <- fit_scp$log_pi_bar
      post_params[[l]][["b_bar"]] <- fit_scp$b_bar
    }
    
    #### fit model ####
    fit <- multi_mich_cpp(y[T:1,], mu_0, fit_intercept, refit = TRUE,
                          max_iter, tol, verbose = verbose & !L_auto,  
                          omega_l, log_pi_l, omega_bar_l, log_omega_bar_l,
                          post_params)
  }
  
  #### return model ####
  fit$y <- y_raw
  
  # rescale
  if (fit_scale) {
    Sigma_sqrt <- Sigma_eigen$vectors %*% diag(sqrt(Sigma_eigen$values)) %*% t(Sigma_eigen$vectors)
    fit$mu_0 <- c(fit$mu_0 %*% Sigma_sqrt)
    fit$mu <- fit$mu %*% Sigma_sqrt
    for (l in seq_len(fit$L)) {
      fit$post_params[[l]][["b_bar"]] <- post_params[[l]][["b_bar"]]  %*% Sigma_sqrt
    }
  }
  
  if (fit$L > 0) fit$pi_bar_l <- sapply(1:fit$L, function(l) fit$post_params[[l]][["pi_bar"]])
  fit$Sigma <- Sigma
  
  return(fit)
}