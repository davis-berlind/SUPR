serial_cov <- function(d, rho, s) {
  cov_mat <- rho^abs(outer(1:d,1:d,"-"))
  for (i in 1:d) {
    for (j in i:d) {
      cov_mat[i,j] <- cov_mat[i,j] * s[i] * s[j]
      cov_mat[j,i] <- cov_mat[i,j]
    }
  }
  return(cov_mat)
} 

rlaplace <- function(n, mu = 0, b = 1) {
  u <- runif(n, min = -0.5, max = 0.5)
  x <- mu - b * sign(u) * log(1 - 2 * abs(u))
  return(x)
}

#' Random Change-Point Generator.
#' 
#' Randomly draws K change-points uniformly distributed between 1:T subject
#' to a minimum spacing criterion
#'
#' @param T         An integer. Number of observations.
#' @param K         An integer. Number of change-points.
#' @param min_space An integer. The minimum spacing criterion between 
#'                  change-points
#'
#' @return K change-points in 1:T.
#'
point_picker <- function(T, K, min_space) {
  if (min_space * (K-1) > T) {
    stop("T must be > (K-1) * min_space")
  }
  valid <- 1:T # current valid change-point locations
  picked <- c() # initialized sampled locations
  for (k in 1:K) {
    # check enough remaining valid points
    if (length(valid) < K - k + 1) return(point_picker(T, K, min_space))
    # sample change point from valid locations
    if (length(valid) == 1) picked <- c(picked, valid)
    else picked <- c(picked, sample(valid, size = 1))
    # update set of valid locations
    valid <- valid[abs(valid - picked[k]) >= min_space]
  }
  return(picked[order(picked)])
}

hsmuce_simulation <- function(T, K, C, min_space, family = "gaussian", df = NULL, theta = NULL) {
  # sample change-points
  chp <- point_picker(T - 2 * min_space, K, min_space) + min_space
  # sample variances
  s <- c(1, 2^runif(K, -2, 2))
  # generate means and signal
  mu <- numeric(K + 1)
  blks <- diff(c(1, chp, T+1))
  jumps <- sapply(2:(K+1), function(k) sqrt(C / min(blks[k] / s[k]^2, blks[k-1] / s[k-1]^2)))
  mu[2:(K+1)] <- cumsum(jumps * sample(c(-1,1), replace = TRUE, K))
  mean_signal = rep(mu, blks) 
  var_signal = rep(s^2, blks)
  if (family == "gaussian") error <- sqrt(var_signal) * rnorm(T)
  else if (family == "laplace") error <- sqrt(var_signal / 2) * rlaplace(T)
  else if (family == "t") error <- sqrt(var_signal * (df - 2) / df) * rt(T, df)
  else if (family == "MA") {
    error <- rnorm(T + 2)
    error <- sqrt(var_signal / (1 + theta^2 + theta^4)) * sapply(1:T, function(t) error[t+2] + theta * error[t+1] + theta^2 * error[t])
  }
  
  y <- mean_signal + error
  
  return(list(y = y, 
              mu = mu, 
              s = s, 
              changepoints = chp, 
              mean_signal = mean_signal, 
              var_signal = var_signal)
  )
}

multi_simulation <- function(T, K, d, p, rho, C, min_space, adapt) {
  # number of active series
  d_0 <- floor(p * d)
  # sample change-points
  chp <- point_picker(T - 2 * min_space, K, min_space) + min_space
  # sample active indices
  active <- sample(1:d, d_0, replace = FALSE)
  # generate means and signal
  s <- 2^runif(d, -2, 2)
  mu <- matrix(0, ncol = d, nrow = K + 1)
  blks <- diff(c(1, chp, T+1))
  mean_signal <- matrix(0, nrow = blks[1], ncol = d)
  for (k in 2:(K+1)) {
    jump <- rnorm(d_0, sqrt(C * s[active] / (ifelse(adapt, d_0, 1) * min(blks[k], blks[k-1]))), sd = 0.1)
    mu[k,active] <- mu[k-1,active] + jump * sample(c(-1,1), d_0, replace = TRUE)
    mean_signal <- rbind(mean_signal, matrix(mu[k,], nrow = blks[k], ncol = d, byrow = TRUE))
  }
  
  Sigma <- serial_cov(d, rho = 0.9, s)
  Sigma_eigen <- eigen(Sigma)
  
  Y <- mean_signal + sapply(1:d, function(x) rnorm(T)) %*% Sigma_eigen$vectors %*% diag(sqrt(Sigma_eigen$values)) %*% t(Sigma_eigen$vectors)
  
  return(list(Y = Y, 
              mu = mu, 
              Sigma = Sigma,
              changepoints = chp, 
              mean_signal = mean_signal))
}

fpsle <- function(true_cp, est_cp) {
  true_cp <- unique(true_cp[order(true_cp)])
  est_cp <- unique(est_cp[order(est_cp)])
  K <- length(est_cp)
  mid <- est_cp[-K] + diff(est_cp) / 2
  nbs <- as.numeric(cut(mid, true_cp, include.lowest = TRUE))
  error <- sum(abs(est_cp[-1] - true_cp[nbs+1]) + abs(est_cp[-K] - true_cp[nbs])) / (2*(K-1))
  return(error)
}

fnsle <- function(true_cp, est_cp) {
  fpsle(est_cp, true_cp)
}

hausdorff <- function(true_cp, est_cp) {
  max(sapply(true_cp, function(x) min(abs(x - est_cp))))
}

