
#' Random Change-Point Generator.
#' 
#' Randomly draws K change-points uniformly distributed between 2:(T-1) subject
#' to a minimum spacing criterion
#'
#' @param T         An integer. Number of observations.
#' @param K         An integer. Number of change-points.
#' @param min_space An integer. The minimum spacing criterion between 
#'                  change-points
#'
#' @return K change-points between 2:(T-1).
#'
point_picker <- function(T, K, min_space) {
  valid <- 2:(T-2) # current valid changepoint locations
  picked <- c()    # initialized sampled locations
  for (k in 1:K) {
    # sample change point from valid locations
    picked <- c(picked, sample(valid, size = 1) )       
    # update set of valid locations
    valid <- valid[valid < picked[k] - min_space | valid > picked[k] + min_space]
    if (length(valid) == 0) {
      picked <- point_picker(T, K, min_space) 
      break
    }
  }
  return(picked)
}

hsmuce_simulation <- function(T, K, C, min_space) {
  # sample changepoints
  chp <- point_picker(T, K, min_space)
  chp <- chp[order(chp)]
  # sample variances
  s <- c(1, 2^runif(K, -2, 2))
  # generate means and signal
  mu <- numeric(K + 1)
  y <- c()
  blks <- diff(c(1, chp, T+1))
  for (k in 1:(K+1)) {
    jump <- sqrt(C / min(blks[k+1] / s[k+1]^2, blks[k] / s[k]^2))
    if (k < K+1) mu[k+1] <- mu[k] + sample(c(-1,1), 1) * jump
    y <- c(y, rnorm(blks[k], mean = mu[k], sd = s[k]))
  }
  return(list(y = y, 
              mu = mu, 
              s = s, 
              changepoints = chp, 
              mean_signal = rep(mu, blks), 
              var_signal = rep(s^2, blks))
         )
}

fpsle <- function(true_cp, est_cp, T) {
  true_cp <- unique(c(1, true_cp[order(true_cp)], T+1))
  est_cp <- unique(c(1, est_cp[order(est_cp)], T+1))
  K <- length(est_cp)
  mid <- est_cp[-K] + diff(est_cp) / 2
  nbs <- as.numeric(cut(mid, true_cp, include.lowest = TRUE))
  error <- sum(abs(est_cp[-1] - true_cp[nbs+1]) + abs(est_cp[-K] - true_cp[nbs])) / (2*(K-1))
  return(error)
}

fnsle <- function(true_cp, est_cp, T) {
  fpsle(est_cp, true_cp, T)
}

hausdorff <- function(true_cp, est_cp, T) {
  if (length(est_cp) == 0) return(T)
  else return(max(sapply(true_cp, function(x) min(abs(x - est_cp)))))
}