
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

hsmuce_simulation <- function(T, K, C, min_space, B_l = min_space, B_r = min_space) {
  # sample change-points
  chp <- point_picker(T-B_l-B_r, K, min_space) + B_l
  # sample variances
  s <- c(1, 2^runif(K, -2, 2))
  # generate means and signal
  mu <- numeric(K + 1)
  blks <- diff(c(1, chp, T+1))
  y <- rnorm(blks[1])
  for (k in 1:K) {
    jump <- sqrt(C / min(blks[k+1] / s[k+1]^2, blks[k] / s[k]^2))
    mu[k+1] <- mu[k] + sample(c(-1,1), 1) * jump
    y <- c(y, rnorm(blks[k+1], mean = mu[k+1], sd = s[k+1]))
  }
  rep(mu, blks)
  return(list(y = y, 
              mu = mu, 
              s = s, 
              changepoints = chp, 
              mean_signal = rep(mu, blks), 
              var_signal = rep(s^2, blks)
              )
         )
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

