mean_prior <- function(T, d) {
  log_pi <- rep(0, T)
  for (t in 1:(T-1)) {
    log_pi[t+1] <- log_pi[t] + 0.5 * d * (log(T-t) - log(T-t+1))
  }
  pi <- exp(log_pi)
  return(pi / sum(pi))
}

var_prior <- function(T) {
  log_pi <- rep(0, T)
  for (t in 1:(T-1)) {
    log_pi[t+1] <- sum(log_pi[t] + 0.5, 
                       lgamma(0.5 * (T-t+1)) - lgamma(0.5 * (T-t)),
                       0.5 * (T-t) * digamma(0.5 * (T-t)),
                       -0.5 * (T-t+1) * digamma(0.5 * (T-t+1)))
  }
  pi <- exp(log_pi)
  return(pi / sum(pi))
}

meanvar_prior <- function(T) {
  log_pi <- rep(0, T-1)
  for (t in 1:(T-2)) {
    log_pi[t+1] <- sum(log_pi[t] + 0.5, 
                       0.5 * (log(T-t) - log(T-t+1)),
                       lgamma(0.5 * (T-t+1)) - lgamma(0.5 * (T-t)),
                       0.5 * (T-t) * digamma(0.5 * (T-t-1)),
                       -0.5 * (T-t+1) * digamma(0.5 * (T-t)))
  }
  pi <- exp(log_pi)
  return(c(pi / sum(pi), 0))
}
