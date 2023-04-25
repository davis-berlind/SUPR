#### Single Change Point Models ####

# The functions in this section return the posterior parameters for the single
# change point models, SMCP, SSCP, SMSCP, and NG.

#' Normal-Gamma Model Posterior
#'
#' Returns the posterior parameters for the model \eqn{y_t \sim N(b, (s\lambda_t)^{-1})} 
#' where \eqn{(b, s) \sim \text{NormalGamma}(0, \tau, u, v)}.
#'
#' @param y      A numeric vector. Length \eqn{T} vector of time-series
#'               observations.
#' @param lambda A numeric vector. Length \eqn{T} vector of positive precision
#'               parameters for \eqn{y_t}.
#' @param tau    A scalar. Prior precision parameter for \eqn{b}. Must be
#'               positive.
#' @param v      A scalar. Prior rate parameter for \eqn{s}. Must be positive.
#'
#' @return Posterior parameters for \eqn{b,s \;|\; \mathbf{y}}.
#'
ng <- function(y, lambda, tau, v) {
  tau <- sum(lambda) + tau
  mu <- sum(y * lambda) / tau
  v <- v - (tau * mu^2) / 2 + sum(lambda * y^2) / 2
  return(list(tau = tau, mu = mu, v = v))
} 

#' Single Mean Change-Point Posterior
#'
#' Returns the posterior parameters for the model \eqn{y_t \sim N(\mu_t, \lambda_t^{-1})} 
#' where \eqn{\mu_t = b \mathbf{1}(t \geq \gamma)}, \eqn{b \sim \text{N}(0, \tau^{-1})}, 
#' and \eqn{\gamma \sim \text{Categorial}(\pi)}.
#'
#' @param y      A numeric vector. Length \eqn{T} vector of time-series
#'               observations.
#' @param lambda A numeric vector. Length \eqn{T} vector positive precision
#'               parameters for \eqn{y_t}.
#' @param tau    A scalar. Prior precision parameter for \eqn{b}. Must be
#'               positive.
#' @param pi     A numeric vector. Length \eqn{T} vector of prior probabilities
#'               that \eqn{\gamma = t}. Elements of \texttt{pi} must some to one.
#'   
#' @return Posterior parameters for \eqn{b \;|\; \mathbf{y}}.
#' 
smcp <- function(y, lambda, tau, pi) {
  tau <- tau + revcumsum(lambda) 
  mu <- revcumsum(lambda * y) / tau
  log_pi <- log(pi) - 0.5 * log(tau) + (tau * mu^2) / 2
  pi <- prop.table(exp(log_pi - max(log_pi)))
  return(list(tau = tau, mu = mu, pi = pi))
}

#' Single Scale Change-Point Posterior
#'
#' Returns the posterior parameters for the model \eqn{y_t \sim N(0, \lambda_t^{-1})} 
#' where \eqn{\lambda_t = s^{\mathbf{1}(t \geq \alpha)}}, \eqn{s \sim \text{Gamma}(u,v)}, 
#' and \eqn{\alpha \sim \text{Categorial}(\omega)}.
#'
#' @param y      A numeric vector. Length \eqn{T} vector of time-series
#'               observations.
#' @param u      A scalar. Prior shape parameter for \eqn{s}. Must be positive.
#' @param v      A scalar. Prior rate parameter for \eqn{s}. Must be positive.
#' @param omega  A numeric vector. Length \eqn{T} vector of prior probabilities
#'               that \eqn{\alpha = t}. Elements of \texttt{omega} must some to one.
#'   
#' @return Posterior parameters for \eqn{s \;|\; \mathbf{y}}.
#' 
sscp <- function(y, u, v, omega) {
  T <- length(y)
  v <- v + revcumsum(y^2) / 2
  u <- u + (T - 1:T + 1) / 2
  log_omega <- log(omega) + lgamma(u) - u * log(v) - cumsum(c(0, y[-T]^2)) / 2 
  omega <- prop.table(exp(log_omega - max(log_omega)))
  return(list(v = v, omega = omega))
}

#' Single Mean-Scale Change-Point Posterior
#'
#' Returns the posterior parameters for the model \eqn{y_t \sim N(\mu_t, (\tau_t \lambda_t)^{-1})} 
#' where \eqn{\mu_t = b \mathbf{1}(t \geq \gamma)}, \eq{\tau_t = \tau_t s^{\mathbf{1}(t \geq \gamma)}},  
#' \eqn{(b,s) \sim \text{NormalGamma}(0,\tau,u,v)}
#' and \eqn{\gamma \sim \text{Categorial}(\pi)}.
#'
#' @param y      A numeric vector. Length \eqn{T} vector of time-series
#'               observations.
#' @param lambda A numeric vector. Length \eqn{T} vector positive precision
#'               parameters for \eqn{y_t}.
#' @param tau    A scalar. Prior precision parameter for \eqn{b}. Must be
#'               positive.
#' @param u      A scalar. Prior shape parameter for \eqn{s}. Must be positive.
#' @param v      A scalar. Prior rate parameter for \eqn{s}. Must be positive.
#' @param pi     A numeric vector. Length \eqn{T} vector of prior probabilities
#'               that \eqn{\gamma = t}. Elements of \texttt{pi} must some to one.
#'   
#' @return Posterior parameters for \eqn{(b,s) \;|\; \mathbf{y}}.
#' 
smscp <- function(y, lambda, tau, u, v, pi) {
  T <- length(y)
  tau <- tau + revcumsum(lambda) 
  mu <- revcumsum(lambda * y) / tau
  v <- v - (tau * mu^2) / 2 + revcumsum(lambda * y^2) / 2
  u <- u + (T - 1:T + 1) / 2
  log_pi <- log(pi) - 0.5 * log(tau) + lgamma(u) - u * log(v) - cumsum(c(0, lambda[-T] * y[-T]^2)) / 2 
  pi <- prop.table(exp(log_pi - max(log_pi)))
  return(list(tau = tau, mu = mu, v = v, pi = pi))
}
