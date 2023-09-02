#### Single Change Point Models ####

# The functions in this section return the posterior parameters for the single
# change point models, SMCP, SSCP, and SMSCP

#' Single Mean Change-Point Posterior
#'
#' Returns the posterior parameters for the model \eqn{y_t \sim N(\mu_t, \lambda_t^{-1})} 
#' where \eqn{\mu_t = b \mathbf{1}(t \geq \gamma)}, \eqn{b \sim \text{N}(0, \tau^{-1})}, 
#' and \eqn{\gamma \sim \text{Categorial}(\pi)}.
#'
#' @param y      A numeric vector. Length \eqn{T+B_r} vector of time-series
#'               observations.
#' @param lambda A numeric vector. Length \eqn{T+B_r} vector positive precision
#'               parameters for \eqn{y_t}.
#' @param tau    A scalar. Prior precision parameter for \eqn{b}. Must be
#'               positive.
#' @param pi     A numeric vector. Length \eqn{T} vector of prior probabilities
#'               that \eqn{\gamma = t}. Elements of \texttt{pi} must some to one.
#' @param B_r    An integer. Tail buffer parameter.
#' 
#' @return Posterior parameters for \eqn{b \;|\; \mathbf{y}}.
#' 
smcp <- function(y, lambda, tau, pi, B_r) {
  T <- length(y) - B_r
  tail <- (T+1):(T+B_r)
  tau <- tau + revcumsum(lambda[1:T])
  if (B_r > 0 ) tau <- tau + sum(lambda[tail])
  mu <- revcumsum(lambda[1:T] * y[1:T]) / tau
  if (B_r > 0 ) mu <- mu + sum(lambda[tail] * y[tail]) / tau
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
#' @param y      A numeric vector. Length \eqn{T+B_r} vector of time-series
#'               observations.
#' @param lambda A numeric vector. Length \eqn{T + B_r} vector positive precision
#'               parameters for \eqn{y_t}.
#' @param u      A scalar. Prior shape parameter for \eqn{s}. Must be positive.
#' @param v      A scalar. Prior rate parameter for \eqn{s}. Must be positive.
#' @param omega  A numeric vector. Length \eqn{T} vector of prior probabilities
#'               that \eqn{\alpha = t}. Elements of \texttt{omega} must some to one.
#' @param B_r    An integer. Tail buffer parameter.
#' 
#' @return Posterior parameters for \eqn{s \;|\; \mathbf{y}}.
#' 
sscp <- function(y, lambda, u, v, omega, B_r) {
  T <- length(y) - B_r
  tail <- (T+1):(T+B_r)
  y2_lambda <- c(0, cumsum(lambda[1:(T-1)] * y[1:(T-1)]^2))
  v <- v + (lambda[T] * y[T]^2 + y2_lambda[T] - y2_lambda) / 2
  if (B_r > 0) v <- v + sum(lambda[tail] * y[tail]^2) / 2
  u <- u + (T + B_r - 1:T + 1) / 2
  log_omega <- log(omega) + lgamma(u) - u * log(v) - y2_lambda / 2 
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
#' @param y      A numeric vector. Length \eqn{T + B_r} vector of time-series
#'               observations.
#' @param lambda A numeric vector. Length \eqn{T + B_r} vector positive precision
#'               parameters for \eqn{y_t}.
#' @param tau    A scalar. Prior precision parameter for \eqn{b}. Must be
#'               positive.
#' @param u      A scalar. Prior shape parameter for \eqn{s}. Must be positive.
#' @param v      A scalar. Prior rate parameter for \eqn{s}. Must be positive.
#' @param theta  A numeric vector. Length \eqn{T} vector of prior probabilities
#'               that \eqn{\beta = t}. Elements of \texttt{\theta} must some to one.
#' @param B_r    An integer. Tail buffer parameter.
#'   
#' @return Posterior parameters for \eqn{(b,s) \;|\; \mathbf{y}}.
#' 
smscp <- function(y, lambda, tau, u, v, theta, B_r) {
  T <- length(y) - B_r
  tail <- (T+1):(T+B_r)
  tau <- tau + revcumsum(lambda[1:T])
  if (B_r > 0 ) tau <- tau + sum(lambda[tail])
  mu <- revcumsum(lambda[1:T] * y[1:T]) / tau
  if (B_r > 0 ) mu <- mu + sum(lambda[tail] * y[tail]) / tau
  y2_lambda <- c(0, cumsum(lambda[1:(T-1)] * y[1:(T-1)]^2))
  v <- v - (tau * mu^2) / 2 + (lambda[T] * y[T]^2 + y2_lambda[T] - y2_lambda) / 2
  if (B_r > 0) v <- v + sum(lambda[tail] * y[tail]^2) / 2
  u <- u + (T + B_r - 1:T + 1) / 2
  log_theta <- log(theta) - 0.5 * log(tau) + lgamma(u) - u * log(v) - y2_lambda / 2 
  theta <- prop.table(exp(log_theta - max(log_theta)))
  return(list(tau = tau, mu = mu, v = v, theta = theta))
}

