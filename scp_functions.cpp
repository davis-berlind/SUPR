#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List smcp(NumericVector y, NumericVector lambda, double tau, NumericVector pi, int B_r) {
  int T = y.length() - B_r;
  NumericVector mu(T), pi_bar(T), tau_bar(T);
  double tot = 0, ls = 0, ys = 0, log_pi_max = R_NegInf;
  
  for (int t = y.length() - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    if (t < T) {
      tau_bar[t] = tau + ls;
      mu[t] = ys / tau_bar[t];
      pi_bar[t] = std::log(pi[t]) - 0.5 * (std::log(tau_bar[t]) - tau_bar[t] * std::pow(mu[t], 2.0));
      if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
    }
  }
  
  for(int t = 0; t < T; t++) {
    pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
    tot += pi_bar[t];
  }
  
  return List::create(_["mu"] = mu, _["tau"] = tau_bar, _["pi"] = pi_bar / tot);
}

// [[Rcpp::export]]
List sscp(NumericVector y, NumericVector tau, double u, NumericVector v, NumericVector pi, int B_r) {
  int T = y.length() - B_r;
  NumericVector pi_bar(T), v_bar(T), u_bar(T);
  double tot = 0, fs = 0, rs = 0, log_pi_max = R_NegInf;
  
  for(int t = 0; t < T; t++) {
    fs += std::pow(y[t], 2.0) * tau[t];
  }
  
  for (int t = y.length() - 1; t >= 0; t--) {
    rs += std::pow(y[t], 2.0) * tau[t];
    if (t < T) {
      fs -= std::pow(y[t], 2.0) * tau[t];
      v_bar[t] = v[t] + 0.5 * rs;
      u_bar[t] = u + 0.5 * (T + B_r - t);
      pi_bar[t] = std::log(pi[t]) + std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]) - 0.5 * fs;
      if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
    }
  }

  for(int t = 0; t < T; t++) {
    pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
    tot += pi_bar[t];
  }
  
  return List::create(_["u"] = u_bar, _["v"] = v_bar, _["pi"] = pi_bar / tot);
}

// [[Rcpp::export]]
List smscp(NumericVector y, NumericVector lambda, double tau, double u, NumericVector v, NumericVector pi, int B_r) {
  int T = y.length() - B_r;
  NumericVector mu(T), pi_bar(T), u_bar(T), v_bar(T), tau_bar(T);
  double tot = 0, fs = 0, rs = 0, ls = 0, ys = 0, log_pi_max = R_NegInf;

  for(int t = 0; t < T; t++) {
    fs += std::pow(y[t], 2.0) * lambda[t];
  }
  
  for (int t = y.length() - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    rs += std::pow(y[t], 2.0) * lambda[t];
    if (t < T) {
      tau_bar[t] = tau + ls;
      mu[t] = ys / tau_bar[t];
      fs -= std::pow(y[t], 2.0) * lambda[t];
      v_bar[t] = v[t] + 0.5 * (rs - std::pow(mu[t], 2.0) * tau_bar[t]);
      u_bar[t] = u + 0.5 * (T + B_r - t);
      pi_bar[t] = std::log(pi[t]) + std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]) - 0.5 * (fs + std::log(tau_bar[t]));
      if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
    }
  }
  
  for(int t = 0; t < T; t++) {
    pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
    tot += pi_bar[t];
  }
  
  return List::create(_["mu"] = mu, _["tau"] = tau_bar, _["u"] = u_bar, _["v"] = v_bar, _["pi"] = pi_bar / tot);
}
