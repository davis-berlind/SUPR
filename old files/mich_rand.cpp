#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <random>
using namespace Rcpp;

inline int randWrapper(const int n) {
  return floor(unif_rand()*n); 
}

// [[Rcpp::export]]
Rcpp::NumericVector random_shuffle(Rcpp::NumericVector a) {
  // clone a into b to leave a alone
  Rcpp::NumericVector b = Rcpp::clone(a);
  int n = b.size();
  int j;
  
  // Fisher-Yates Shuffle Algorithm
  for (int i = 0; i < n - 1; i++) {
    j = i + randWrapper(n - i);
    std::swap(b[i], b[j]);
  }
  return b;
}

// [[Rcpp::export]]
double elbo_fn(int T, int B_l, int B_r, double mu_0, double lambda_0,
               NumericVector y_0, NumericVector r_tilde, NumericVector lambda_bar, NumericVector delta,
               NumericMatrix b_bar_j, NumericMatrix tau_bar_j, NumericMatrix u_bar_j, NumericMatrix v_bar_j, NumericMatrix pi_bar_j,
               double tau_j, double u_j, double v_j, NumericMatrix log_pi_j,
               NumericMatrix b_bar_l, NumericMatrix tau_bar_l, NumericMatrix pi_bar_l, 
               double tau_l, NumericMatrix log_pi_l,
               NumericMatrix u_bar_k, NumericMatrix v_bar_k, NumericMatrix pi_bar_k,
               double u_k, double v_k, NumericMatrix log_pi_k) {
  
  int J = pi_bar_j.ncol(), L = pi_bar_l.ncol(), K = pi_bar_k.ncol();
  
  // calculate E[log p(y)]
  double elbo = 0.5 * (T + B_l + B_r) * (std::log(lambda_0) - std::log(2) - std::log(M_PI));
  elbo -= 0.5 * Rcpp::sum(Rcpp::pow(y_0 - Rcpp::rep(mu_0, B_l), 2));
  elbo -= 0.5 * Rcpp::sum(lambda_bar * (Rcpp::pow(r_tilde, 2) + delta));
  
  for (int t = 0; t < T; t++) {
    for (int j = 0; j < J; j++) {
      // E[log lambda_j] component of E[log p(y)]
      elbo += 0.5 * (T - t + B_r) * (R::digamma(u_bar_j(t,j)) - std::log(v_bar_j(t,j))) * pi_bar_j(t,j);
      // E[log p_j] - E[log q_j]
      // adding log odds with numerical stability for small pi_bar_j values
      if (pi_bar_j(t,j) > 1e-100) elbo += pi_bar_j(t,j) * (log_pi_j(t,j) - std::log(pi_bar_j(t,j)));
      // mean component
      elbo += 0.5 * pi_bar_j(t,j) * (std::log(tau_j) - std::log(tau_bar_j(t,j)));
      elbo += 0.5 * pi_bar_j(t,j) * (1 - tau_j * (std::pow(b_bar_j(t,j), 2) * (u_bar_j(t,j) / v_bar_j(t,j)) + 1 / tau_bar_j(t,j)));
      // // precision component
      elbo += pi_bar_j(t,j) * (u_j * (std::log(v_j) - std::log(v_bar_j(t,j))) + std::lgamma(u_bar_j(t,j)) - std::lgamma(u_j));
      elbo += pi_bar_j(t,j) * ((u_j - u_bar_j(t,j)) * R::digamma(u_bar_j(t,j)) + u_bar_j(t,j) * (1 - v_j / v_bar_j(t,j)));
    }
    for (int l = 0; l < L; l++) {
      // E[log p_l] - E[log q_l]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_l(t,l) > 1e-100) elbo += pi_bar_l(t,l) * (log_pi_l(t,l) - std::log(pi_bar_l(t,l)));
      elbo += 0.5 * pi_bar_l(t,l) * (std::log(tau_l) - std::log(tau_bar_l(t,l)));
      elbo += 0.5 * pi_bar_l(t,l) * (1 - tau_l * (std::pow(b_bar_l(t,l), 2) + 1 / tau_bar_l(t,l)));
    }
    for (int k = 0; k < K; k++) {
      // E[log lambda_k] component of E[log p(y)]
      elbo += 0.5 * (T - t + B_r) * (R::digamma(u_bar_k(t,k)) - std::log(v_bar_k(t,k))) * pi_bar_k(t,k);
      // E[log p_j] - E[log q_j]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_k(t,k) > 1e-100) elbo += pi_bar_k(t,k) * (log_pi_k(t,k) - std::log(pi_bar_k(t,k)));
      elbo += pi_bar_k(t,k) * (u_k * (std::log(v_k) - std::log(v_bar_k(t,k))) + std::lgamma(u_bar_k(t,k)) - std::lgamma(u_k));
      elbo += pi_bar_k(t,k) * ((u_k - u_bar_k(t,k)) * R::digamma(u_bar_k(t,k)) + u_bar_k(t,k) * (1 - v_k / v_bar_k(t,k)));
    }
  }
  return elbo;
}

// [[Rcpp::export]]
NumericVector lambda_bar_fn(NumericVector u, NumericVector v, 
                            NumericVector prob, int B_r) {
  int T = prob.length();
  double fwd_sum = 0.0, rev_sum = 1.0;
  NumericVector lambda_bar (T + B_r, 0.0);
  for (int t = 0; t < T + B_r; t++) {
    if (t < T) {
      rev_sum -= prob[t];
      fwd_sum += (u[t] / v[t]) * prob[t];
      lambda_bar[t] = fwd_sum + rev_sum;
    } else {
      lambda_bar[t] = lambda_bar[T-1];
    }
  }
  return lambda_bar;
}

// [[Rcpp::export]]
NumericVector mu_bar_fn(NumericVector b, NumericVector prob, int B_r) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu_bar (T + B_r, 0.0);
  for (int t = 0; t < T + B_r; t++) {
    if (t < T) {
      fwd_sum += b[t] * prob[t];
      mu_bar[t] = fwd_sum;
    } else {
      mu_bar[t] = mu_bar[T-1];
    }
  }
  return mu_bar;
}

// [[Rcpp::export]]
NumericVector mu2_bar_fn(NumericVector b, NumericVector tau, 
                         NumericVector prob, int B_r) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu2_bar (T + B_r, 0.0);
  for (int t = 0; t < T + B_r; t++) {
    if (t < T) {
      fwd_sum += (std::pow(b[t], 2) + 1 / tau[t]) * prob[t];
      mu2_bar[t] = fwd_sum;
    } else {
      mu2_bar[t] = mu2_bar[T-1];
    }
  }
  return mu2_bar;
}

// [[Rcpp::export]]
NumericVector mu_lambda_fn(NumericVector b, NumericVector u, NumericVector v, 
                           NumericVector prob, int B_r) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu_lambda (T + B_r, 0.0);
  for (int t = 0; t < T + B_r; t++) {
    if (t < T) {
      fwd_sum += b[t] * (u[t] / v[t]) * prob[t];
      mu_lambda[t] = fwd_sum;
    } else {
      mu_lambda[t] = mu_lambda[T-1];
    }
  }
  return mu_lambda;
}

// [[Rcpp::export]]
NumericVector mu2_lambda_fn(NumericVector b, NumericVector tau, NumericVector u,
                            NumericVector v, NumericVector prob, int B_r) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu2_lambda (T + B_r, 0.0);
  for (int t = 0; t < T + B_r; t++) {
    if (t < T) {
      fwd_sum += (std::pow(b[t], 2) * (u[t] / v[t]) + 1 / tau[t]) * prob[t];
      mu2_lambda[t] = fwd_sum;
    } else {
      mu2_lambda[t] = mu2_lambda[T-1];
    }
  }
  return mu2_lambda;
}

// [[Rcpp::export]]
List mean_scp(NumericVector y, NumericVector lambda, double tau, 
              NumericVector log_pi, int B_r) {
  
  int T = y.length() - B_r;
  NumericVector b_bar (T), pi_bar (T), tau_bar (T);
  double tot = 0, ls = 0, ys = 0, log_pi_max = R_NegInf;
  
  for (int t = T + B_r - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    if (t < T) {
      tau_bar[t] = tau + ls;
      b_bar[t] = ys / tau_bar[t];
      pi_bar[t] = log_pi[t] - 0.5 * (std::log(tau_bar[t]) - tau_bar[t] * std::pow(b_bar[t], 2.0));
      if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
    }
  }
  
  for(int t = 0; t < T; t++) {
    pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
    tot += pi_bar[t];
  }
  return List::create(_["b_bar"] = b_bar, _["tau_bar"] = tau_bar, _["pi_bar"] = pi_bar / tot);
}

// [[Rcpp::export]]
List var_scp(NumericVector y, NumericVector tau, double u, 
             NumericVector v, NumericVector log_pi, int B_r) {
  
  int T = y.length() - B_r;
  double half = T / 2.0;
  NumericVector pi_bar (T, 0.0), v_bar (T), u_bar (T);
  double tot = 0, fs = 0, rs = 0, log_pi_max = R_NegInf;
  
  for (int t = T + B_r - 1; t >= 0; t--) {
    rs += std::pow(y[t], 2.0) * tau[t];
    if (t < T) {
      v_bar[t] = v[t] + 0.5 * rs;
      u_bar[t] = u + 0.5 * (T + B_r - t);
      pi_bar[t] += log_pi[t] + std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]);
      pi_bar[T-t-1] -= 0.5 * fs;
      if (t < half) {
        if (pi_bar[T-t-1] > log_pi_max) log_pi_max = pi_bar[T-t-1];
        if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
      }
      fs += std::pow(y[T-t-1], 2.0) * tau[T-t-1];
    }
  }
  
  for(int t = 0; t < T; t++) {
    pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
    tot += pi_bar[t];
  }
  return List::create(_["v_bar"] = v_bar, _["pi_bar"] = pi_bar / tot);
}

// [[Rcpp::export]]
List meanvar_scp(NumericVector y, NumericVector lambda, 
                 double tau, double u, NumericVector v, 
                 NumericVector log_pi, int B_r) {
  
  int T = y.length() - B_r;
  double half = T / 2.0;
  NumericVector b_bar (T), pi_bar (T, 0.0), u_bar (T), v_bar (T), tau_bar (T);
  double tot = 0, fs = 0, rs = 0, ls = 0, ys = 0, log_pi_max = R_NegInf;
  
  for (int t = T + B_r - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    rs += std::pow(y[t], 2.0) * lambda[t];
    if (t < T) {
      tau_bar[t] = tau + ls;
      b_bar[t] = ys / tau_bar[t];
      v_bar[t] = v[t] + 0.5 * (rs - std::pow(b_bar[t], 2.0) * tau_bar[t]);
      u_bar[t] = u + 0.5 * (T + B_r - t);
      pi_bar[t] += log_pi[t] + std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]) - 0.5 * std::log(tau_bar[t]);
      pi_bar[T-t-1] -= 0.5 * fs;
      if (t < half) {
        if (pi_bar[T-t-1] > log_pi_max) log_pi_max = pi_bar[T-t-1];
        if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
      }
      fs += std::pow(y[T-t-1], 2.0) * lambda[T-t-1];
    }
  }
  
  for(int t = 0; t < T; t++) {
    pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
    tot += pi_bar[t];
  }
  return List::create(_["b_bar"] = b_bar, _["tau_bar"] = tau_bar,
                      _["v_bar"] = v_bar, _["pi_bar"] = pi_bar / tot);
}

// [[Rcpp::export]]
List mich_cpp_rand(NumericVector y_0, NumericVector y, 
              int J, int L, int K, double mu_0, double lambda_0,
              int B_l, int B_r, bool fit_intercept, bool fit_scale, 
              double max_iter, bool verbose, double tol, 
              double tau_j, double u_j, double v_j, NumericMatrix log_pi_j,
              NumericMatrix pi_bar_j, NumericMatrix b_bar_j, NumericMatrix tau_bar_j, NumericMatrix u_bar_j, NumericMatrix v_bar_j,
              double tau_l, NumericMatrix log_pi_l,
              NumericMatrix pi_bar_l, NumericMatrix b_bar_l, NumericMatrix tau_bar_l, 
              double u_k, double v_k, NumericMatrix log_pi_k,
              NumericMatrix pi_bar_k, NumericMatrix u_bar_k, NumericMatrix v_bar_k) {
  
  int T = y.length() - B_r;
  
  double y_0_sum = Rcpp::sum(y_0);
  
  // initialize residual 
  NumericVector r_tilde = clone(y);
  
  // initialize precision
  NumericVector lambda_bar (T + B_r, lambda_0);
  
  // initialize correction term
  NumericVector delta (T + B_r, 0.0);
  
  // initialize J expected mean-variance parameters
  NumericMatrix mu_lambda_j (T + B_r, J), mu2_lambda_j (T + B_r, J), lambda_bar_j (T + B_r, J);
  for (int j = 0; j < J; j++) {
    mu_lambda_j(_,j) = mu_lambda_fn(b_bar_j(_,j), u_bar_j(_,j), v_bar_j(_,j), pi_bar_j(_,j), B_r);
    mu2_lambda_j(_,j) = mu2_lambda_fn(b_bar_j(_,j), tau_bar_j(_,j),  u_bar_j(_,j), v_bar_j(_,j), pi_bar_j(_,j), B_r);
    lambda_bar_j(_,j) = lambda_bar_fn(u_bar_j(_,j), v_bar_j(_,j), pi_bar_j(_,j), B_r);
  }

  // initialize L expected mean parameters
  NumericMatrix mu_bar_l (T + B_r, L), mu2_bar_l (T + B_r, L);
  for (int l = 0; l < L; l++) {
    mu_bar_l(_,l) = mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l), B_r);
    mu2_bar_l(_,l) = mu2_bar_fn(b_bar_l(_,l), tau_bar_l(_,l), pi_bar_l(_,l), B_r);
  }

  // initialize K expected variance parameters
  NumericMatrix lambda_bar_k (T + B_r, K);
  for (int k = 0; k < K; k++) {
    lambda_bar_k(_,k) = lambda_bar_fn(u_bar_k(_,k), v_bar_k(_,k), pi_bar_k(_,k), B_r);
  }

  // initialize combined residual, variance, and correction term
  for (int t = 0; t < T + B_r; t++) {
    r_tilde[t] -= mu_0;

    for (int j = 0; j < J; j++) {
      r_tilde[t] -= mu_lambda_j(t,j) / lambda_bar_j(t,j);
      lambda_bar[t] *= lambda_bar_j(t,j);
      delta[t] -= (std::pow(mu_lambda_j(t,j), 2) / lambda_bar_j(t,j) - mu2_lambda_j(t,j)) / lambda_bar_j(t,j);
    }

    for (int l = 0; l < L; l++) {
      r_tilde[t] -= mu_bar_l(t,l);
      delta[t] -= std::pow(mu_bar_l(t,l), 2) - mu2_bar_l(t,l);
    }

    for (int k = 0; k < K; k++) {
      lambda_bar[t] *= lambda_bar_k(t,k);
    }
  }

  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);

  // initialize model fits
  List mean_scp_fit, var_scp_fit, meanvar_scp_fit;

  // initialize correction priors
  NumericVector v_tilde (T, 0.0), log_pi_tilde (T, 0.0);
  
  // vb algorithm
  NumericVector order_j (J), order_l (L), order_k (K);
  int id_j, id_l, id_k;
  
  if (J > 0) order_j = Rcpp::seq(0, J-1);
  if (L > 0) order_l = Rcpp::seq(0, L-1);
  if (K > 0) order_k = Rcpp::seq(0, K-1);
  
  int iter = 0;
  
  while (iter < max_iter) {

    // updating q(b_j, s_j, gamma_j)
    order_j = random_shuffle(order_j); // randomize update order
    for (int j = 0; j < J; j++) {
      id_j = order_j[j];
       
      v_tilde[T-1] = 0.0;    // initialize corrected prior rate
      log_pi_tilde[0] = 0.0; // initialize corrected prior change-point probs

      for (int t = T + B_r -1; t >= 0; t--) {
        // add back j^th component of residual
        r_tilde[t] += mu_lambda_j(t, id_j) / lambda_bar_j(t, id_j);
        // divide out j^th component of precision
        lambda_bar[t] /= lambda_bar_j(t, id_j);
        // subtract j^th component of variance correction
        delta[t] += (std::pow(mu_lambda_j(t, id_j), 2) / lambda_bar_j(t, id_j) - mu2_lambda_j(t, id_j)) / lambda_bar_j(t, id_j);

        // calculate corrected priors
        if (t >= T-1) {
          v_tilde[T-1] += 0.5 * lambda_bar[t] * delta[t];
        } else {
          v_tilde[t] = 0.5 * lambda_bar[t] * delta[t] + v_tilde[t+1];
          log_pi_tilde[T-t-1] = -0.5 * lambda_bar[T-t-2] * delta[T-t-2] + log_pi_tilde[T-t-2];
        }
      }

      // fit single meanvar_scp model on residuals
      meanvar_scp_fit = meanvar_scp(r_tilde, lambda_bar, tau_j, u_j, v_tilde + v_j, log_pi_tilde + log_pi_j(_, id_j), B_r);

      // store updated posterior parameters
      pi_bar_j(_, id_j) =  as<NumericVector>(meanvar_scp_fit["pi_bar"]);
      b_bar_j(_, id_j) =  as<NumericVector>(meanvar_scp_fit["b_bar"]);
      tau_bar_j(_, id_j) =  as<NumericVector>(meanvar_scp_fit["tau_bar"]);
      v_bar_j(_, id_j) =  as<NumericVector>(meanvar_scp_fit["v_bar"]);

      // update j^th component of mean and precision
      mu_lambda_j(_, id_j) = mu_lambda_fn(b_bar_j(_, id_j), u_bar_j(_, id_j), v_bar_j(_, id_j), pi_bar_j(_, id_j), B_r);
      mu2_lambda_j(_, id_j) = mu2_lambda_fn(b_bar_j(_, id_j), tau_bar_j(_, id_j),  u_bar_j(_, id_j), v_bar_j(_, id_j), pi_bar_j(_, id_j), B_r);
      lambda_bar_j(_, id_j) = lambda_bar_fn(u_bar_j(_, id_j), v_bar_j(_, id_j), pi_bar_j(_, id_j), B_r);

      for (int t = 0; t < T + B_r; t++) {
        // subtract out j^th component of residual
        r_tilde[t] -= mu_lambda_j(t, id_j) / lambda_bar_j(t, id_j);
        // multiply back j^th component or precision
        lambda_bar[t] *= lambda_bar_j(t, id_j);
        // add back j^th component of variance correction
        delta[t] -= (std::pow(mu_lambda_j(t, id_j), 2) / lambda_bar_j(t, id_j) - mu2_lambda_j(t, id_j)) / lambda_bar_j(t, id_j);
      }
    }

    // updating q(b_l, gamma_l)
    order_l = random_shuffle(order_l); // randomize update order
    for (int l = 0; l < L; l++) {
      id_l = order_l[l];
      for (int t = T + B_r - 1; t >= 0; t--) {
        // add back l^th component of residual
        r_tilde[t] += mu_bar_l(t, id_l);
        // subtract l^th component of variance correction
        delta[t] += std::pow(mu_bar_l(t, id_l), 2) - mu2_bar_l(t, id_l);
      }

      // fit single mean_scp model on residuals
      mean_scp_fit = mean_scp(r_tilde, lambda_bar, tau_l, log_pi_l(_, id_l), B_r);

      // store updated posterior parameters
      pi_bar_l(_, id_l) =  as<NumericVector>(mean_scp_fit["pi_bar"]);
      b_bar_l(_, id_l) =  as<NumericVector>(mean_scp_fit["b_bar"]);
      tau_bar_l(_, id_l) =  as<NumericVector>(mean_scp_fit["tau_bar"]);

      // update l^th component of mean
      mu_bar_l(_, id_l) = mu_bar_fn(b_bar_l(_, id_l), pi_bar_l(_, id_l), B_r);
      mu2_bar_l(_, id_l) = mu2_bar_fn(b_bar_l(_, id_l), tau_bar_l(_, id_l), pi_bar_l(_, id_l), B_r);

      for (int t = 0; t < T + B_r; t++) {
        // subtract out l^th component of residual
        r_tilde[t] -= mu_bar_l(t, id_l);
        // add back l^th component of variance correction
        delta[t] -= std::pow(mu_bar_l(t, id_l), 2) - mu2_bar_l(t, id_l);
      }
    }

    // updating q(s_k, gamma_k)
    order_k = random_shuffle(order_k); // randomize update order
    for (int k = 0; k < K; k++) {
      id_k = order_k[k];

      v_tilde[T-1] = 0.0;    // initialize corrected prior rate
      log_pi_tilde[0] = 0.0; // initialize corrected prior change-point probs
      for (int t = T + B_r - 1; t >= 0; t--) {
        // divide out k^th component of precision
        lambda_bar[t] /= lambda_bar_k(t, id_k);

        // calculate corrected priors
        if (t >= T-1) {
          v_tilde[T-1] += 0.5 * lambda_bar[t] * delta[t];
        } else {
          v_tilde[t] = 0.5 * lambda_bar[t] * delta[t] + v_tilde[t+1];
          log_pi_tilde[T-t-1] = -0.5 * lambda_bar[T-t-2] * delta[T-t-2] + log_pi_tilde[T-t-2];
        }
      }

      // fit single var_scp model on residuals
      var_scp_fit = var_scp(r_tilde, lambda_bar, u_k, v_tilde + v_k, log_pi_tilde + log_pi_k(_, id_k), B_r);

      // store updated posterior parameters
      pi_bar_k(_, id_k) = as<NumericVector>(var_scp_fit["pi_bar"]);
      v_bar_k(_, id_k) = as<NumericVector>(var_scp_fit["v_bar"]);

      // updated k^th component precision
      lambda_bar_k(_, id_k) = lambda_bar_fn(u_bar_k(_, id_k), v_bar_k(_, id_k), pi_bar_k(_, id_k), B_r);

      for (int t = 0; t < T + B_r; t++) {
        // multiply back k^th component or precision
        lambda_bar[t] *= lambda_bar_k(t, id_k);
      }
    }

    // updating mu_0 and lambda_0
    if (fit_intercept) {
      r_tilde += Rcpp::rep(mu_0, T + B_r);
      mu_0 = (y_0_sum + Rcpp::sum(lambda_bar * r_tilde)) / (B_l + Rcpp::sum(lambda_bar));
      r_tilde += Rcpp::rep(-mu_0, T + B_r);
    }

    if (fit_scale) {
      for (int t = 0; t < T + B_r; t++) {
        lambda_bar[t] /= lambda_0;
      }

      lambda_0 = (T + B_l + B_r) / (Rcpp::sum(Rcpp::pow(y_0 - Rcpp::rep(mu_0, B_l), 2)) + Rcpp::sum(lambda_bar * (Rcpp::pow(r_tilde, 2) + delta)));

      for (int t = 0; t < T + B_r; t++) {
        lambda_bar[t] *= lambda_0;
      }
    }

    iter++;
    elbo.push_back(0.0);
    elbo[iter] = elbo_fn(T, B_l, B_r, mu_0, lambda_0,
                         y_0, r_tilde, lambda_bar, delta,
                         b_bar_j, tau_bar_j, u_bar_j, v_bar_j, pi_bar_j,
                         tau_j, u_j, v_j, log_pi_j,
                         b_bar_l, tau_bar_l, pi_bar_l,
                         tau_l, log_pi_l,
                         u_bar_k, v_bar_k, pi_bar_k,
                         u_k, v_k, log_pi_k);
    if (verbose & (iter % 1000 == 0)) Rcout << "Iteration: " << iter << " elbo: " << elbo[iter] << "\n";
    if (elbo[iter] - elbo[iter - 1] < tol) break;
  }
  
  // construct mean signal
  NumericVector mu (T + B_r, mu_0);
  for (int j = 0; j < J; j++) {
    mu += mu_bar_fn(b_bar_j(_,j), pi_bar_j(_,j), B_r);
  }
  for (int l = 0; l < L; l++) {
    mu += mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l), B_r);
  }

  for (int t = 0; t < B_l; t++) {
    r_tilde.push_front(y_0[t] - mu_0);
    mu.push_front(mu_0); // construct mean signal
    lambda_bar.push_front(lambda_0); // construct scale signal
  }

  // creating lists of posterior parameters
  List J_model, L_model, K_model;
  
    List result = List::create(_["residual"] = r_tilde, _["mu"] = mu,
                             _["lambda"] = lambda_bar,
                             _["converged"] = (max_iter > iter),
                             _["elbo"] = elbo,
                             _["mu_0"] = mu_0, _["lambda_0"] = lambda_0,
                             _["J"] = J, _["L"] = L, _["K"] = K);

  if (J > 0) {
    J_model =  List::create(_["pi_bar"] = pi_bar_j, _["b_bar"] = b_bar_j,
                            _["tau_bar"] = tau_bar_j, _["v_bar"] = v_bar_j,
                            _["u_bar"] = u_bar_j);
    result["J_model"] = J_model;
  }
  if (L > 0) {
    L_model =  List::create(_["pi_bar"] = pi_bar_l, _["b_bar"] = b_bar_l,
                            _["tau_bar"] = tau_bar_l);
    result["L_model"] = L_model;
  }
  if (K > 0) {
    K_model =  List::create(_["pi_bar"] = pi_bar_k, _["v_bar"] = v_bar_k,
                            _["u_bar"] = u_bar_k);
    result["K_model"] = K_model;
  }

  return result;
}