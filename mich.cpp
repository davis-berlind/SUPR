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

const double log_2_pi = std::log(2 * M_PI);

// [[Rcpp::export]]
double elbo_fn(int T, double mu_0, double lambda_0,
               NumericVector r_tilde, NumericVector lambda_bar, NumericVector delta,
               NumericMatrix b_bar_j, NumericMatrix omega_bar_j, 
               NumericVector u_bar_j, NumericMatrix v_bar_j, NumericMatrix pi_bar_j, NumericMatrix log_pi_bar_j,
               NumericVector lgamma_u_bar_j, NumericVector digamma_u_bar_j, 
               double omega_j, double u_j, double v_j, 
               double log_omega_j, double log_u_j, double lgamma_u_j, double log_v_j, NumericMatrix log_pi_j,
               NumericMatrix b_bar_l, NumericMatrix omega_bar_l, NumericMatrix pi_bar_l, NumericMatrix log_pi_bar_l,
               double omega_l, double log_omega_l, NumericMatrix log_pi_l,
               NumericVector u_bar_k, NumericMatrix v_bar_k, NumericMatrix pi_bar_k,  
               NumericVector lgamma_u_bar_k, NumericVector digamma_u_bar_k, NumericMatrix log_pi_bar_k,
               double u_k, double v_k, 
               double log_u_k, double lgamma_u_k, double log_v_k, NumericMatrix log_pi_k) {
  
  int J = pi_bar_j.ncol(), L = pi_bar_l.ncol(), K = pi_bar_k.ncol();
  
  // calculate E[log p(y)]
  double elbo = 0.5 * T * (std::log(lambda_0) - log_2_pi);
  elbo -= 0.5 * Rcpp::sum(lambda_bar * (Rcpp::pow(r_tilde, 2) + delta));
  
  for (int t = 0; t < T; t++) {
    for (int j = 0; j < J; j++) {
      // E[log lambda_j] component of E[log p(y)]
      elbo += 0.5 * (T - t) * (digamma_u_bar_j[t] - std::log(v_bar_j(t,j))) * pi_bar_j(t,j);
      // E[log p_j] - E[log q_j]
      // adding log odds with numerical stability for small pi_bar_j values
      if (pi_bar_j(t,j) > 1e-20) elbo += pi_bar_j(t,j) * (log_pi_j(t,j) - log_pi_bar_j(t,j));
      // mean component
      elbo += 0.5 * pi_bar_j(t,j) * (log_omega_j - std::log(omega_bar_j(t,j)));
      elbo += 0.5 * pi_bar_j(t,j) * (1 - omega_j * (b_bar_j(t,j) * b_bar_j(t,j) * (u_bar_j[t] / v_bar_j(t,j)) + 1 / omega_bar_j(t,j)));
      // // precision component
      elbo += pi_bar_j(t,j) * (u_j * (log_v_j - std::log(v_bar_j(t,j))) + lgamma_u_bar_j[t] - lgamma_u_j);
      elbo += pi_bar_j(t,j) * ((u_j - u_bar_j[t]) * digamma_u_bar_j[t] + u_bar_j[t] * (1 - v_j / v_bar_j(t,j)));
    }
    for (int l = 0; l < L; l++) {
      // E[log p_l] - E[log q_l]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_l(t,l) > 1e-20) elbo += pi_bar_l(t,l) * (log_pi_l(t,l) - log_pi_bar_l(t,l));
      elbo += 0.5 * pi_bar_l(t,l) * (log_omega_l - std::log(omega_bar_l(t,l)));
      elbo += 0.5 * pi_bar_l(t,l) * (1 - omega_l * (b_bar_l(t,l) * b_bar_l(t,l) + 1 / omega_bar_l(t,l)));
    }
    for (int k = 0; k < K; k++) {
      // E[log lambda_k] component of E[log p(y)]
      elbo += 0.5 * (T - t) * (digamma_u_bar_k[t] - std::log(v_bar_k(t,k))) * pi_bar_k(t,k);
      // E[log p_j] - E[log q_j]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_k(t,k) > 1e-20) elbo += pi_bar_k(t,k) * (log_pi_k(t,k) - log_pi_bar_k(t,k));
      elbo += pi_bar_k(t,k) * (u_k * (log_v_k - std::log(v_bar_k(t,k))) + lgamma_u_bar_k[t] - lgamma_u_k);
      elbo += pi_bar_k(t,k) * ((u_k - u_bar_k[t]) * digamma_u_bar_k[t] + u_bar_k[t] * (1 - v_k / v_bar_k(t,k)));
    }
  }
  return elbo;
}

// [[Rcpp::export]]
NumericVector lambda_bar_fn(NumericVector u, NumericVector v, 
                            NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0, rev_sum = 1.0;
  NumericVector lambda_bar (T, 0.0);
  for (int t = 0; t < T; t++) {
    rev_sum -= prob[t];
    fwd_sum += (u[t] / v[t]) * prob[t];
    lambda_bar[t] = fwd_sum + rev_sum;
  }
  return lambda_bar;
}

// [[Rcpp::export]]
NumericVector mu_bar_fn(NumericVector b, NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu_bar (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += b[t] * prob[t];
    mu_bar[t] = fwd_sum;
  }
  return mu_bar;
}

// [[Rcpp::export]]
NumericVector mu2_bar_fn(NumericVector b, NumericVector omega, 
                         NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu2_bar (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += (b[t] * b[t] + 1 / omega[t]) * prob[t];
    mu2_bar[t] = fwd_sum;
  }
  return mu2_bar;
}

// [[Rcpp::export]]
NumericVector mu_lambda_fn(NumericVector b, NumericVector u, NumericVector v, 
                           NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu_lambda (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += b[t] * (u[t] / v[t]) * prob[t];
    mu_lambda[t] = fwd_sum;
  }
  return mu_lambda;
}

// [[Rcpp::export]]
NumericVector mu2_lambda_fn(NumericVector b, NumericVector omega, NumericVector u,
                            NumericVector v, NumericVector prob) {
  int T = prob.length();
  double fwd_sum = 0.0;
  NumericVector mu2_lambda (T, 0.0);
  for (int t = 0; t < T; t++) {
    fwd_sum += (b[t] * b[t] * (u[t] / v[t]) + 1 / omega[t]) * prob[t];
    mu2_lambda[t] = fwd_sum;
  }
  return mu2_lambda;
}

// [[Rcpp::export]]
List mean_scp(NumericVector y, NumericVector lambda, double omega, 
              NumericVector log_pi) {
  
  int T = y.length();
  NumericVector b_bar (T, 0.0), pi_bar (T, 0.0), log_pi_bar (T, 0.0), omega_bar (T, 0.0);
  double tot = 0, ls = 0, ys = 0, log_pi_max = R_NegInf;
  
  for (int t = T - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    omega_bar[t] = omega + ls;
    b_bar[t] = ys / omega_bar[t];
    log_pi_bar[t] = log_pi[t] - 0.5 * (std::log(omega_bar[t]) - omega_bar[t] * b_bar[t] * b_bar[t]);
    log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
  }
  
  for(int t = 0; t < T; t++) {
    log_pi_bar[t] -= log_pi_max;
    if (log_pi[t] > R_NegInf) pi_bar[t] = std::exp(log_pi_bar[t]);
    else pi_bar[t] = 0.0;
    tot += pi_bar[t];
  }
  
  return List::create(_["b_bar"] = b_bar, 
                      _["omega_bar"] = omega_bar, 
                      _["pi_bar"] = pi_bar / tot, 
                      _["log_pi_bar"] = log_pi_bar);
}

// [[Rcpp::export]]
List var_scp(NumericVector y, NumericVector omega, 
             NumericVector u_bar, NumericVector lgamma_u_bar,
             NumericVector v, NumericVector log_pi) {
  
  int T = y.length();
  NumericVector pi_bar (T, 0.0), log_pi_bar (T, 0.0), v_bar (T, 0.0);
  double tot = 0, fs = 0, rs = 0, log_pi_max = R_NegInf;
  
  for (int t = T - 1; t >= 0; t--) {
    rs += y[t] * y[t] * omega[t];
    v_bar[t] = v[t] + 0.5 * rs;
    log_pi_bar[t] += log_pi[t] + lgamma_u_bar[t] - u_bar[t] * std::log(v_bar[t]);
    log_pi_bar[T-t-1] -= 0.5 * fs;
    if (t <= T-t-1) {
      log_pi_max = std::max(log_pi_max, log_pi_bar[T-t-1]);
      log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
    }
    fs += y[T-t-1] * y[T-t-1] * omega[T-t-1];
  }
  
  for(int t = 0; t < T; t++) {
    log_pi_bar[t] -= log_pi_max;
    if (log_pi[t] > R_NegInf) pi_bar[t] = std::exp(log_pi_bar[t]);
    else pi_bar[t] = 0.0;
    tot += pi_bar[t];
  }
  
  return List::create(_["v_bar"] = v_bar, 
                      _["pi_bar"] = pi_bar / tot, 
                      _["log_pi_bar"] = log_pi_bar);
}

// [[Rcpp::export]]
List meanvar_scp(NumericVector y, NumericVector lambda, 
                 double omega, NumericVector u_bar, NumericVector lgamma_u_bar, NumericVector v, 
                 NumericVector log_pi) {
  
  int T = y.length();
  NumericVector b_bar (T, 0.0), pi_bar (T, 0.0), log_pi_bar (T, 0.0), v_bar (T, 0.0), omega_bar (T, 0.0);
  double tot = 0, fs = 0, rs = 0, ls = 0, ys = 0, log_pi_max = R_NegInf;
  
  for (int t = T - 1; t >= 0; t--) {
    ys += y[t] * lambda[t];
    ls += lambda[t];
    rs += y[t] * y[t] * lambda[t];
    omega_bar[t] = omega + ls;
    b_bar[t] = ys / omega_bar[t];
    v_bar[t] = v[t] + 0.5 * (rs - b_bar[t] * b_bar[t] * omega_bar[t]);
    log_pi_bar[t] += log_pi[t] + lgamma_u_bar[t] - u_bar[t] * std::log(v_bar[t]) - 0.5 * std::log(omega_bar[t]);
    log_pi_bar[T-t-1] -= 0.5 * fs;
    if (t <= T-t-1) {
      log_pi_max = std::max(log_pi_max, log_pi_bar[T-t-1]);
      log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
    }
    fs += y[T-t-1] * y[T-t-1] * lambda[T-t-1];
  }
  
  for(int t = 0; t < T; t++) {
    log_pi_bar[t] -= log_pi_max;
    if (log_pi[t] > R_NegInf) pi_bar[t] = std::exp(log_pi_bar[t]);
    else pi_bar[t] = 0.0;
    tot += pi_bar[t];
  }
  
  return List::create(_["b_bar"] = b_bar, 
                      _["omega_bar"] = omega_bar, 
                      _["v_bar"] = v_bar,
                      _["pi_bar"] = pi_bar / tot, 
                      _["log_pi_bar"] = log_pi_bar);
}

// [[Rcpp::export]]
List mich_cpp(NumericVector y, 
              int J, int L, int K, double mu_0, double lambda_0,
              bool fit_intercept, bool fit_scale, 
              bool refit, double max_iter, bool verbose, double tol, 
              double omega_j, double u_j, double v_j, NumericMatrix log_pi_j,
              NumericMatrix pi_bar_j, NumericMatrix log_pi_bar_j, 
              NumericMatrix b_bar_j, NumericMatrix omega_bar_j, 
              NumericVector u_bar_j, NumericMatrix v_bar_j,
              NumericVector lgamma_u_bar_j, NumericVector digamma_u_bar_j,
              double omega_l, NumericMatrix log_pi_l,
              NumericMatrix pi_bar_l, NumericMatrix log_pi_bar_l,
              NumericMatrix b_bar_l, NumericMatrix omega_bar_l, 
              double u_k, double v_k, NumericMatrix log_pi_k,
              NumericMatrix pi_bar_k, NumericMatrix log_pi_bar_k,
              NumericVector u_bar_k, NumericMatrix v_bar_k,
              NumericVector lgamma_u_bar_k, NumericVector digamma_u_bar_k) {
  
  int T = y.length();

  // log parameters
  double log_omega_j = std::log(omega_j); 
  double log_v_j = std::log(v_j); double log_u_j = std::log(u_j);
  double lgamma_u_j = std::lgamma(u_j);
  double log_omega_l = std::log(omega_l);
  double log_u_k = std::log(u_k); double log_v_k = std::log(v_k);
  double lgamma_u_k = std::lgamma(u_k);

  // initialize residual 
  NumericVector r_tilde = clone(y);
  
  // initialize precision
  NumericVector lambda_bar (T, lambda_0);
  
  // initialize correction term
  NumericVector delta (T, 0.0);
  
  // initialize J expected mean-variance parameters
  NumericMatrix mu_lambda_j (T, J), mu2_lambda_j (T, J), lambda_bar_j (T, J);
  double mu_bar_jt, mu_bar_jT;
  if (refit) {
    for (int j = 0; j < J; j++) {
      mu_lambda_j(_,j) = mu_lambda_fn(b_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      mu2_lambda_j(_,j) = mu2_lambda_fn(b_bar_j(_,j), omega_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      lambda_bar_j(_,j) = lambda_bar_fn(u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
    }
  } else {
    lambda_bar_j.fill(1.0);
  }

  // initialize L expected mean parameters
  NumericMatrix mu_bar_l (T, L), mu2_bar_l (T, L);
  if (refit) {
    for (int l = 0; l < L; l++) {
      mu_bar_l(_,l) = mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l));
      mu2_bar_l(_,l) = mu2_bar_fn(b_bar_l(_,l), omega_bar_l(_,l), pi_bar_l(_,l));
    }
  }

  // initialize K expected variance parameters
  NumericMatrix lambda_bar_k (T, K);
  if (refit) {
    for (int k = 0; k < K; k++) {
      lambda_bar_k(_,k) = lambda_bar_fn(u_bar_k, v_bar_k(_,k), pi_bar_k(_,k));
    }
  } else {
    lambda_bar_k.fill(1.0);
  }

  // initialize combined residual, variance, and correction term
  for (int t = 0; t < T; t++) {
    r_tilde[t] -= mu_0;
    
    if (refit) {
      for (int j = 0; j < J; j++) {
        mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
        r_tilde[t] -= mu_bar_jt;
        lambda_bar[t] *= lambda_bar_j(t,j);
        delta[t] -= (mu_lambda_j(t,j) * mu_bar_jt - mu2_lambda_j(t,j)) / lambda_bar_j(t,j);
      }
      
      for (int l = 0; l < L; l++) {
        r_tilde[t] -= mu_bar_l(t,l);
        delta[t] -= mu_bar_l(t,l) * mu_bar_l(t,l) - mu2_bar_l(t,l);
      }
      
      for (int k = 0; k < K; k++) {
        lambda_bar[t] *= lambda_bar_k(t,k);
      }
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
  double v_tilde_sum, log_pi_tilde_sum, log_pi_tilde_max;

  // vb algorithm
  int iter = 0;
  
  while (iter < max_iter) {
    
    // updating q(b_l, gamma_l)
    for (int l = 0; l < L; l++) {
      for (int t = T - 1; t >= 0; t--) {
        // add back l^th component of residual
        r_tilde[t] += mu_bar_l(t,l);
        // subtract l^th component of variance correction
        delta[t] += mu_bar_l(t,l) * mu_bar_l(t,l) - mu2_bar_l(t,l);
      }
      
      // fit single mean_scp model on residuals
      mean_scp_fit = mean_scp(r_tilde, lambda_bar, omega_l, log_pi_l(_,l));
      
      // store updated posterior parameters
      pi_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["pi_bar"]);
      log_pi_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["log_pi_bar"]);
      b_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["b_bar"]);
      omega_bar_l(_,l) =  as<NumericVector>(mean_scp_fit["omega_bar"]);
      
      // update l^th component of mean
      mu_bar_l(_,l) = mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l));
      mu2_bar_l(_,l) = mu2_bar_fn(b_bar_l(_,l), omega_bar_l(_,l), pi_bar_l(_,l));
      
      for (int t = 0; t < T; t++) {
        // subtract out l^th component of residual
        r_tilde[t] -= mu_bar_l(t,l);
        // add back l^th component of variance correction
        delta[t] -= mu_bar_l(t,l) * mu_bar_l(t,l) - mu2_bar_l(t,l);
      }
    }
    
    // updating q(s_k, gamma_k)
    for (int k = 0; k < K; k++) {
      v_tilde_sum = 0.0;
      log_pi_tilde_sum = log_pi_k(0,k);
      log_pi_tilde_max = R_NegInf;
      
      for (int t = T - 1; t >= 0; t--) {
        // divide out k^th component of precision
        if (t > T-t-1) {
          lambda_bar[t] /= lambda_bar_k(t,k);
          lambda_bar[T-t-1] /= lambda_bar_k(T-t-1,k);
        } else if (t == T-t-1){
          lambda_bar[t] /= lambda_bar_k(t,k);
        }
        
        // calculate corrected priors
        v_tilde_sum += 0.5 * lambda_bar[t] * delta[t];
        v_tilde[t] = v_k + v_tilde_sum;
        log_pi_tilde[T-t-1] = log_pi_k(T-t-1,k) + log_pi_tilde_sum;
        log_pi_tilde_sum += -0.5 * lambda_bar[T-t-1] * delta[T-t-1];
        log_pi_tilde_max = std::max(log_pi_tilde_max, log_pi_tilde[T-t-1]);
      }
      
      for (int t = 1; t < T; t++) {
        log_pi_tilde[t] -= log_pi_tilde_max;
      }
      
      // fit single var_scp model on residuals
      var_scp_fit = var_scp(r_tilde, lambda_bar, u_bar_k, lgamma_u_bar_k, 
                            v_tilde, log_pi_tilde);
      
      // store updated posterior parameters
      pi_bar_k(_,k) = as<NumericVector>(var_scp_fit["pi_bar"]);
      log_pi_bar_k(_,k) = as<NumericVector>(var_scp_fit["log_pi_bar"]);
      v_bar_k(_,k) = as<NumericVector>(var_scp_fit["v_bar"]);
      
      // updated k^th component precision
      lambda_bar_k(_,k) = lambda_bar_fn(u_bar_k, v_bar_k(_,k), pi_bar_k(_,k));
      
      for (int t = 0; t < T; t++) {
        // multiply back k^th component or precision
        lambda_bar[t] *= lambda_bar_k(t,k);
      }
    }

    // updating q(b_j, s_j, gamma_j)
    for (int j = 0; j < J; j++) {
      v_tilde_sum = 0.0;
      log_pi_tilde_sum = log_pi_j(0,j);
      log_pi_tilde_max = R_NegInf;
      
      // for (int t = T-1; t >= 0; t--) {
      //   // add back j^th component of residual
      //   mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
      //   r_tilde[t] += mu_bar_jt;
      //   // divide out j^th component of precision
      //   lambda_bar[t] /= lambda_bar_j(t,j);
      //   // subtract j^th component of variance correction
      //   delta[t] += mu_bar_jt * mu_bar_jt - mu2_lambda_j(t,j) / lambda_bar_j(t,j);
      // 
      //   // calculate corrected priors
      //   v_tilde_sum += 0.5 * lambda_bar[t] * delta[t];
      //   v_tilde[t] = v_j + v_tilde_sum;
      // }
      
      // for (int t = 1; t < T; t++) {
      //   log_pi_tilde_sum += -0.5 * lambda_bar[t-1] * delta[t-1];
      //   log_pi_tilde[t] = log_pi_j(t,j) + log_pi_tilde_sum;
      //   log_pi_tilde_max = std::max(log_pi_tilde_max, log_pi_tilde[t]);
      // }
      
      // for (int t = 1; t < T; t++) {
      //   log_pi_tilde[t] -= log_pi_tilde_max;
      // }
      
      for (int t = T-1; t >= 0; t--) {
        if (t > T-t-1) {
          mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
          mu_bar_jT = mu_lambda_j(T-t-1,j) / lambda_bar_j(T-t-1,j);
          // add back j^th component of residual
          r_tilde[t] += mu_bar_jt;
          r_tilde[T-t-1] += mu_bar_jT;
          // divide out j^th component of precision
          lambda_bar[t] /= lambda_bar_j(t,j);
          lambda_bar[T-t-1] /= lambda_bar_j(T-t-1,j);
          // subtract j^th component of variance correction
          delta[t] +=  mu_bar_jt * mu_bar_jt - mu2_lambda_j(t,j) / lambda_bar_j(t,j);
          delta[T-t-1] += mu_bar_jT * mu_bar_jT - mu2_lambda_j(T-t-1,j) / lambda_bar_j(T-t-1,j);
        } else if (t == T-t-1) {
          mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
          // add back j^th component of residual
          r_tilde[t] += mu_bar_jt;
          // divide out j^th component of precision
          lambda_bar[t] /= lambda_bar_j(t,j);
          // subtract j^th component of variance correction
          delta[t] += mu_bar_jt * mu_bar_jt - mu2_lambda_j(t,j) / lambda_bar_j(t,j);
        }

        // calculate corrected priors
        v_tilde_sum += 0.5 * lambda_bar[t] * delta[t];
        v_tilde[t] = v_j + v_tilde_sum;
        log_pi_tilde[T-t-1] = log_pi_j(T-t-1,j) + log_pi_tilde_sum;
        log_pi_tilde_sum += -0.5 * lambda_bar[T-t-1] * delta[T-t-1];
        log_pi_tilde_max = std::max(log_pi_tilde_max, log_pi_tilde[T-t-1]);
      }
      
      for (int t = 1; t < T; t++) {
        log_pi_tilde[t] -= log_pi_tilde_max;
      }

      // fit single meanvar_scp model on residuals
      meanvar_scp_fit = meanvar_scp(r_tilde, lambda_bar, omega_j, u_bar_j, 
                                    lgamma_u_bar_j, v_tilde, log_pi_tilde);

      // store updated posterior parameters
      pi_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["pi_bar"]);
      log_pi_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["log_pi_bar"]);
      b_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["b_bar"]);
      omega_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["omega_bar"]);
      v_bar_j(_,j) =  as<NumericVector>(meanvar_scp_fit["v_bar"]);

      // update j^th component of mean and precision
      mu_lambda_j(_,j) = mu_lambda_fn(b_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      mu2_lambda_j(_,j) = mu2_lambda_fn(b_bar_j(_,j), omega_bar_j(_,j), u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));
      lambda_bar_j(_,j) = lambda_bar_fn(u_bar_j, v_bar_j(_,j), pi_bar_j(_,j));

      for (int t = 0; t < T; t++) {
        mu_bar_jt = mu_lambda_j(t,j) / lambda_bar_j(t,j);
        // subtract out j^th component of residual
        r_tilde[t] -= mu_bar_jt;
        // multiply back j^th component or precision
        lambda_bar[t] *= lambda_bar_j(t,j);
        // add back j^th component of variance correction
        delta[t] -= mu_bar_jt * mu_bar_jt - mu2_lambda_j(t,j) / lambda_bar_j(t,j);
      }
    }

    // updating mu_0 and lambda_0
    if (fit_intercept) {
      r_tilde += Rcpp::rep(mu_0, T);
      mu_0 = Rcpp::sum(lambda_bar * r_tilde) / Rcpp::sum(lambda_bar);
      r_tilde += Rcpp::rep(-mu_0, T);
    }

    if (fit_scale) {
      for (int t = 0; t < T; t++) {
        lambda_bar[t] /= lambda_0;
      }

      lambda_0 = T / Rcpp::sum(lambda_bar * (Rcpp::pow(r_tilde, 2) + delta));

      for (int t = 0; t < T; t++) {
        lambda_bar[t] *= lambda_0;
      }
    }

    iter++;
    elbo.push_back(0.0);
    
    elbo[iter] = elbo_fn(T, mu_0, lambda_0,
                         r_tilde, lambda_bar, delta,
                         b_bar_j, omega_bar_j, u_bar_j, v_bar_j, pi_bar_j, log_pi_bar_j,
                         lgamma_u_bar_j, digamma_u_bar_j,
                         omega_j, u_j, v_j, log_omega_j, log_u_j, lgamma_u_j, log_v_j, log_pi_j,
                         b_bar_l, omega_bar_l, pi_bar_l, log_pi_bar_l,
                         omega_l, log_omega_l, log_pi_l, 
                         u_bar_k, v_bar_k, pi_bar_k,
                         lgamma_u_bar_k, digamma_u_bar_k, log_pi_bar_k,
                         u_k, v_k, log_u_k, lgamma_u_k, log_v_k, log_pi_k);
    if (verbose & (iter % 5000 == 0)) Rcout << "Iteration: " << iter << "; elbo: " << elbo[iter] << ";\n";
    if (std::abs((elbo[iter] - elbo[iter - 1]) / elbo[iter - 1]) < tol) break;
  }
  
  // construct mean signal
  NumericVector mu (T, mu_0);
  for (int j = 0; j < J; j++) {
    mu += mu_bar_fn(b_bar_j(_,j), pi_bar_j(_,j));
  }
  for (int l = 0; l < L; l++) {
    mu += mu_bar_fn(b_bar_l(_,l), pi_bar_l(_,l));
  }

  // creating lists of posterior parameters
  List J_model, L_model, K_model;
  
  List result = List::create(_["y"] = y,
                             _["residual"] = r_tilde, 
                             _["mu"] = mu,
                             _["lambda"] = lambda_bar,
                             _["delta"] = delta,
                             _["converged"] = (max_iter > iter),
                             _["elbo"] = elbo,
                             _["mu_0"] = mu_0, 
                             _["lambda_0"] = lambda_0,
                             _["J"] = J, 
                             _["L"] = L, 
                             _["K"] = K);

  if (J > 0) {
    J_model =  List::create(_["pi_bar"] = pi_bar_j, 
                            _["b_bar"] = b_bar_j,
                            _["omega_bar"] = omega_bar_j, 
                            _["v_bar"] = v_bar_j,
                            _["u_bar"] = u_bar_j,
                            _["mu_lambda_bar"] = mu_lambda_j,
                            _["mu2_lambda_bar"] = mu2_lambda_j,
                            _["lambda_bar"] = lambda_bar_j);
    result["meanvar_model"] = J_model;
  }
  if (L > 0) {
    L_model =  List::create(_["pi_bar"] = pi_bar_l, 
                            _["b_bar"] = b_bar_l,
                            _["omega_bar"] = omega_bar_l, 
                            _["mu_bar"] = mu_bar_l);
    result["mean_model"] = L_model;
  }
  if (K > 0) {
    K_model =  List::create(_["pi_bar"] = pi_bar_k, 
                            _["v_bar"] = v_bar_k,
                            _["u_bar"] = u_bar_k,
                            _["lambda_bar"] = lambda_bar_k);
    result["var_model"] = K_model;
  }

  return result;
}