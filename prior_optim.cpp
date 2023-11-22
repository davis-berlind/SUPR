#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector prior_j(int T, int B_r, int N, double tau, double u, double v) {
  
  double half = T / 2.0;
  NumericVector moment (T, 0.0), y (T + B_r);
  NumericVector pi_bar (T), mu_bar (T), tau_bar (T), v_bar (T), u_bar (T);
  double tot, ys, ls, fs, rs, log_pi_max = R_NegInf;
  
  for (int i = 0; i < N; i++) {
    ys = 0.0, ls = 0.0, fs = 0.0, rs = 0.0, tot = 0.0;
    pi_bar.fill(0.0);
    y = Rcpp::rnorm(T + B_r);
    
    for (int t = T + B_r - 1; t >= 0; t--) {
      ys += y[t];
      ls += 1;
      rs += std::pow(y[t], 2.0);
      if (t < T) {
        tau_bar[t] = tau + ls;
        mu_bar[t] = ys / tau_bar[t];
        v_bar[t] = v + 0.5 * (rs - std::pow(mu_bar[t], 2.0) * tau_bar[t]);
        u_bar[t] = u + 0.5 * (T + B_r - t);
        pi_bar[t] += std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]) - 0.5 * std::log(tau_bar[t]);
        pi_bar[T-t-1] -= 0.5 * fs;
        if (t < half) {
          if (pi_bar[T-t-1] > log_pi_max) log_pi_max = pi_bar[T-t-1];
          if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
        }
        fs += std::pow(y[T-t-1], 2.0);
      }
    }
    
    for(int t = 0; t < T; t++) {
      pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
      tot += pi_bar[t];
    }
    
    log_pi_max = R_NegInf;
    
    for(int t = 0; t < T; t++) {
      moment[t] += pi_bar[t] / tot;
    }
  }
  
  moment = moment / N;
  moment = 1.0 / moment;
  moment = moment / Rcpp::sum(moment);
  return moment;
}

// [[Rcpp::export]]
NumericVector prior_l(int T, int B_r, int N, double tau) {
  
  NumericVector moment (T, 0.0), y (T + B_r);
  NumericVector pi_bar (T, 0.0), mu_bar (T, 0.0), tau_bar (T, 0.0);
  double tot, ls, ys, log_pi_max = R_NegInf;
  
  for (int i = 0; i < N; i++) {
    ys = 0.0, ls = 0.0, tot = 0.0;
    y = Rcpp::rnorm(T + B_r);
    
    for (int t = T + B_r - 1; t >= 0; t--) {
      ys += y[t];
      ls += 1;
      if (t < T) {
        tau_bar[t] = tau + ls;
        mu_bar[t] = ys / tau_bar[t];
        pi_bar[t] = -0.5 * (std::log(tau_bar[t]) - tau_bar[t] * std::pow(mu_bar[t], 2.0));
        if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
      }
    }
    
    for(int t = 0; t < T; t++) {
      pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
      tot += pi_bar[t];
    }
    
    log_pi_max = R_NegInf;
    
    for(int t = 0; t < T; t++) {
      moment[t] += pi_bar[t] / tot;
    }
  }
  
  moment = moment / N;
  moment = 1.0 / moment;
  moment = moment / Rcpp::sum(moment);
  return moment;
}

// [[Rcpp::export]]
NumericVector prior_k(int T, int B_r, int N, double u, double v) {
  
  double half = T / 2.0;
  NumericVector moment (T, 0.0), y (T + B_r);
  NumericVector pi_bar (T), v_bar (T), u_bar (T);
  double tot, fs, rs, log_pi_max = R_NegInf;
  
  for (int i = 0; i < N; i++) {
    fs = 0.0, rs = 0.0, tot = 0.0;
    pi_bar.fill(0.0);
    y = Rcpp::rnorm(T + B_r);
    
    for (int t = T + B_r - 1; t >= 0; t--) {
      rs += std::pow(y[t], 2.0);
      if (t < T) {
        v_bar[t] = v + 0.5 * rs;
        u_bar[t] = u + 0.5 * (T + B_r - t);
        pi_bar[t] += std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]);
        pi_bar[T-t-1] -= 0.5 * fs;
        if (t < half) {
          if (pi_bar[T-t-1] > log_pi_max) log_pi_max = pi_bar[T-t-1];
          if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
        }
        fs += std::pow(y[T-t-1], 2.0);
      }
    }
    
    for(int t = 0; t < T; t++) {
      pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
      tot += pi_bar[t];
    }
    
    log_pi_max = R_NegInf;
    
    for(int t = 0; t < T; t++) {
      moment[t] += pi_bar[t] / tot;
    }
  }
  
  moment = moment / N;
  moment = 1.0 / moment;
  moment = moment / Rcpp::sum(moment);
  return moment;
}

// [[Rcpp::export]]
double objective_j(NumericVector log_pi, NumericMatrix mat, NumericMatrix y, 
                   double tau, double u, double v, int B_r) {
  
  int N = y.ncol();
  int T = y.nrow() - B_r;
  double half = T / 2.0;
  NumericVector moment (T, 0.0), pi_bar (T), mu_bar (T), tau_bar (T), v_bar (T), u_bar (T);
  double ret = 0.0, tot, ys, ls, fs, rs, log_pi_max = R_NegInf;
  
  for (int i = 0; i < N; i++) {
    ys = 0.0, ls = 0.0, fs = 0.0, rs = 0.0, tot = 0.0;
    pi_bar.fill(0.0);
    
    for (int t = T + B_r - 1; t >= 0; t--) {
      ys += y(t, i);
      ls += 1;
      rs += std::pow(y(t,i), 2.0);
      if (t < T) {
        tau_bar[t] = tau + ls;
        mu_bar[t] = ys / tau_bar[t];
        v_bar[t] = v + 0.5 * (rs - std::pow(mu_bar[t], 2.0) * tau_bar[t]);
        u_bar[t] = u + 0.5 * (T + B_r - t);
        pi_bar[t] += log_pi[t] + std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]) - 0.5 * std::log(tau_bar[t]);
        pi_bar[T-t-1] -= 0.5 * fs;
        if (t < half) {
          if (pi_bar[T-t-1] > log_pi_max) log_pi_max = pi_bar[T-t-1];
          if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
        }
        fs += std::pow(y(T-t-1,i), 2.0);
      }
    }
    
    for(int t = 0; t < T; t++) {
      pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
      tot += pi_bar[t];
    }
    
    log_pi_max = R_NegInf;
    
    for(int t = 0; t < T; t++) {
      moment[t] += (pi_bar[t] / tot - std::pow(T, -1)) / N;
    }
  }
  
  for (int i = 0; i < T; i++) {
    for (int j = i; j < T; j++) {
      if (mat(i,j) == 0.0) continue;
      if (i == j) ret += std::pow(moment[i], 2.0) * mat(i,i);
      else ret += 2 * moment[i] * mat(i,j) * moment[j];
    }
  }
  
  return ret;
}

// [[Rcpp::export]]
double objective_l(NumericVector log_pi, NumericMatrix mat, NumericMatrix y, 
                   double tau, int B_r) {
  
  int N = y.ncol();
  int T = y.nrow() - B_r;
  NumericVector moment (T, 0.0), mu_bar (T, 0.0), tau_bar (T, 0.0), pi_bar (T, 0.0);
  double ret = 0.0, tot, ls, ys, log_pi_max = R_NegInf;
  
  for (int i = 0; i < N; i++) {
    ys = 0.0;
    ls = 0.0;
    tot = 0.0;
    for (int t = T + B_r - 1; t >= 0; t--) {
      ys += y(t,i);
      ls += 1;
      if (t < T) {
        tau_bar[t] = tau + ls;
        mu_bar[t] = ys / tau_bar[t];
        pi_bar[t] = log_pi[t] - 0.5 * (std::log(tau_bar[t]) - tau_bar[t] * std::pow(mu_bar[t], 2.0));
        if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
      }
    }
    
    for(int t = 0; t < T; t++) {
      pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
      tot += pi_bar[t];
    }
    
    log_pi_max = R_NegInf;
    
    for(int t = 0; t < T; t++) {
      moment[t] += (pi_bar[t] / tot - std::pow(T, -1)) / N;
    }
  }
  
  for (int i = 0; i < T; i++) {
    for (int j = i; j < T; j++) {
      if (mat(i,j) == 0.0) continue;
      if (i == j) ret += std::pow(moment[i], 2.0) * mat(i,i);
      else ret += 2 * moment[i] * mat(i,j) * moment[j];
    }
  }
  
  return ret;
}

// [[Rcpp::export]]
double objective_k(NumericVector log_pi, NumericMatrix mat, NumericMatrix y, 
                   double u, double v, int B_r) {
  
  int N = y.ncol();
  int T = y.nrow() - B_r;
  double half = T / 2.0;
  NumericVector moment (T, 0.0), pi_bar (T), v_bar (T), u_bar (T);
  double ret = 0.0, tot, fs, rs, log_pi_max = R_NegInf;
  
  for (int i = 0; i < N; i++) {
    fs = 0.0, rs = 0.0, tot = 0.0;
    pi_bar.fill(0.0);
    
    for (int t = T + B_r - 1; t >= 0; t--) {
      rs += std::pow(y(t,i), 2.0);
      if (t < T) {
        v_bar[t] = v + 0.5 * rs;
        u_bar[t] = u + 0.5 * (T + B_r - t);
        pi_bar[t] += log_pi[t] + std::lgamma(u_bar[t]) - u_bar[t] * std::log(v_bar[t]);
        pi_bar[T-t-1] -= 0.5 * fs;
        if (t < half) {
          if (pi_bar[T-t-1] > log_pi_max) log_pi_max = pi_bar[T-t-1];
          if (pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
        }
        fs += std::pow(y(T-t-1, i), 2.0);
      }
    }
    
    for(int t = 0; t < T; t++) {
      pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
      tot += pi_bar[t];
    }
    
    log_pi_max = R_NegInf;
    
    for(int t = 0; t < T; t++) {
      moment[t] += (pi_bar[t] / tot - std::pow(T, -1)) / N;
    }
  }
  
  for (int i = 0; i < T; i++) {
    for (int j = i; j < T; j++) {
      if (mat(i,j) == 0.0) continue;
      if (i == j) ret += std::pow(moment[i], 2.0) * mat(i,i);
      else ret += 2 * moment[i] * mat(i,j) * moment[j];
    }
  }
  
  return ret;
}