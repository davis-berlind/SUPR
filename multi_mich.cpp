#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

const double log_2_pi = std::log(2 * M_PI);

// [[Rcpp::export]]
double multi_elbo_fn(NumericVector mu_0, NumericMatrix r_bar, 
                     double omega_l, double log_omega_l, NumericMatrix log_pi_l,
                     NumericVector mu_var, List post_params, 
                     NumericVector omega_bar_l, NumericVector log_omega_bar_l) {
  
  int T = r_bar.nrow();
  int d = r_bar.ncol();
  int L = post_params.length();
  
  double elbo = -0.5 * d * T * log_2_pi;
  List post_params_l;
  NumericMatrix b_bar_l (T, d);
  NumericVector pi_bar_l (T), log_pi_bar_l (T);
  
  for (int t = 0; t < T; t++) {
    elbo += mu_var[t]; 
    for (int i = 0; i < d; i++) {
      elbo -= 0.5 * r_bar(t,i) * r_bar(t,i);
    }
    
    for (int l = 0; l < L; l++) {
      post_params_l = post_params[l];
      pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
      log_pi_bar_l = as<NumericVector>(post_params_l["log_pi_bar"]);
      b_bar_l = as<NumericMatrix>(post_params_l["b_bar"]);

      // E[log p_l] - E[log q_l]
      // adding log odds with numerical stability for small pi_bar_l values
      if (pi_bar_l[t] > 1e-20) elbo += pi_bar_l[t] * (log_pi_l(t,l) - log_pi_bar_l[t]);
      elbo += 0.5 * d * pi_bar_l[t]  * (log_omega_l - log_omega_bar_l[t]);
      
      for (int i = 0; i < d; i++) {
        elbo += 0.5 * pi_bar_l[t] * (1 - omega_l * (b_bar_l(t,i) * b_bar_l(t,i) + 1 / omega_bar_l[t]));
      }
    }
  }
  return elbo;
}

// [[Rcpp::export]]
List multi_mean_scp(NumericMatrix y, NumericVector omega_bar, 
                    NumericVector log_omega_bar, NumericVector log_pi) {
  int T = y.nrow();
  int d = y.ncol();
  NumericVector pi_bar (T, 0.0), log_pi_bar (T, 0.0), ys (d, 0.0);
  NumericMatrix b_bar (T, d);
  double tot = 0, log_pi_max = R_NegInf;
  
  for (int t = T - 1; t >= 0; t--) {
    log_pi_bar[t] = log_pi[t];
    for (int i = 0; i < d; i++) {
      ys[i] += y(t, i);
      b_bar(t, i) = ys[i] / omega_bar[t];
      log_pi_bar[t] += 0.5 *(ys[i] * ys[i] / omega_bar[t] - log_omega_bar[t]);
    }
    log_pi_max = std::max(log_pi_max, log_pi_bar[t]);
  }
  
  for(int t = 0; t < T; t++) {
    log_pi_bar[t] -= log_pi_max;
    if (log_pi[t] > R_NegInf) pi_bar[t] = std::exp(log_pi_bar[t]);
    else pi_bar[t] = 0.0;
    tot += pi_bar[t];
  }
  
  return List::create(_["b_bar"] = b_bar, 
                      _["pi_bar"] = pi_bar / tot, 
                      _["log_pi_bar"] = log_pi_bar);
}

// [[Rcpp::export]]
NumericMatrix multi_mu_bar_fn(NumericMatrix b, NumericVector prob) {
  int T = prob.length();
  int d = b.ncol();
  NumericVector fwd_sum (d);
  NumericMatrix mu_bar (T, d);
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < d; i++) {
      fwd_sum[i] += b(t, i) * prob[t];
      mu_bar(t, i) = fwd_sum[i];
    }
  }
  return mu_bar;
}

// [[Rcpp::export]]
NumericMatrix multi_mu2_bar_fn(NumericMatrix b, NumericVector omega, 
                               NumericVector prob) {
  int T = prob.length();
  int d = b.ncol();
  NumericVector fwd_sum (d);
  NumericMatrix mu2_bar (T, d);
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < d; i++) {
      fwd_sum[i] += (b(t,i) * b(t,i) + 1 / omega[t]) * prob[t];
      mu2_bar(t, i) = fwd_sum[i];
    }
  }
  return mu2_bar;
}

// [[Rcpp::export]]
List multi_mich_cpp(NumericMatrix y, NumericVector mu_0,
                    bool fit_intercept, bool refit,
                    double max_iter, double tol, bool verbose,  
                    double omega_l, NumericMatrix log_pi_l, 
                    NumericVector omega_bar_l, NumericVector log_omega_bar_l,
                    List post_params) {
  
  int T = y.nrow();
  int d = y.ncol();
  int L = post_params.length();
  
  // log parameters
  double log_omega_l = std::log(omega_l);
  
  // initialize residual
  NumericMatrix r_bar = clone(y);
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < d; i++) {
      r_bar(t,i) -= mu_0[i];
    }
  }
  
  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);
  
  // initialize L expected mean parameters
  List post_params_l, mu_bar, mu2_bar, mean_scp_fit;
  NumericVector mu_var (T), mu_0_new (d);
  NumericMatrix mu_bar_l (T, d), mu2_bar_l (T, d), mu_var_l (T, L);
  double cntr = 0;
  
  for (int l = 0; l < L; l++) {
    post_params_l = post_params[l];
    if (refit) {
      mu_bar_l = multi_mu_bar_fn(post_params_l["b_bar"], post_params_l["pi_bar"]);
      mu2_bar_l = multi_mu2_bar_fn(post_params_l["b_bar"], omega_bar_l, post_params_l["pi_bar"]);
    }
    mu_bar.push_back(mu_bar_l);
    mu2_bar.push_back(mu2_bar_l);
    
    // initialize residual and variance terms
    for (int t = 0; t < T ; t++) {
      cntr += 1;
      for (int i = 0; i < d; i++) { 
        if (refit) {
          r_bar(t,i) -= mu_bar_l(t,i);
          mu_var_l(t,l) += mu2_bar_l(t,i) - mu_bar_l(t,i) * mu_bar_l(t,i);
        }
      }
    }
  }
  
  // vb algorithm
  int iter = 0;
  
  while (iter < max_iter) {
    
    // updating q(b_l, gamma_l)
    for (int l = 0; l < L; l++) {
      post_params_l = post_params[l];
      mu_bar_l = as<NumericMatrix>(mu_bar[l]);
      mu2_bar_l = as<NumericMatrix>(mu2_bar[l]);
      
      for (int t = 0; t < T; t++) {
        for (int i = 0; i < d; i++) {
          // add back l^th component of residual
          r_bar(t,i) += mu_bar_l(t,i);
          // add back l^th component of variance of mu
          mu_var_l(t,l) -= mu2_bar_l(t,i) - mu_bar_l(t,i) * mu_bar_l(t,i);
        }
      }
      
      // fit single multi_mean_scp model on residuals
      mean_scp_fit = multi_mean_scp(r_bar, omega_bar_l, log_omega_bar_l, log_pi_l(_,l));
      
      // store updated prior parameters
      post_params_l["b_bar"] = as<NumericMatrix>(mean_scp_fit["b_bar"]);
      post_params_l["pi_bar"] = as<NumericVector>(mean_scp_fit["pi_bar"]);
      post_params_l["log_pi_bar"] = as<NumericVector>(mean_scp_fit["log_pi_bar"]);
      post_params[l] = post_params_l;
      
      // update l^th component of mean
      mu_bar_l = multi_mu_bar_fn(post_params_l["b_bar"], post_params_l["pi_bar"]);
      mu_bar[l] = mu_bar_l;
      mu2_bar_l = multi_mu2_bar_fn(post_params_l["b_bar"], omega_bar_l, post_params_l["pi_bar"]);
      mu2_bar[l] = mu2_bar_l;
      
      for (int t = T - 1; t >= 0; t--) {
        for (int i = 0; i < d; i++){
          // subtract out l^th component of residual
          r_bar(t, i) -= mu_bar_l(t,i);
          // subtract out l^th component of variance of mu
          mu_var_l(t, l) += mu2_bar_l(t,i) - mu_bar_l(t,i) * mu_bar_l(t,i);
        }
      }
    }
    
    // update mu_0
    if (fit_intercept) {
      mu_0_new.fill(0.0);
      for (int i = 0; i < d; i++) {
        for (int t = 0; t < T; t++) {
          r_bar(t, i) += mu_0[i];
          mu_0_new[i] += r_bar(t, i);
        }
        
        mu_0[i] = mu_0_new[i] / T;
        for (int t = 0; t < T; t++) {
          r_bar(t, i) -= mu_0[i];
        }
      }
    }
    
    // calculate var(mu_t)
    mu_var = Rcpp::rowSums(mu_var_l);
    
    iter++;
    elbo.push_back(0.0);
    elbo[iter] = multi_elbo_fn(mu_0, r_bar, omega_l, log_omega_l, log_pi_l,
                               mu_var, post_params, omega_bar_l, log_omega_bar_l); 
    if (verbose & (iter % 5000 == 0)) Rcout << "Iteration: " << iter << " elbo: " << elbo[iter] << "\n"; 
    if (std::abs((elbo[iter] - elbo[iter - 1]) / elbo[iter - 1]) < tol) break;
  }
  
  // construct mean signal
  NumericMatrix mu (T, d);
  for (int t = 0; t < T; t++) {
    mu(t,_) = y(t,_) - r_bar(t,_);
  }
  
  // creating lists of posterior parameters
  List result = List::create(_["residual"] = r_bar, 
                             _["mu"] = mu, 
                             _["mu_0"] = mu_0, 
                             _["L"] = L, 
                             _["post_params"] = post_params,
                             _["elbo"] = elbo, 
                             _["converged"] = (max_iter > iter));
  return result;
}