#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
double multi_elbo_fn(int T, int B_l, int B_r, NumericVector mu_0, 
                     NumericMatrix y_0, NumericMatrix r_bar, 
                     double omega_l, NumericMatrix log_pi_l,
                     NumericVector mu_var, List post_params) {
  
  int d = r_bar.ncol();
  int L = post_params.length();
  
  double elbo = -0.5 * d * (T + B_l + B_r) * std::log(2 * M_PI);
  List post_params_l;
  NumericMatrix b_bar_l (T, d);
  NumericVector pi_bar_l (T), omega_bar_l (T);
  
  for (int t = 0; t < B_l; t++) {
    for (int i = 0; i < d; i++) {
      elbo -= 0.5 * std::pow(y_0(t,i) - mu_0[i], 2);
    }
  }
  
  for (int t = 0; t < T + B_r; t++) {
    elbo += mu_var[t]; 
    for (int i = 0; i < d; i++) {
      elbo -= 0.5 * std::pow(r_bar(t,i), 2);
    }
    
    if (t < T) {
      for (int l = 0; l < L; l++) {
        post_params_l = post_params[l];
        pi_bar_l = as<NumericVector>(post_params_l["pi_bar"]);
        b_bar_l = as<NumericMatrix>(post_params_l["b_bar"]);
        omega_bar_l = as<NumericVector>(post_params_l["omega_bar"]);
        
        // E[log p_l] - E[log q_l]
        // adding log odds with numerical stability for small pi_bar_l values
        if (pi_bar_l[t] > 1e-100) elbo += pi_bar_l[t] * (log_pi_l(t,l) - std::log(pi_bar_l[t]));
        elbo += 0.5 * d * pi_bar_l[t]  * (std::log(omega_l) - std::log(omega_bar_l[t]));
        
        for (int i = 0; i < d; i++) {
          elbo += 0.5 * pi_bar_l[t] * (1 - omega_l * (std::pow(b_bar_l(t,i), 2) + 1 / omega_bar_l[t]));
        }
      }  
    }
  }
  return elbo;
}

// [[Rcpp::export]]
List multi_mean_scp(NumericMatrix y, double lambda_0, 
                    NumericVector log_pi, int B_r) {
  int T = y.nrow() - B_r;
  int d = y.ncol();
  NumericVector pi_bar (T), ys (d), omega_bar (T);
  NumericMatrix b_bar (T, d);
  double tot = 0, cntr = 0, log_pi_max = R_NegInf;
  
  for (int t = T + B_r - 1; t >= 0; t--) {
    cntr += 1;
    if (t < T) {
      omega_bar[t] = cntr + lambda_0;
      pi_bar[t] = log_pi[t];
    }
    for (int i = 0; i < d; i++) {
      ys[i] += y(t, i);
      if (t < T) {
        b_bar(t, i) = ys[i] / omega_bar[t];
        pi_bar[t] += 0.5 *(std::pow(ys[i], 2.0) / omega_bar[t] - std::log(omega_bar[t]));
      }
    }
    if (t < T && pi_bar[t] > log_pi_max) log_pi_max = pi_bar[t];
  }
  
  for(int t = 0; t < T; t++) {
    pi_bar[t] = std::exp(pi_bar[t] - log_pi_max);
    tot += pi_bar[t];
  }
  return List::create(_["b_bar"] = b_bar, _["omega_bar"] = omega_bar,
                      _["pi_bar"] = pi_bar / tot);
}

// [[Rcpp::export]]
NumericMatrix multi_mu_bar_fn(NumericMatrix b, NumericVector prob, int B_r) {
  int T = prob.length();
  int d = b.ncol();
  NumericVector fwd_sum (d);
  NumericMatrix mu_bar (T + B_r, d);
  for (int t = 0; t < T + B_r; t++) {
    for (int i = 0; i < d; i++) {
      if (t < T) {
        fwd_sum[i] += b(t, i) * prob[t];
        mu_bar(t, i) = fwd_sum[i];
      } else {
        mu_bar(t, i) = mu_bar(T-1, i);
      }
    }
  }
  return mu_bar;
}

// [[Rcpp::export]]
NumericMatrix multi_mu2_bar_fn(NumericMatrix b, NumericVector omega, 
                               NumericVector prob, int B_r) {
  int T = prob.length();
  int d = b.ncol();
  NumericVector fwd_sum (d);
  NumericMatrix mu2_bar (T + B_r, d);
  for (int t = 0; t < T + B_r; t++) {
    for (int i = 0; i < d; i++) {
      if (t < T) {
        fwd_sum[i] += (std::pow(b(t,i), 2) + 1 / omega[t]) * prob[t];
        mu2_bar(t, i) = fwd_sum[i];
      } else {
        mu2_bar(t, i) = mu2_bar(T-1, i);
      }
    }
  }
  return mu2_bar;
}

// [[Rcpp::export]]
List multi_mich_cpp(NumericMatrix y_0, NumericMatrix y, NumericVector mu_0,
                    int B_l, int B_r, bool fit_intercept, bool refit,
                    double max_iter, double tol, bool verbose,  
                    double omega_l, NumericMatrix log_pi_l, 
                    List post_params) {
  
  int T = y.nrow() - B_r;
  int d = y.ncol();
  int L = post_params.length();
  
  NumericVector y_0_sum = Rcpp::colSums(y_0);
  
  // initialize residual
  NumericMatrix r_bar = clone(y);
  for (int t = 0; t < T + B_r; t++) {
    for (int i = 0; i < d; i++) {
      r_bar(t,i) -= mu_0[i];
    }
  }
  
  // initialize ELBO
  std::vector<double> elbo;
  elbo.reserve(std::min(max_iter, 1e6));
  elbo.push_back(R_NegInf);
  
  List post_params_l, mu_bar, mu2_bar, mean_scp_fit;
  NumericVector mu_var (T + B_r), mu_0_new (d);
  NumericMatrix mu_bar_l (T + B_r, d), mu2_bar_l (T + B_r, d), mu_var_l (T + B_r, L);
  
  for (int l = 0; l < L; l++) {
    post_params_l = post_params[l];
    if (refit) {
      mu_bar_l = multi_mu_bar_fn(post_params_l["b_bar"], post_params_l["pi_bar"], B_r);
      mu2_bar_l = multi_mu2_bar_fn(post_params_l["b_bar"], post_params_l["omega_bar"], post_params_l["pi_bar"], B_r);
    }
    mu_bar.push_back(mu_bar_l);
    mu2_bar.push_back(mu2_bar_l);
    for (int t = 0; t < T + B_r; t++) {
      for (int i = 0; i < d; i++) { 
        if (refit) {
          r_bar(t,i) -= mu_bar_l(t,i);
          mu_var_l(t,l) += mu2_bar_l(t,i) - std::pow(mu_bar_l(t,i), 2);
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
      
      for (int t = 0; t < T + B_r; t++) {
        for (int i = 0; i < d; i++) {
          // add back l^th component of residual
          r_bar(t,i) += mu_bar_l(t,i);
          // add back l^th component of variance of mu
          mu_var_l(t,l) -= mu2_bar_l(t,i) - std::pow(mu_bar_l(t,i), 2);
        }
      }
      
      // fit single multi_mean_scp model on residuals
      mean_scp_fit = multi_mean_scp(r_bar, omega_l, log_pi_l, B_r);
      
      // store updated prior parameters
      post_params_l["b_bar"] = as<NumericMatrix>(mean_scp_fit["b_bar"]);
      // post_params_l["omega_bar"] = as<NumericVector>(mean_scp_fit["omega_bar"]);
      post_params_l["pi_bar"] = as<NumericVector>(mean_scp_fit["pi_bar"]);
      post_params[l] = post_params_l;
      
      // update l^th component of mean
      mu_bar_l = multi_mu_bar_fn(post_params_l["b_bar"], post_params_l["pi_bar"], B_r);
      mu_bar[l] = mu_bar_l;
      mu2_bar_l = multi_mu2_bar_fn(post_params_l["b_bar"], post_params_l["omega_bar"], post_params_l["pi_bar"], B_r);
      mu2_bar[l] = mu2_bar_l;
      
      for (int t = T + B_r -1; t >= 0; t--) {
        for (int i = 0; i < d; i++){
          // subtract out l^th component of residual
          r_bar(t, i) -= mu_bar_l(t,i);
          // subtract out l^th component of variance of mu
          mu_var_l(t, l) += mu2_bar_l(t,i) - std::pow(mu_bar_l(t,i), 2);
        }
      }
    }
    
    // update mu_0
    if (fit_intercept) {
      mu_0_new = clone(y_0_sum);
      for (int i = 0; i < d; i++) {
        for (int t = 0; t < T + B_r; t++) {
          r_bar(t, i) += mu_0[i];
          mu_0_new[i] += r_bar(t, i);
        }
        
        mu_0[i] = mu_0_new[i] / (T + B_l + B_r);
        for (int t = 0; t < T + B_r; t++) {
          r_bar(t, i) -= mu_0[i];
        }
      }
    }
    
    // calculate var(mu_t)
    mu_var = Rcpp::rowSums(mu_var_l);
    
    iter++;
    elbo.push_back(0.0);
    elbo[iter] = multi_elbo_fn(T, B_l, B_r, mu_0, y_0, r_bar, omega_l, log_pi_l, mu_var, post_params); 
    if (verbose & (iter % 1000 == 0)) Rcout << "Iteration: " << iter << " elbo: " << elbo[iter] << "\n"; 
    if (elbo[iter] - elbo[iter - 1] < tol) break;
  }
  
  // construct mean signal
  NumericMatrix mu (T + B_l + B_r, d), resid (T + B_l + B_r, d);
  for (int t = 0; t < T + B_l + B_r; t++) {
    if (t < B_l) {
      mu(t,_) = mu_0;
      resid(t,_) = y_0(t,_) - mu_0;
    } else {
      mu(t,_) = y(t - B_l,_) - r_bar(t - B_l,_);
      resid(t,_) = r_bar(t - B_l,_);
    }
  }
  
  // creating lists of posterior parameters
  List result = List::create(_["residual"] = resid, _["mu"] = mu, _["mu_0"] = mu_0, 
                             _["L"] = L, _["post_params"] = post_params,
                             _["elbo"] = elbo, _["converged"] = (max_iter > iter));
  return result;
}