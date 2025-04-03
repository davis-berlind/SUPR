#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List mich(NumericVector y, int J, int L, int K, 
          int B_l, int B_r, bool fit_intercept, bool fit_scale, 
          int max_iter, bool verbose, 
          NumericVector tau_j, NumericVector u_j, NumericVector v_j, NumericMatrix pi_j,
          NumericVector tau_l, NumericMatrix pi_l,
          NumericVector u_k, NumericVector v_k, NumericMatrix pi_k) {
  
  // initialize intercept terms
  double mu_0 = 0.0, lambda_0 = 1.0;
  
  // extract observations in start buffer
  NumericVector y_0 (B_l);
  if (B_l > 0) {
    y_0 = y[Rcpp::Range(0, B_l - 1)];
    y.erase(0,B_l);
  }
  int T = y.length() - B_r;
  double half = T / 2.0;
  
  // initialize J components
  NumericMatrix b_bar_j (T, J);
  NumericMatrix tau_bar_j (T, J);
  NumericMatrix u_bar_j (T, J);
  NumericMatrix v_bar_j (T, J);
  NumericMatrix pi_bar_j (T, J);
  NumericMatrix mu_lambda_j (T + B_r, J);
  NumericMatrix mu2_lambda_j (T + B_r, J);
  NumericMatrix lambda_bar_j (T + B_r, J);

  for (int j = 0; j < J; j++) {
    for (int t = 0; t < T + B_r; t++) {
      lambda_bar_j(t,j) = 1.0;
      if (t < T) {
        tau_bar_j(t,j) = 1.0;
        u_bar_j(t,j) = u_j[j] + (T + B_r - t) / 2.0;
        v_bar_j(t,j) = u_bar_j(t,j);
        pi_bar_j(t,j) = (double) 1 / T;
      }
    }
  }
  
  // initialize L components
  NumericMatrix b_bar_l (T, L);
  NumericMatrix tau_bar_l (T, L);
  NumericMatrix pi_bar_l (T, L);
  NumericMatrix mu_bar_l (T + B_r, L);
  NumericMatrix mu2_bar_l (T + B_r, L);

  for (int l = 0; l < L; l++) {
    for (int t = 0; t < T; t++) {
      tau_bar_l(t,l) = 1.0;
      pi_bar_l(t,l) = (double) 1 / T;
    }
  }
  
  // initialize K components
  NumericMatrix u_bar_k (T, K);
  NumericMatrix v_bar_k (T, K);
  NumericMatrix pi_bar_k (T, K);
  NumericMatrix lambda_bar_k (T + B_r, K);

  for (int k = 0; k < K; k++) {
    for (int t = 0; t < T + B_r; t++) { 
      lambda_bar_k(t,k) = 1.0;
      if (t < T) {
        u_bar_k(t,k) = u_k[k] + (T + B_r - t) / 2.0;
        v_bar_k(t,k) = u_bar_k(t,k);
        pi_bar_k(t,k) =  (double) 1 / T;
      }
    }
  }
  
  // initialize mu_0
  if (fit_intercept) {
    if (B_l > 0) mu_0 = Rcpp::mean(y_0);
    else mu_0 = Rcpp::mean(y[Rcpp::Range(0,4)]);
  } 
  
  // initialize lambda_0
  if (fit_scale) {
    if (B_l > 1) lambda_0 = 1 / Rcpp::var(y_0);
    else lambda_0 = 1 / Rcpp::var(y[Rcpp::Range(0,4)]);
  } 
  
  // initialize residual and scale vector
  NumericVector r_tilde = y - mu_0;
  NumericVector lambda_bar (T + B_r, lambda_0);
  
  // initialize correction term
  NumericVector delta (T + B_r, 0.0);
      
  // initialize ELBO
  NumericVector elbo (max_iter, R_NegInf);
  
  // tracking terms 
  double signal_rev_sum = 0.0, signal_fwd_sum = 0.0;
  double scale_rev_sum = 0.0;
  double tot_prob = 0.0, log_pi_max = R_NegInf;
  
  // vb algorithm
  int iter = 0;
  while (iter < max_iter) {
    
    iter++;
    if (verbose & (iter % 1000 == 0)) Rcout << "Iteration: " << iter << "\n"; 
      
    // updating q(b_j, s_j, gamma_j)
    // if (J > 0) {
    //   for (int j = 0; j < J; j++) {
    //     // deleting j^{th} component from residual terms
    //     for (int t = 0; t < T + B_r; t++) {
    //       r_tilde[t] += mu_lambda_j(t,j) / lambda_bar_j(t,j);
    //       lambda_bar[t] /= lambda_bar_j(t,j);
    //       delta[t] += (mu_lambda_j(t,j) / lambda_bar_j(t,j) - mu2_lambda_j(t,j) / lambda_bar_j(t,j)) / lambda_bar_j(t,j);
    //     }
    // 
    //   }
    // }
    
    for (int l = 0; l < L; l++) {
      for (int t = T + B_r - 1; t >= 0; t--) {
        // adding back l^{th} component from residual terms
        r_tilde[t] += mu_bar_l(t,l); 
        delta[t] += std::pow(mu_bar_l(t,l), 2) - mu2_bar_l(t,l);
        
        signal_rev_sum += lambda_bar[t] * r_tilde[t];
        scale_rev_sum += lambda_bar[t];
        
        if (t < T) {
          tau_bar_l(t,l) = tau_l[l] + scale_rev_sum;
          b_bar_l(t,l) = signal_rev_sum / tau_bar_l(t,l);
          pi_bar_l(t,l) = std::log(pi_l(t,l)) - 0.5 * (std::log(tau_bar_l(t,l)) - tau_bar_l(t,l) * std::pow(b_bar_l(t,l), 2));
          if (pi_bar_l(t,l) > log_pi_max) log_pi_max = pi_bar_l(t,l);
        }
      }
      
      for(int t = 0; t < T; t++) {
        pi_bar_l(t,l) = std::exp(pi_bar_l(t,l) - log_pi_max);
        tot_prob += pi_bar_l(t,l);
      }
      
      pi_bar_l(0,l) /= tot_prob;
      mu_bar_l(0,l) = b_bar_l(0,l) * pi_bar_l(0,l);
      mu2_bar_l(0,l) = (std::pow(b_bar_l(0,l), 2) + 1 / tau_bar_l(0,l)) * pi_bar_l(0,l);
      
      for (int t = 1; t < T + B_r; t++) {
        if (t < T) {
          pi_bar_l(t,l) /= tot_prob;
          mu_bar_l(t,l) = b_bar_l(t,l) * pi_bar_l(t,l) + mu_bar_l(t-1,l);
          mu2_bar_l(t,l) = (std::pow(b_bar_l(t,l), 2) + 1 / tau_bar_l(t,l)) * pi_bar_l(t,l) + mu2_bar_l(t-1,l);
        } else {
          mu_bar_l(t,l) = mu_bar_l(T-1,l);
          mu2_bar_l(t,l) = mu2_bar_l(T-1,l);
        }
        // subtracting l^{th} component from residual terms
        r_tilde[t] -= mu_bar_l(t,l); 
        delta[t] -= std::pow(mu_bar_l(t,l), 2) - mu2_bar_l(t,l);
      }
      
      // reset tracking variables
      signal_rev_sum = 0.0, tot_prob = 0.0, scale_rev_sum = 0.0, log_pi_max = R_NegInf;
    }
    
    for (int k = 0; k < K; k++) {
      
      pi_bar_k(_,k) = Rcpp::rep(T, 0.0);
      
      for (int t = T + B_r - 1; t >= 0; t--) {
        
        // dividing out k^{th} scale component
        lambda_bar[t] /= lambda_bar_k(t,k);
        
        signal_rev_sum += lambda_bar[t] * (std::pow(r_tilde[t], 2) + delta[t]);
        
        if (t < T) {
          v_bar_k(t,k) = v_k[k] + 0.5 * signal_rev_sum;
          pi_bar_k(t,k) += std::log(pi_k(t,k)) + std::lgamma(u_bar_k(t,k)) - u_bar_k(t,k) * std::log(v_bar_k(t,k));
          pi_bar_k(T-t-1,k) -= 0.5 * signal_fwd_sum;
          if (t < half) {
            if (pi_bar_k(T-t-1,k) > log_pi_max) log_pi_max = pi_bar_k(T-t-1,k);
            if (pi_bar_k(t,k) > log_pi_max) log_pi_max = pi_bar_k(t,k);
          }
          signal_fwd_sum += lambda_bar[T-t-1] * (std::pow(r_tilde[T-t-1], 2) + delta[T-t-1]);
        }
      }
      
      for(int t = 0; t < T; t++) {
        pi_bar_k(t,k) = std::exp(pi_bar_k(t,k) - log_pi_max);
        tot_prob += pi_bar_k(t,k);
      }
      
      // reusing trackers
      signal_fwd_sum = 0;
      signal_rev_sum = 1.0;
      
      for (int t = 0; t < T + B_r; t++) {
        if (t < T) {
          pi_bar_k(t,k) /= tot_prob;
          signal_rev_sum -= pi_bar_k(t,k);
          signal_fwd_sum += (u_bar_k(t,k) / v_bar_k(t,k)) * pi_bar_k(t,k);
          lambda_bar_k(t,k) = signal_fwd_sum + signal_rev_sum;
        } else {
          lambda_bar_k(t,k) = lambda_bar_k(T-1,k);
        }
        // multiplying back k^{th} scale component 
        lambda_bar[t] *= lambda_bar_k(t,k);
      }
      
      // reset tracking variables
      signal_rev_sum = 0.0, signal_fwd_sum = 0.0, tot_prob = 0.0, log_pi_max = R_NegInf;
    }
  }
    
  
  List J_model;
  List L_model;
  List K_model;
  
  if (J > 0) {
    J_model =  List::create(_["pi"] = pi_bar_j, _["b"] = b_bar_j, _["tau"] = tau_bar_j, _["v"] = v_bar_j, _["u"] = u_bar_j);
  }
  if (L > 0) {
    L_model =  List::create(_["pi"] = pi_bar_l, _["b"] = b_bar_l, _["tau"] = tau_bar_l);
  }
  if (K > 0) {
    K_model =  List::create(_["pi"] = pi_bar_k, _["v"] = v_bar_k, _["u"] = u_bar_k);
  }
  
  List result = List::create(_["elbo"] = elbo, _["lambda"] = lambda_bar, _["J_model"] = J_model, _["L_model"] = L_model, _["K_model"] = K_model);
  return result;
}


/*** R
T = 100; B_l = 1; B_r = 1;
lambda = 1 / sqrt(c(rep(1, 50), rep(4, 33), rep(0.25, 36)))
mu = c(rep(0, 33), rep(4, 33), rep(1, 36))
y = rnorm(T+B_r+B_l, mean = mu, sd = lambda); J = 0; L = 2; K = 2; 
fit_intercept = FALSE; fit_scale = FALSE; max_iter = 10000;
tau_j = rep(0.1, J); u_j = tau_j; v_j = u_j; pi_j = matrix(1/T, T, J);
tau_l = rep(0.1, L); pi_l = matrix(1/T, T, L);
u_k = rep(0.1, K); v_k = u_k; pi_k = matrix(1/T, T, K);

fit = mich(y, J, L, K, 
           B_l, B_r, fit_intercept, fit_scale, 
           max_iter, verbose = TRUE,
           tau_j, u_j, v_j, pi_j,
           tau_l, pi_l,
           u_k, v_k, pi_k)
tmp=michR(y, J, L, K, fit.intercept = FALSE, fit.scale = FALSE, max_iter = max_iter, tol = -Inf, B_l = B_l, B_r=B_r)
matplot(fit$K_model$pi)
points(tmp$scale.model$probs[2:101,1],col="red")
points(tmp$scale.model$probs[2:101,2])
# microbenchmark(
# mich(y, J, L, K,
#      B_l, B_r, fit_intercept, fit_scale, max_iter, verbose = FALSE,
#      tau_j, u_j, v_j, pi_j,
#      tau_l, pi_l,
#      u_k, v_k, pi_k),
# michR(y, J, L, K, max_iter = max_iter, tol = -Inf))
*/
