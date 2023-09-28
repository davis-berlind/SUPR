#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List mich(NumericVector y, int J, int L, int K, 
          int B_l, int B_r, bool fit_intercept, bool fit_scale, int max_iter, 
          NumericVector tau_j, NumericVector u_j, NumericVector v_j, NumericMatrix pi_j,
          NumericVector tau_l, NumericMatrix pi_l,
          NumericVector u_k, NumericVector v_k, NumericMatrix pi_k) {
  
  // initialize intercept terms
  double mu_0 = 0, lambda_0 = 1;
  
  // extract observations in start buffer
  NumericVector y_0 (B_l);
  if (B_l > 0) {
    y_0 = y[Rcpp::Range(0, B_l - 1)];
    y.erase(0,B_l);
  }
  int T = y.length() - B_r;
  
  // initialize J components
  List smscp_fit;
  NumericMatrix b_bar_j (T, J);
  NumericMatrix tau_bar_j (T, J);
  NumericMatrix u_bar_j (T, J);
  NumericMatrix v_bar_j (T, J);
  NumericMatrix pi_bar_j (T, J);
  NumericMatrix mu_lambda_j (T + B_r, J);
  NumericMatrix mu2_lambda_j (T + B_r, J);
  NumericMatrix lambda_bar_j (T + B_r, J);
  
  if (J > 0) {
    for (int t = 0; t < T + B_r; t++) {
      for (int j = 0; j < J; j++) {
        tau_bar_j(t,j) = 1;
        lambda_bar_j(t,j) = 1;
        if (t < T) {
          u_bar_j(t,j) = u_j[j] + (T + B_r - t) / 2;
          v_bar_j(t,j) = u_bar_j(t,j); 
          pi_bar_j(t,j) = 1 / T;
        }
      }
    }
  }
  
  // initialize L components
  NumericMatrix b_bar_l (T, L);
  NumericMatrix tau_bar_l (T, L);
  NumericMatrix pi_bar_l (T, L);
  NumericMatrix mu_bar_l (T + B_r, L);
  NumericMatrix mu2_bar_l (T + B_r, L);

  if (L > 0) {
    for (int t = 0; t < T + B_r; t++) {
      for (int l = 0; l < L; l++) {
        tau_bar_l(t,l) = 1;
        if (t < T) {
          pi_bar_l(t,l) = 1 / T;
        }
      }
    }
  }
  
  // initialize K components
  NumericMatrix u_bar_k (T, K);
  NumericMatrix v_bar_k (T, K);
  NumericMatrix pi_bar_k (T, K);
  NumericMatrix lambda_bar_k (T + B_r, K);
  
  if (K > 0) {
    for (int t = 0; t < T + B_r; t++) {
      for (int k = 0; k < K; k++) {
        lambda_bar_k(t,k) = 1;
        if (t < T) {
          u_bar_k(t,k) = u_k[k] + (T + B_r - t) / 2;
          v_bar_k(t,k) = u_bar_k(t,k); 
          pi_bar_k(t,k) = 1 / T;
        }
      }
    }
  }
  
  // initialize mu_0
  if (fit_intercept) {
    if (B_l > 0) mu_0 = Rcpp::mean(y_0);
    else mu_0 = y[0];
  } 
  
  // initialize lambda_0
  if (fit_scale) {
    if (B_l == 1) lambda_0 = 1;
    else if (B_l > 1) lambda_0 = 1 / Rcpp::var(y_0);
    else lambda_0 = 1 / Rcpp::var(y[Rcpp::Range(0,1)]);
  } 
  
  // initialize residual and scale vector
  NumericVector r_tilde = y - mu_0;
  NumericVector lambda_bar (T + B_r, lambda_0);
  
  // initialize correction term
  NumericVector delta (T + B_r, 0.0);
      
  // initialize ELBO
  NumericVector elbo (max_iter, R_NegInf);
  
  // vb algorithm
  int iter = 0;
  while (iter < max_iter) {
    // updating q(b_j, s_j, gamma_j)
    if (J > 0) {
      for (int j = 0; j < J; j++) {
        // deleting j^{th} component from residual terms
        for (int t = 0; t < T + B_r; t++) {
          r_tilde[t] += mu_lambda_j(t,j) / lambda_bar_j(t,j); 
          lambda_bar[t] /= lambda_bar_j(t,j);
          delta[t] += (mu_lambda_j(t,j) / lambda_bar_j(t,j) - mu2_lambda_j(t,j) / lambda_bar_j(t,j)) / lambda_bar_j(t,j);
        }

      }
    }
  }
  
  List result = List::create(_["elbo"] = elbo);
  return result;
}


/*** R
T = 100; B_l = 1; B_r = 1;
y = rnorm(T + B_l + B_r); J = 1; L = 1; K = 1; 
fit_intercept = TRUE; fit_scale = TRUE; max_iter = 100;
tau_j = rep(0.1, J); u_j = tau_j; v_j = u_j; pi_j = matrix(1/T, T, J);
tau_l = rep(0.1, L); pi_l = matrix(1/T, T, L);
u_k = rep(0.1, K); v_k = u_k; pi_k = matrix(1/T, T, K);

mich(y, J, L, K, 
     B_l, B_r, fit_intercept, fit_scale, max_iter,
     tau_j, u_j, v_j, pi_j,
     tau_l, pi_l,
     u_k, v_k, pi_k)
*/
