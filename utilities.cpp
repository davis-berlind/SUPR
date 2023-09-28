#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector revcumsum(NumericVector x) {
  for(int i = x.length()-1; i > 0; i--) {
    x[i-1] += x[i];
  }
  return x;
}

// [[Rcpp::export]]
NumericVector test(NumericVector y, int B_l) {
  NumericVector y_0(B_l);
  if (B_l > 0) {
    y_0 = y[Rcpp::Range(0, B_l - 1)];
    y.erase(0,B_l);
  }
  return y;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test(1:10, 0)
*/
