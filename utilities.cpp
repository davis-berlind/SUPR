#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector revcumsum(NumericVector x) {
  for(int i = x.length()-1; i > 0; i--) {
    x[i-1] += x[i];
  }
  return x;
}

// [[Rcpp::export]]
double test(NumericVector y_0) {
  double out = Rcpp::sum(y_0);
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test(numeric(0)) 
*/
