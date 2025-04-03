#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalMatrix force_adj(LogicalMatrix M) {
  int N = M.ncol();
  if (N <= 2) return M;
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        if (M(i, k) & M(k, j)) {
          M(i,j) = true;
        }
      }
    }
  }
  
  return M;
}