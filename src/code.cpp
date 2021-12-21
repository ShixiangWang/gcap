#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calc_dist(NumericMatrix x) {
  int n = x.nrow();
  NumericMatrix out(n, n); // output a cosine value matrix

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (x(i, 0) == x(j, 0)) {
        out(i, j) = abs(x(i, 1) - x(j, 1));
      } else {
        // Set a hard distance for genes not in same chromosome
        out(i, j) = 1e8;
      }
    }
  }

  return out;
}
