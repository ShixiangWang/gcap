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

// [[Rcpp::export]]
NumericVector vec_shift(NumericVector arr) {
  int SIZE = arr.size();

  double last = arr[SIZE - 1];
  for (int i = SIZE - 1; i > 0; i--)
    arr[i] = arr[i - 1];

  arr[0] = last;

  return arr;
}


// [[Rcpp::export]]
NumericVector conv_binomial(NumericVector p) {
  int n_p = p.size();
  NumericVector z(n_p + 1, 0.0);
  NumericVector z2(n_p + 1, 0.0);
  z[0] = 1.0;
  double pi = 0.0;

  //Rcout << "The value of z : " << z << "\n";

  IntegerVector rg = Range(0, n_p - 1);
  for (int i=0; i<n_p; i++) {
    pi = p[i];
    //z2 = clone(z);
    //z2 = vec_shift(z2);
    z2[0] = 0.0;
    z2[rg+1] = z[rg];
    //Rcout << "The value of z : " << z << "\n";
    //Rcout << "The value of z2 : " << z2 << "\n";
    z = (1 - pi) * z + pi * z2;
  }
  return z;
}

