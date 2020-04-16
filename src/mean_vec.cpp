#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Vector of means
//'
//' @param x \code{numeric vector}
//'
//' @description Auxiliar function to calculate \eqn{\hat{H_{12}}_i(d)}.
//'
//' @return \code{numeric vector}
//'
// [[Rcpp::export]]
arma::vec mean_vec(arma::vec& x) {
  arma::uvec idx(1);
  arma::uword n = x.size();
  arma::vec y(n);
  arma::vec idx_aux = arma::linspace(0, (n - 1), n);

  for(arma::uword i = 0; i < n; i++) {
    idx = arma::find(idx_aux != idx_aux.at(i));
    y.at(i) = arma::mean(x.elem(idx));
  }
  return y;
}
