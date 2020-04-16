#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Auxiliar function to calculate the Covariogram from Lavancier et al. 2019
//'
//' @param h \code{numeric scalar}
//' @param dists \code{numeric vector}
//' @param resp \code{numeric vector}
//' @param tol \code{numeric scalar}
//' @param n_xy \code{numeric scalar}
//'
//' @description Auxiliar function to calculate \eqn{\hat{H_{12}}_i(d)}.
//'
//' @return \code{numeric vector}
//'
// [[Rcpp::export]]
double covariog_aux(const double h, const arma::vec dists,
                    const arma::vec resp, const double tol, const int n_xy) {
  double out = 0;
  double p_hat = arma::mean(resp);
  int id_dist = 0;
  int n_h = 0;

  if(h == 0) {
    n_h = resp.n_elem;
    for(int j = 0; j < n_xy; j++) {
      out += (resp.at(j) - p_hat)*(resp.at(j) - p_hat);
    }
  } else {
    for(int j = 1; j < n_xy; j++) {
      for(int i = 0; i < j; i++) {
        id_dist = (n_xy - 1)*(i) - (((i + 1)*i)/2) + j - 1;
        if(dists(id_dist) >= h - tol && dists(id_dist) <= h + tol)
          // n_h = i + 1;
          n_h += 1;
        out += (resp.at(j) - p_hat)*(resp.at(i) - p_hat);
      }
    }
  }

  if(n_h > 0)
    out = out/n_h;
  else
    out = 0;

  return out;
}
