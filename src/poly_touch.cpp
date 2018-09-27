// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

//' Polygons that touch a bbox
//'
//' @param x a \code{SpatialPolygon}
//' @param bbox a \code{numeric matrix}
//'
//' @export
//'
// [[Rcpp::export]]
S4 poly_touch(S4& x, arma::mat& bbox) {
  S4 output = clone(x);
  List aux = output.slot("polygons");
  List out_poly;
  int n = aux.length();

  for(int i = 0; i < n; i++) {
    S4 poly_aux = aux[i];
    arma::vec cent_aux = poly_aux.slot("labpt");
    bool logic1 = false;
    if(cent_aux.at(0) >= bbox(0, 0) && cent_aux.at(0) <= bbox(0, 1) &&
       cent_aux.at(1) >= bbox(1, 0) && cent_aux.at(1) <= bbox(1, 1)) {
      logic1 = true;
      S4 aux2 = aux[i];
      out_poly.push_front(aux2);
    }
    if(logic1 == false) {
      List aux_list = poly_aux.slot("Polygons");
      S4 aux_polygons = aux_list[0];
      arma::mat coord_aux = aux_polygons.slot("coords");
      arma::vec x_aux = coord_aux.col(0);
      arma::vec y_aux = coord_aux.col(1);
      bool logic2 = (arma::any(x_aux >= bbox(0, 0)) &&
                     arma::any(x_aux <= bbox(0, 1)) &&
                     arma::any(y_aux >= bbox(1, 0)) &&
                     arma::any(y_aux <= bbox(1, 1)));
      if(logic2 == true) {
        S4 aux2 = aux[i];
        out_poly.push_front(aux2);
      }
    }
  }

  output.slot("polygons") = out_poly;

  return output;
}

