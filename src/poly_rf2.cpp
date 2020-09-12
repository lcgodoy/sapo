#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Polygons' Random Shift - 2
//'
//' @param objsp object from class \code{SpatialPolygons}
//'
//' @return an object from class \code{SpatialPolygons} randomly translated
//'
// [[Rcpp::export]]
S4 poly_rf2(const S4& objsp) {

  NumericMatrix bbox_obj = objsp.slot("bbox");

  S4 output = clone(objsp);

  double range_x = bbox_obj(0,1);
  double range_y = bbox_obj(1,1);

  range_x = (range_x - bbox_obj(0, 0));
  range_y = (range_y - bbox_obj(1, 0));

  double n_x = Rcpp::runif(1)[0];
  double n_y = Rcpp::runif(1)[0];

  n_x = bbox_obj(0, 0) + (bbox_obj(0, 1) - bbox_obj(0, 0))*n_x;
  n_y = bbox_obj(1, 0) + (bbox_obj(1, 1) - bbox_obj(1, 0))*n_y;

  NumericMatrix bbox_new(2,2);

  IntegerVector side_aux = seq(1, 4);
  int side = sample(side_aux, 1, false)[0];

  double jump_x, jump_y;
  
  if(side == 1) {
    bbox_new(0, 0) = n_x;
    bbox_new(0, 1) = n_x + range_x;
    bbox_new(1, 0) = n_y;
    bbox_new(1, 1) = n_y + range_y;
  } else if(side == 2) {
     bbox_new(0, 0) = n_x - range_x;
     bbox_new(0, 1) = n_x;
     bbox_new(1, 0) = n_y;
     bbox_new(1, 1) = n_y + range_y;
  } else if(side == 3) {
     bbox_new(0, 0) = n_x - range_x;
     bbox_new(0, 1) = n_x;
     bbox_new(1, 0) = n_y - range_y;
     bbox_new(1, 1) = n_y;
  } else {
     bbox_new(0, 0) = n_x;
     bbox_new(0, 1) = n_x + range_x;
     bbox_new(1, 0) = n_y - range_y;
     bbox_new(1, 1) = n_y;
  }

  jump_x = bbox_new(0, 0) - bbox_obj(0, 0);
  jump_y = bbox_new(1, 0) - bbox_obj(1, 0);
  
  CharacterVector rNames = CharacterVector::create("x", "y");
  rownames(bbox_new) = rNames;

  CharacterVector cNames = CharacterVector::create("min", "max");
  colnames(bbox_new) = cNames;

  output.slot("bbox") = bbox_new;

  List Aux = output.slot("polygons");

  int n_poly = Aux.size();

  for(int i = 0; i < n_poly; i++) {
    S4 out = Aux[i];
    List Aux2 = out.slot("Polygons");

    for(int j = 0; j < Aux2.size(); j++) {
      S4 out2 = Aux2[j];

      NumericVector vec_aux = out2.slot("labpt");

      vec_aux(0) = vec_aux(0) + jump_x;
      vec_aux(1) = vec_aux(1) + jump_y;

      NumericMatrix mat_aux = out2.slot("coords");

      mat_aux(_,0) = mat_aux(_,0) + jump_x;
      mat_aux(_,1) = mat_aux(_,1) + jump_y;

      out2.slot("coords") = mat_aux;
      out2.slot("labpt") = vec_aux;

      Aux2[j] = out2;
    }

    NumericVector vec_aux = out.slot("labpt");

    vec_aux(0) = vec_aux(0) + jump_x;
    vec_aux(1) = vec_aux(1) + jump_y;

    out.slot("Polygons") = Aux2;

    out.slot("labpt") = vec_aux;

    Aux[i] = out;
  }

  // objsp.slot("polygons") = Aux;

  output.slot("polygons") = Aux;

  return(output);
}
