#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Auxiliar function
//'
//' @param obj_sp a spatial object
//' @param obj_sp2 a spatial object
//' @param obj_sp3 a spatial object
//' @param obj_sp4 a spatial object
//' @param n_poly an integer
//' @param range_x a numeric
//' @param range_y a numeric
//'
// [[Rcpp::export]]
S4 shift_aux(S4 obj_sp, S4 obj_sp2, S4 obj_sp3, S4 obj_sp4, int n_poly,
             double range_x, double range_y) {

  // object to be changed and features

  List Aux = obj_sp.slot("polygons");

  IntegerVector pOrder = seq_len(n_poly*4);

  S4 out = Aux[0];

  String ID = 1;
  out.slot("ID") = ID;

  List Aux_Final = List::create(out);

  // Setting the polygons IDs

  for(int i = 1; i < n_poly; i++)
  {
    S4 out1 = Aux[i];
    String ID = (i + 1);
    out1.slot("ID") = ID;
    Aux_Final.push_back(out1);
  }

  Aux = obj_sp2.slot("polygons");

  // Changing obj_sp2

  for(int i = 0; i < n_poly; i++)
  {
    S4 out2 = Aux[i];

    String ID = (i + 1);
    out2.slot("ID") = ID;

    Aux_Final.push_back(out2);
  }

  Aux = obj_sp3.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out3 = Aux[i];
    String ID = (i + 1);

    out3.slot("ID") = ID;

    Aux_Final.push_back(out3);
  }

  Aux = obj_sp4.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out4 = Aux[i];
    String ID = (i + 1);

    out4.slot("ID") = ID;

    Aux_Final.push_back(out4);
  }

  // Defining the Final Output

  NumericMatrix mat_aux = obj_sp.slot("bbox");
  NumericMatrix box_aux(2,2);

  box_aux = clone(mat_aux);

  box_aux(0,1) = box_aux(0,1) + range_x;
  box_aux(1,1) = box_aux(1,1) + range_y;

  // CharacterVector rNames("x", "y");
  // StringVector rNames(2);
  // rNames[0] = "x";
  // rNames[1] = "y";
  // rownames(box_aux) = CharacterVector::create("x", "y");

  CharacterVector rNames = CharacterVector::create("x", "y");
  rownames(box_aux) = rNames;

  // CharacterVector cNames("min", "max");
  // StringVector cNames(2);
  // cNames[0] = "min";
  // cNames[1] = "max";
  // colnames(box_aux) = CharacterVector::create("min", "max");

  CharacterVector cNames = CharacterVector::create("min", "max");
  colnames(box_aux) = cNames;

  S4 output = clone(obj_sp);

  S4 projString("CRS");

  String charac_aux = NA_STRING;

  projString.slot("projargs") = charac_aux;

  output.slot("proj4string") = projString;

  output.slot("plotOrder") =  pOrder;

  output.slot("bbox") = box_aux;

  output.slot("polygons") = Aux_Final;

  return output;
}
