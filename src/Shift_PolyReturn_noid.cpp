#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Create copies of a set of polygons
//'
//' @param obj_sp object from class \code{SpatialPolygons}
//' @param bbox_tot Boundary box from class \code{matrix}
//'
//' @description Auxiliar function.
//'
//' @importFrom stats runif
//' @return an object from class \code{SpatialPolygons}
//'
// [[Rcpp::export]]
S4 poly_shift_noid(S4& obj_sp, NumericMatrix& bbox_tot) {

  // Verifying if there are a bbox


  obj_sp.slot("bbox") = bbox_tot;


  // Setting object bbox

  NumericMatrix bbox_obj = obj_sp.slot("bbox");

  // Calculating the ranges of the bbox

  double range_x = bbox_obj(0,1);
  double range_y = bbox_obj(1,1);

  range_x = (range_x - bbox_obj(0,0));
  range_y = (range_y - bbox_obj(1,0));

  // object to be changed and features

  List Aux = obj_sp.slot("polygons");

  int n_poly = Aux.size();

  IntegerVector pOrder = seq_len(n_poly*4);

  S4 out = Aux[0];

  String ID = (1);

  // ID.push_front("a ");

  out.slot("ID") = ID;

  List Aux_Final = List::create(out);

  // Setting the polygons IDs

  for(int i = 1; i < n_poly; i++)
  {
    S4 out1 = Aux[i];

    String ID = (i + 1);

    // ID.push_front("a ");

    out1.slot("ID") = ID;

    Aux_Final.push_back(out1);
  }

  // obj_sp2 = copy of obj_sp to be X translated

  S4 obj_sp2 = clone(obj_sp);

  List Aux2 = obj_sp2.slot("polygons");

  // Changing obj_sp2

  for(int i = 0; i < n_poly; i++)
  {
    S4 out2 = Aux2[i];
    List Aux2_2 = out2.slot("Polygons");
    for(int j = 0; j < Aux2_2.size(); j++) {
      S4 out2_2 = Aux2_2[j];

      NumericMatrix mat_aux = out2_2.slot("coords");

      mat_aux(_,0) = mat_aux(_,0) + range_x;

      out2_2.slot("coords") = mat_aux;

      NumericVector vec_aux = out2_2.slot("labpt");

      vec_aux(0) = vec_aux(0) + range_x;

      out2_2.slot("labpt") = vec_aux;

      Aux2_2[j] = out2_2;
    }

    out2.slot("Polygons") = Aux2_2;

    NumericVector vec_aux = out2.slot("labpt");

    vec_aux(0) = vec_aux(0) + range_x;

    out2.slot("labpt") = vec_aux;

    String ID = 2*(i + 1);

    // ID.push_front("a ");

    out2.slot("ID") = ID;

    Aux2[i] = out2;

    Aux_Final.push_back(out2);
  }

  obj_sp2.slot("polygons") = Aux2;

  // obj_sp3 = copy of obj_sp - Y translated

  S4 obj_sp3 = clone(obj_sp);

  List Aux3 = obj_sp3.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out3 = Aux3[i];
    List Aux2_3 = out3.slot("Polygons");
    for(int j = 0; j < Aux2_3.size(); j++) {
      S4 out2_3 = Aux2_3[j];

      NumericMatrix mat_aux = out2_3.slot("coords");

      mat_aux(_,1) = mat_aux(_,1) + range_y;

      out2_3.slot("coords") = mat_aux;

      NumericVector vec_aux = out2_3.slot("labpt");

      out2_3.slot("labpt") = vec_aux;

      vec_aux(1) = vec_aux(1) + range_y;

      Aux2_3[j] = out2_3;
    }

    NumericVector vec_aux = out3.slot("labpt");

    vec_aux(1) = vec_aux(1) + range_y;

    out3.slot("Polygons") = Aux2_3;

    out3.slot("labpt") = vec_aux;

    String ID = 3*(i + 1);

    // ID.push_front("a ");

    out3.slot("ID") = ID;

    Aux_Final.push_back(out3);

    Aux3[i] = out3;
  }

  obj_sp3.slot("polygons") = Aux3;

  // obj_sp4 = copy of obj_sp - X and Y translated

  S4 obj_sp4 = clone(obj_sp);

  List Aux4 = obj_sp4.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out4 = Aux4[i];
    List Aux2_4 = out4.slot("Polygons");
    for(int j = 0; j < Aux2_4.size(); j++) {
      S4 out2_4 = Aux2_4[j];

      NumericMatrix mat_aux = out2_4.slot("coords");


      mat_aux(_,0) = mat_aux(_,0) + range_x;
      mat_aux(_,1) = mat_aux(_,1) + range_y;

      out2_4.slot("coords") = mat_aux;

      NumericVector vec_aux = out2_4.slot("labpt");

      vec_aux(0) = vec_aux(0) + range_x;
      vec_aux(1) = vec_aux(1) + range_y;


      out2_4.slot("labpt") = vec_aux;

      Aux2_4[j] = out2_4;
    }

    NumericVector vec_aux = out4.slot("labpt");

    vec_aux(0) = vec_aux(0) + range_x;
    vec_aux(1) = vec_aux(1) + range_y;

    out4.slot("Polygons") = Aux2_4;

    out4.slot("labpt") = vec_aux;

    String ID = 4*(i + 1);

    // ID.push_front("a ");

    out4.slot("ID") = ID;

    Aux_Final.push_back(out4);

    Aux4[i] = out4;
  }

  obj_sp4.slot("polygons") = Aux4;

  // Defining the Final Output

  NumericMatrix mat_aux = obj_sp.slot("bbox");
  NumericMatrix box_aux(2,2);

  box_aux = clone(mat_aux);

  box_aux(0,1) = box_aux(0,1) + range_x;
  box_aux(1,1) = box_aux(1,1) + range_y;

  CharacterVector rNames = CharacterVector::create("x", "y");
  rownames(box_aux) = rNames;

  CharacterVector cNames = CharacterVector::create("min", "max");
  colnames(box_aux) = cNames;

  S4 output = clone(obj_sp);

  // S4 projString("CRS");

  // String charac_aux = NA_STRING;

  // projString.slot("projargs") = charac_aux;

  // output.slot("proj4string") = projString;

  output.slot("plotOrder") =  pOrder;

  output.slot("bbox") = box_aux;

  output.slot("polygons") = Aux_Final;

  return output;
}

