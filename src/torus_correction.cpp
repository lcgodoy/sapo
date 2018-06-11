#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Toroidal edge correction
//'
//' @param objsp a \code{SpatialPolygon}
//' @param bbox_tot a \code{numeric matrix}
//'
//' @export
//'
// [[Rcpp::export]]
S4 torus_corr(S4& objsp, Rcpp::Nullable<Rcpp::NumericMatrix> bbox_tot = R_NilValue) {

  // Verifying if there are a bbox

  if (bbox_tot.isNotNull())
  {
    objsp.slot("bbox") = bbox_tot;
  }

  // Setting object bbox

  NumericMatrix bbox_obj = objsp.slot("bbox");

  // Calculating the ranges of the bbox

  double range_x = bbox_obj(0,1);
  double range_y = bbox_obj(1,1);

  range_x = (range_x - bbox_obj(0,0));
  range_y = (range_y - bbox_obj(1,0));

  // object to be changed and features

  List Aux = objsp.slot("polygons");

  int n_poly = Aux.size();

  IntegerVector pOrder = seq_len(n_poly*9);

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

  // Objsp2 = copy of objsp to be X translated

  S4 objsp2 = clone(objsp);

  List Aux2 = objsp2.slot("polygons");

  // Changing objsp2

  for(int i = 0; i < n_poly; i++)
  {
    S4 out2 = Aux2[i];
    List Aux2_2 = out2.slot("Polygons");
    S4 out2_2 = Aux2_2[0];

    NumericMatrix mat_aux = out2_2.slot("coords");

    mat_aux(_,0) = mat_aux(_,0) + range_x;

    out2_2.slot("coords") = mat_aux;

    NumericVector vec_aux = out2.slot("labpt");

    vec_aux(0) = vec_aux(0) + range_x;

    Aux2_2[0] = out2_2;

    out2.slot("Polygons") = Aux2_2;

    out2.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out2.slot("ID") = ID;

    Aux2[i] = out2;

    Aux_Final.push_back(out2);
  }

  objsp2.slot("polygons") = Aux2;

  // Objsp3 = copy of objsp - Y translated

  S4 objsp3 = clone(objsp);

  List Aux3 = objsp3.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out3 = Aux3[i];
    List Aux2_3 = out3.slot("Polygons");
    S4 out2_3 = Aux2_3[0];

    NumericMatrix mat_aux = out2_3.slot("coords");

    mat_aux(_,1) = mat_aux(_,1) + range_y;

    out2_3.slot("coords") = mat_aux;

    NumericVector vec_aux = out3.slot("labpt");

    vec_aux(1) = vec_aux(1) + range_y;

    Aux2_3[0] = out2_3;

    out3.slot("Polygons") = Aux2_3;

    out3.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out3.slot("ID") = ID;

    Aux_Final.push_back(out3);

    Aux3[i] = out3;
  }

  objsp3.slot("polygons") = Aux3;


  // Objsp4 = copy of objsp - X and Y translated

  S4 objsp4 = clone(objsp);

  List Aux4 = objsp4.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out4 = Aux4[i];
    List Aux2_4 = out4.slot("Polygons");
    S4 out2_4 = Aux2_4[0];

    NumericMatrix mat_aux = out2_4.slot("coords");


    mat_aux(_,0) = mat_aux(_,0) + range_x;
    mat_aux(_,1) = mat_aux(_,1) + range_y;

    out2_4.slot("coords") = mat_aux;

    NumericVector vec_aux = out4.slot("labpt");

    vec_aux(0) = vec_aux(0) + range_x;
    vec_aux(1) = vec_aux(1) + range_y;

    Aux2_4[0] = out2_4;

    out4.slot("Polygons") = Aux2_4;

    out4.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out4.slot("ID") = ID;

    Aux_Final.push_back(out4);

    Aux4[i] = out4;
  }

  objsp4.slot("polygons") = Aux4;

  S4 objsp5 = clone(objsp);

  List Aux5 = objsp5.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out5 = Aux5[i];
    List Aux2_5 = out5.slot("Polygons");
    S4 out2_5 = Aux2_5[0];

    NumericMatrix mat_aux = out2_5.slot("coords");


    mat_aux(_,0) = mat_aux(_,0) - range_x;
    mat_aux(_,1) = mat_aux(_,1) - range_y;

    out2_5.slot("coords") = mat_aux;

    NumericVector vec_aux = out5.slot("labpt");

    vec_aux(0) = vec_aux(0) - range_x;
    vec_aux(1) = vec_aux(1) - range_y;

    Aux2_5[0] = out2_5;

    out5.slot("Polygons") = Aux2_5;

    out5.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out5.slot("ID") = ID;

    Aux_Final.push_back(out5);

    Aux5[i] = out5;
  }

  objsp5.slot("polygons") = Aux5;

  S4 objsp6 = clone(objsp);

  List Aux6 = objsp6.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out6 = Aux6[i];
    List Aux2_6 = out6.slot("Polygons");
    S4 out2_6 = Aux2_6[0];

    NumericMatrix mat_aux = out2_6.slot("coords");

    mat_aux(_,1) = mat_aux(_,1) - range_y;

    out2_6.slot("coords") = mat_aux;

    NumericVector vec_aux = out6.slot("labpt");

    vec_aux(1) = vec_aux(1) - range_y;

    Aux2_6[0] = out2_6;

    out6.slot("Polygons") = Aux2_6;

    out6.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out6.slot("ID") = ID;

    Aux_Final.push_back(out6);

    Aux6[i] = out6;
  }

  objsp6.slot("polygons") = Aux6;

  S4 objsp7 = clone(objsp);

  List Aux7 = objsp7.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out7 = Aux7[i];
    List Aux2_7 = out7.slot("Polygons");
    S4 out2_7 = Aux2_7[0];

    NumericMatrix mat_aux = out2_7.slot("coords");

    mat_aux(_,0) = mat_aux(_,0) - range_x;

    out2_7.slot("coords") = mat_aux;

    NumericVector vec_aux = out7.slot("labpt");

    vec_aux(0) = vec_aux(0) - range_x;

    Aux2_7[0] = out2_7;

    out7.slot("Polygons") = Aux2_7;

    out7.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out7.slot("ID") = ID;

    Aux_Final.push_back(out7);

    Aux7[i] = out7;
  }

  objsp7.slot("polygons") = Aux7;

  S4 objsp8 = clone(objsp);

  List Aux8 = objsp8.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out8 = Aux8[i];
    List Aux2_8 = out8.slot("Polygons");
    S4 out2_8 = Aux2_8[0];

    NumericMatrix mat_aux = out2_8.slot("coords");

    mat_aux(_,0) = mat_aux(_,0) + range_x;
    mat_aux(_,1) = mat_aux(_,1) - range_y;

    out2_8.slot("coords") = mat_aux;

    NumericVector vec_aux = out8.slot("labpt");

    vec_aux(0) = vec_aux(0) + range_x;
    vec_aux(1) = vec_aux(1) - range_y;

    Aux2_8[0] = out2_8;

    out8.slot("Polygons") = Aux2_8;

    out8.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out8.slot("ID") = ID;

    Aux_Final.push_back(out8);

    Aux8[i] = out8;
  }

  objsp8.slot("polygons") = Aux8;

  S4 objsp9 = clone(objsp);

  List Aux9 = objsp9.slot("polygons");

  for(int i = 0; i < n_poly; i++)
  {
    S4 out9 = Aux9[i];
    List Aux2_9 = out9.slot("Polygons");
    S4 out2_9 = Aux2_9[0];

    NumericMatrix mat_aux = out2_9.slot("coords");

    mat_aux(_,0) = mat_aux(_,0) - range_x;
    mat_aux(_,1) = mat_aux(_,1) + range_y;

    out2_9.slot("coords") = mat_aux;

    NumericVector vec_aux = out9.slot("labpt");

    vec_aux(0) = vec_aux(0) - range_x;
    vec_aux(1) = vec_aux(1) + range_y;

    Aux2_9[0] = out2_9;

    out9.slot("Polygons") = Aux2_9;

    out9.slot("labpt") = vec_aux;

    String ID = (i + 1);

    // ID.push_front("a ");

    out9.slot("ID") = ID;

    Aux_Final.push_back(out9);

    Aux9[i] = out9;
  }

  objsp9.slot("polygons") = Aux9;

  // Defining the Final Output

  NumericMatrix mat_aux = objsp.slot("bbox");
  NumericMatrix box_aux(2,2);

  box_aux = clone(mat_aux);

  box_aux(0,1) = box_aux(0,1) + range_x;
  box_aux(0,0) = box_aux(0,0) - range_x;
  box_aux(1,1) = box_aux(1,1) + range_y;
  box_aux(1,0) = box_aux(1,0) - range_y;

  CharacterVector rNames(2);
  rNames(0) = 'x';
  rNames(1) = 'y';
  rownames(box_aux) = rNames;

  CharacterVector cNames(2);
  cNames(0) = "min";
  cNames(1) = "max";
  colnames(box_aux) = cNames;

  S4 output = clone(objsp);

  S4 projString("CRS");

  String charac_aux = NA_STRING;

  projString.slot("projargs") = charac_aux;

  output.slot("proj4string") = projString;

  output.slot("plotOrder") =  pOrder;

  output.slot("bbox") = box_aux;

  output.slot("polygons") = Aux_Final;

  return output;
}

