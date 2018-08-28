#include <RcppArmadillo.h>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// //' Poly dist C
// //'
// //' @param objsp1 object from class \code{SpatialPolygons}
// //' @param objsp2 object from class \code{SpatialPolygons}
// //'
// //' @return an object from class \code{SpatialPolygons} randomly translated
// //'
// // [[Rcpp::export]]
// List PolyDist(S4 objsp1, S4 objsp2) {
//
//   // Using rpackage rgeos
//
//   Rcpp::Environment rgeos_env = Environment::namespace_env("rgeos");
//   Rcpp::Function Dist = rgeos_env["gDistance"];
//
//   List Aux1 = objsp1.slot("polygons");
//   List Aux2 = objsp2.slot("polygons");
//
//   StringVector ID1(Aux1.size());
//   StringVector ID2(Aux2.size());
//
//   for(int i = 0; i < Aux1.size(); i++)
//   {
//     S4 Aux_ = Aux1[i];
//     String aux_ = Aux_.slot("ID");
//     ID1[i] = aux_;
//   }
//
//   for(int i = 0; i < Aux1.size(); i++)
//   {
//     S4 Aux_ = Aux2[i];
//     String aux_ = Aux_.slot("ID");
//     ID2[i] = aux_;
//   }
//
//   // Removing Duplicate IDs
//
//   ID1 = unique(ID1);
//   ID1.sort();
//   ID2 = unique(ID2);
//   ID2.sort();
//
//   // Declaring the output and the auxiliar matrices
//
//   NumericMatrix aux, out(ID2.size(), ID1.size());
//
//   // Calculate the auxiliar distance Matrix
//
//   aux = Dist(objsp1, objsp2, true);
//
//   // rownames and colnames
//
//   colnames(out) = ID1;
//   rownames(out) = ID2;
//
//   StringVector ID1_aux = colnames(aux);
//   StringVector ID2_aux = rownames(aux);
//
//   // Calculate the output Matrix
//
//   // for(int i = 0; i < out.ncol(); i++)
//   // {
//   //   for(int j = 0; j < aux.ncol(); j++)
//   //   {
//   //     if(ID1_aux[j] == ID1[i])
//   //     {
//   //       if(is_true(any(out(_,i) != 0)))
//   //       {
//   //         for(int k = 0; k < out.nrow(); k ++)
//   //         {
//   //           double o, ax;
//   //           o = out(k, i);
//   //           ax = aux(k, j);
//   //
//   //           if(o <= ax)
//   //           {
//   //             out(k, i) = o;
//   //           }
//   //
//   //           else
//   //           {
//   //             out(k, i) = ax;
//   //           }
//   //         }
//   //       }
//   //
//   //       else
//   //       {
//   //         out(_, i) = aux(_, j);
//   //       }
//   //
//   //     }
//   //   }
//   // }
//
//   arma::mat aux2 = as<arma::mat>(aux);
//   arma::mat out2 = as<arma::mat>(out);
//
//
//   List output = List::create(Named("aux") = aux2, Named("out") = out2);
//
//   return output;
// }

template <typename T>
void remove_duplicates(std::vector<T>& vec)
{
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

//' Poly dist3 C
//'
//' @param objsp1 object from class \code{SpatialPolygons}
//' @param objsp2 object from class \code{SpatialPolygons}
//'
//' @return an object from class \code{SpatialPolygons} randomly translated
//'
// [[Rcpp::export]]
List PolyDist(S4& objsp1, S4& objsp2) {

  // Using rpackage rgeos

  Rcpp::Environment rgeos_env = Environment::namespace_env("rgeos");
  Rcpp::Function Dist = rgeos_env["gDistance"];

  List Aux1 = objsp1.slot("polygons");
  List Aux2 = objsp2.slot("polygons");

  std::vector<std::string> ID1(Aux1.size());
  std::vector<std::string> ID2(Aux2.size());

  for(int i = 0; i < Aux1.size(); i++)
  {
    S4 Aux_ = Aux1[i];
    char aux_ = Aux_.slot("ID");
    ID1[i] = aux_;
  }

  for(int i = 0; i < Aux2.size(); i++)
  {
    S4 Aux_ = Aux2[i];
    char aux_ = Aux_.slot("ID");
    ID2[i] = aux_;
  }

  // Removing Duplicate IDs

  std::vector<std::string> ID1_aux = ID1;
  remove_duplicates(ID1_aux);
  std::vector<std::string> ID2_aux = ID2;
  remove_duplicates(ID2_aux);

  // Declaring the output and the auxiliar matrices

  arma::mat aux, out(ID2_aux.size(), ID1_aux.size());

  // Calculate the auxiliar distance Matrix

  aux = as<arma::mat>(Dist(objsp1, objsp2, true));

  if(aux.n_cols == out.n_cols || aux.n_rows == out.n_rows) {
    out = aux;
  }
  // else if (aux.n_cols != out.n_cols && aux.n_rows == out.n_rows) {
  //   for(int i = 0; i < out.n_cols; i++) {
  //     arma::mat aux_col = aux.cols(arma::find(std::find(ID1.begin(), ID1.end(), ID1_aux[i])));
  //     arma::vec min_aux = arma::zeros<arma::vec>(aux_col.n_rows);
  //     if(aux_col.n_cols > 1) {
  //       for(int j = 0; j < min_aux.size(); j++) {
  //         arma::vec min_aux2 = aux_col.row(j);
  //         min_aux(j) = min_aux2.min();
  //       }
  //     }
  //     out.col(i) = min_aux;
  //   }
  // }
  // else if (aux.n_cols == out.n_cols && aux.n_rows != out.n_rows) {
  //   for(int i = 0; i < out.n_rows; i++) {
  //     arma::mat aux_row = aux.rows(arma::find(std::find(ID2.begin(), ID2.end(), ID2_aux[i])));
  //     arma::vec min_aux = arma::zeros<arma::vec>(aux_row.n_rows);
  //     if(aux_row.n_rows > 1) {
  //       for(int j = 0; j < min_aux.size(); j++) {
  //         arma::vec min_aux2 = aux_row.col(j);
  //         min_aux(j) = min_aux2.min();
  //       }
  //     }
  //     out.row(i) = min_aux;
  //   }
  // }

  List output = List::create(Named("aux") = aux, Named("out") = out,
                             Named("rownames_aux") = ID2, Named("rownames_out") = ID2_aux,
                             Named("colnames_aux") = ID1, Named("colnames_out") = ID1_aux);

  return output;
}

