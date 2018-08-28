// #include <Rcpp.h>
//
// using namespace Rcpp;
//
// S4 PolyInter(S4 objsp1, S4 objsp2) {
//
//   // Using rpackage rgeos
//
//   Rcpp::Environment base_env("package:rgeos");
//   Rcpp::Function Intersec = base_env["gIntersection"];
//
//   // Cloning R objects
//
//   S4 Poly1("SpatialPolygons");
//   Poly1 = clone(objsp1);
//   S4 Poly2("SpatialPolygons");
//   Poly2 = clone(objsp2);
//
//   // Declaring the output variable
//
//   S4 out = clone(Poly1);
//
//   // List that contais the polygons
//
//   out = Intersec(Poly1, Poly2, true);
//
//   return out;
// }

// NumericMatrix sp_Poly_ID(NumericMatrix mat) {
//
//   // Cloning R objects
//
//   S4 Poly1("SpatialPolygons");
//   Poly1 = clone(objsp1);
//   S4 Poly2("SpatialPolygons");
//   Poly2 = clone(objsp2);
//
//   // Getting polygons IDs
//
//   List Aux1 = Poly1.slot("polygons");
//   List Aux2 = Poly2.slot("polygons");
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
//   for(int i = 0; i < out.ncol(); i++)
//   {
//     for(int j = 0; j < aux.ncol(); j++)
//     {
//       if(ID1_aux[j] == ID1[i])
//       {
//         if(is_true(any(out(_,i) != 0)))
//         {
//           for(int k = 0; k < out.nrow(); k ++)
//           {
//             double o = 0.0, ax = 0.0;
//             o = out(k, i);
//             ax = aux(k, j);
//
//             if(o <= ax)
//             {
//               out(k, i) = o;
//             }
//
//             else
//             {
//               out(k, i) = ax;
//             }
//           }
//         }
//
//         else
//         {
//           out(_, i) = aux(_, j);
//         }
//
//       }
//     }
//   }
//
//   return out;
// }


// // [[Rcpp::export]]
// NumericVector sp_poly_mc_C(S4 objsp1, S4 objsp2, NumericMatrix bbox, int nsim) {
//
//   S4 Poly1 = clone(objsp1);
//   S4 Poly2 = clone(objsp2);
//   NumericMatrix BBox = clone(bbox);
//
//   Rcpp::Environment env = Environment::global_env();
//   Rcpp::Function limits_to_sp = env["limits_to_sp"];
//   Rcpp::Function pam = env["pam"];
//
//   Poly1.slot("bbox") = BBox;
//
//   Poly2.slot("bbox") = BBox;
//
//   S4 obj1_shift = sp_poly_shift(Poly1);
//
//   NumericVector out(nsim);
//
//   for(int i = 0; i < nsim; i++)
//   {
//     S4 obj2_rshift = sp_poly_random_shift(Poly2, obj1_shift.slot("bbox"));
//
//     S4 obj1_aux = PolyInter(obj1_shift, limits_to_sp(obj2_rshift.slot("bbox")));
//
//     obj1_aux.slot("bbox") = obj2_rshift.slot("bbox");
//
//     out[i] = as<double>(pam(obj1_aux, obj2_rshift));
//   }
//
//   return out;
// }
