#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
S4 PolyInter(S4& objsp1, S4& objsp2) {

  // Using rpackage rgeos

  Rcpp::Environment base_env("package:rgeos");
  Rcpp::Function Intersec = base_env["gIntersection"];

  // Cloning R objects

  S4 Poly1("SpatialPolygons");
  Poly1 = clone(objsp1);
  S4 Poly2("SpatialPolygons");
  Poly2 = clone(objsp2);

  // Declaring the output variable

  S4 out = clone(Poly1);

  // List that contais the polygons

  out = Intersec(Poly1, Poly2, true);

  return out;
}

