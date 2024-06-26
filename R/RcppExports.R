# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Vector of means
#'
#' @param x \code{numeric vector}
#'
#' @description Auxiliar function to calculate \eqn{\hat{H_{12}}_i(d)}.
#'
#' @return \code{numeric vector}
#'
mean_vec <- function(x) {
    .Call('_sapo_mean_vec', PACKAGE = 'sapo', x)
}

#' Polygons' Random Shift - 2
#'
#' @param objsp object from class \code{SpatialPolygons}
#'
#' @return an object from class \code{SpatialPolygons} randomly translated
#'
poly_rf2 <- function(objsp) {
    .Call('_sapo_poly_rf2', PACKAGE = 'sapo', objsp)
}

#' Create copies of a set of polygons
#'
#' @param obj_sp object from class \code{SpatialPolygons}
#' @param bbox_tot Boundary box from class \code{matrix}
#'
#' @description Auxiliar function.
#'
#' @importFrom stats runif
#' @return an object from class \code{SpatialPolygons}
#'
poly_shift <- function(obj_sp, bbox_tot) {
    .Call('_sapo_poly_shift', PACKAGE = 'sapo', obj_sp, bbox_tot)
}

#' Create copies of a set of polygons
#'
#' @param obj_sp object from class \code{SpatialPolygons}
#' @param bbox_tot Boundary box from class \code{matrix}
#'
#' @description Auxiliar function.
#'
#' @importFrom stats runif
#' @return an object from class \code{SpatialPolygons}
#'
poly_shift_noid <- function(obj_sp, bbox_tot) {
    .Call('_sapo_poly_shift_noid', PACKAGE = 'sapo', obj_sp, bbox_tot)
}

#' Polygons that touch a bbox
#'
#' @param x a \code{SpatialPolygon}
#' @param bbox a \code{numeric matrix}
#'
#' @export
#'
poly_touch <- function(x, bbox) {
    .Call('_sapo_poly_touch', PACKAGE = 'sapo', x, bbox)
}

#' Auxiliar function
#'
#' @param obj_sp a spatial object
#' @param obj_sp2 a spatial object
#' @param obj_sp3 a spatial object
#' @param obj_sp4 a spatial object
#' @param n_poly an integer
#' @param range_x a numeric
#' @param range_y a numeric
#'
shift_aux <- function(obj_sp, obj_sp2, obj_sp3, obj_sp4, n_poly, range_x, range_y) {
    .Call('_sapo_shift_aux', PACKAGE = 'sapo', obj_sp, obj_sp2, obj_sp3, obj_sp4, n_poly, range_x, range_y)
}

#' Toroidal edge correction
#'
#' @param objsp a \code{SpatialPolygon}
#' @param bbox_tot a \code{numeric matrix}
#'
#' @export
#'
torus_corr <- function(objsp, bbox_tot) {
    .Call('_sapo_torus_corr', PACKAGE = 'sapo', objsp, bbox_tot)
}

