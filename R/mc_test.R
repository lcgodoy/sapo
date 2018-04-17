#' Polygons Spatial Association Test (Old version)
#'
#' @description A Monte Carlo test to verify if two sets of polygons are
#'    or not.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param n_sim an \code{integer} corresponding to the number of Monte Carlo simulations
#'  for the test
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains the both sets
#' @param same_bbox a \code{boolean} - (desnecessario)
#' @param bbox_1 a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  of the first spatial object
#' @param bbox_2 a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  of the second spatial object
#' @param alpha a \code{numeric} indicating the confidence level
#' @param alternative a string indicating the alternative hypothesis, it can be: "independece",
#' "repulsion", or "attraction" if you interest is  only check if the sets are independent
#'           or not, if the two sets repulses each other, or if the two sets attracts each other,
#'           respectively.
#'
#' @importFrom rgeos gIntersection
#' @importFrom methods slot
#' @import sp
#'
#' @return a list from class \code{\link{psa_test}}, with values: \describe{
#'     \item{p_value}{a \code{numeric} scalar giving the p-value of the test}
#'     \item{sample_ts}{a \code{numeric} scalar giving the test statistic calculated in the original sample}
#'     \item{mc_ts}{a \code{numeric} vector giving the test statistic for each of the Monte Carlo simulations}
#'     \item{alternative}{a \code{character} giving the alternative hypothesis}
#'     \item{alpha}{a \code{numeric} scalar giving the significance level used on the test}
#'     \item{rejects}{a \code{logical} scalar, TRUE if the null hypothesis is reject}
#'   }
#'
#' @export
#'
psat_mc2 <- function(obj_sp1, obj_sp2, n_sim = 100L, unique_bbox = NULL,
                    same_bbox = T, bbox_1 = NULL, bbox_2 = NULL,
                    alpha = 0.05, alternative = "two_sided") {

  if(!(class(obj_sp1) %in% c("SpatialPolygons",
                             "SpatialPointsDataFrame") &
       class(obj_sp2) %in% c("SpatialPolygons",
                             "SpatialPointsDataFrame"))) {
    stop("obj_sp1 and obj_sp2 must be from class SpatialPolygons or SpatialPointsDataFrame")
  }

  if((class(obj_sp1) %in% "SpatialPointsDataFrame" &
      class(obj_sp2) %in% "SpatialPointsDataFrame")) {
    warning("if obj_sp1 and obj_sp2 are from class SpatialPolygons,
            then an approach for multitype point patterns
            would be more appropiated.")
  }

  if(is.null(bbox_1)) {
    bbox_1 <- obj_sp1@bbox
  }

  if(is.null(bbox_2)) {
    bbox_2 <- obj_sp2@bbox
  }

  if(is.null(unique_bbox) & isTRUE(same_bbox)) {
    unique_bbox <- matrix(c(min(c(bbox_1[1,1], bbox_2[1,1])),
                            max(c(bbox_1[1,2], bbox_2[1,2])),
                            min(c(bbox_1[2,1], bbox_2[2,1])),
                            max(c(bbox_1[2,2], bbox_2[2,2]))),
                          ncol = 2, byrow = T)
  }

  if(!isTRUE(same_bbox)) {
    unique_bbox <- matrix(c(max(c(bbox_1[1,1], bbox_2[1,1])),
                            min(c(bbox_1[1,2], bbox_2[1,2])),
                            max(c(bbox_1[2,1], bbox_2[2,1])),
                            min(c(bbox_1[2,2], bbox_2[2,2]))),
                          ncol = 2, byrow = T)

    rm(bbox_1, bbox_2)
    gc()
  }

  if(!alternative %in% c("attraction", "two_sided", "repulsion")) {
    stop('alternative value must be: attraction, two_sided or repulsion!')
  }

  if(class(obj_sp1) %in% 'SpatialPolygons') {
    obj_sp1 <- gIntersection(obj_sp1, limits_to_sp(unique_bbox), byid = T,
                             id = suppressWarnings(names(obj_sp1)))
  } else {
    # id_aux <- slot(obj_sp1, 'data')
    # obj_sp1 <- gIntersection(obj_sp1, limits_to_sp(unique_bbox),
    #                          byid = T,
    #                          id = as.character(id_aux$id))
    # obj_sp1 <- sp::SpatialPointsDataFrame(obj_sp1, data = id_aux)
    k <- which(rgeos::gWithin(obj_sp1, limits_to_sp(unique_bbox), byid = T))
    obj_sp1 <- obj_sp1[k, ]
    rm(k)
  }


  if(class(obj_sp2) %in% 'SpatialPolygons') {
    obj_sp2 <- gIntersection(obj_sp2, limits_to_sp(unique_bbox), byid = T,
                             id = suppressWarnings(names(obj_sp2)))
  } else {
    # id_aux <- slot(obj_sp2, 'data')
    # obj_sp2 <- gIntersection(obj_sp2, limits_to_sp(unique_bbox),
    #                          byid = T,
    #                          id = as.character(id_aux$id))
    # obj_sp2 <- sp::SpatialPointsDataFrame(obj_sp2, data = id_aux)

    k <- which(rgeos::gWithin(obj_sp2, limits_to_sp(unique_bbox), byid = T))
    obj_sp2 <- obj_sp2[k, ]
    rm(k)
  }

  attr(obj_sp1, "bbox") <- unique_bbox

  attr(obj_sp2, "bbox") <- unique_bbox

  obj1_shift <- poly_shift(obj_sp = obj_sp1, bbox_tot = unique_bbox)
  obj2_shift <- poly_shift(obj_sp = obj_sp2, bbox_tot = unique_bbox)

  output <- vector(mode = "list", length = 6)
  names(output) <- c("p_value", "rejects",
                     "sample_ts", "mc_ts",
                     "alternative", "alpha")

  mc_values <- vector(mode = "numeric", length = n_sim)

  output$sample_ts <- psam(obj_sp1, obj_sp2)
  output$alternative <- alternative
  output$alpha <- alpha
  output$rejects <- FALSE

  for(i in 1:round(n_sim/2)) {
    obj2_rshift <- poly_rf(obj_sp2, obj1_shift@bbox)

    obj1_aux <- obj1_shift

    if(class(obj1_aux) %in% 'SpatialPolygons') {
      obj1_aux <- gIntersection(obj1_aux, limits_to_sp(obj2_rshift@bbox), byid = T,
                                id = suppressWarnings(names(obj1_aux)))
    } else {
      # id_aux <- slot(obj1_aux, 'data')
      # obj1_aux <- gIntersection(obj1_aux, limits_to_sp(obj2_rshift@bbox), byid = T,
      #                           id = as.character(id_aux$id))
      # obj1_aux <- sp::SpatialPointsDataFrame(obj1_aux,
      #                                        data = data.frame(id = row.names(obj1_aux),
      #                                                          row.names = row.names(obj1_aux)))
      k <- which(rgeos::gWithin(obj1_aux, limits_to_sp(obj2_rshift@bbox), byid = T))
      obj1_aux <- obj1_aux[k, ]
      rm(k)
    }

    attr(obj1_aux, "bbox") <- obj2_rshift@bbox

    mc_values[i] <- psam(obj1_aux, obj2_rshift)
  }

  for(i in (round(n_sim/2) + 1):n_sim) {
    obj1_rshift <- poly_rf(obj_sp1, obj2_shift@bbox)

    obj2_aux <- obj2_shift

    if(class(obj2_aux) %in% 'SpatialPolygons') {
      obj2_aux <- gIntersection(obj2_aux, limits_to_sp(obj1_rshift@bbox), byid = T,
                                id = suppressWarnings(names(obj2_aux)))
    } else {
      # id_aux <- slot(obj2_aux, 'data')
      # obj2_aux <- gIntersection(obj2_aux, limits_to_sp(obj1_rshift@bbox), byid = T,
      #                           id = as.character(id_aux$id))
      # obj2_aux <- sp::SpatialPointsDataFrame(obj2_aux,
      #                                        data = data.frame(id = row.names(obj2_aux),
      #                                                          row.names = row.names(obj2_aux)))
      k <- which(rgeos::gWithin(obj2_aux, limits_to_sp(obj1_rshift@bbox), byid = T))
      obj2_aux <- obj2_aux[k, ]
      rm(k)
    }

    attr(obj2_aux, "bbox") <- obj1_rshift@bbox

    mc_values[i] <- psam(obj2_aux, obj1_rshift)
  }

  output$mc_ts <- mc_values

  if(alternative == "two_sided") {
    output$p_value <- mean(output$sample_ts < mc_values | output$sample_ts > mc_values)
  }
  if(alternative == "attraction") {
    output$p_value <- mean(output$sample_ts > mc_values)
  } else {
    output$p_value <- mean(output$sample_ts < mc_values)
  }

  if(output$p_value <= output$alpha) output$rejects <- TRUE

  # class(output) <- "psa_test"
  class(output) <- psa_psam(output)
  class(output) <- psa_test(output)

  rm(list = ls()[!ls() %in% c("output", "mc_values")])

  return(output)

}


#' @useDynLib tpsa
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import magrittr
NULL
