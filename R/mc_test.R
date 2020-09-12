#' Polygons Spatial Association Test - PSAM
#'
#' @description A Monte Carlo test based on the test statistic \code{\link{psam}}, to verify if
#' two sets of polygons are associated.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param n_sim an \code{integer} corresponding to the number of Monte Carlo simulations
#'  for the test
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains both sets
#' @param alpha a \code{numeric} indicating the confidence level
#' @param alternative a \code{character} indicating the alternative hypothesis, it can be: "two_sided",
#' "repulsion", or "attraction" if you interest is  only check if the sets are independent
#'           or not, if the two sets repulses each other, or if the two sets attracts each other,
#'           respectively.
#' @param fixed a \code{boolean} indicating if the first pattern should be fixed on the toroidal
#' shift or the first will be fixed in half of iterations and then the other one. \code{TRUE} or
#' \code{FALSE}, respectively.
#' @param hausdorff a \code{boolean}. If \code{TRUE}, then the Hausdorff distance is used.
#' Otherwise the Euclidean distance is used. Default is \code{FALSE}.
#' @param ... parameters for test statistics functions
#'
#' @importFrom rgeos gIntersection
#' @importFrom methods slot
#' @importFrom stats quantile
#' @import sp
#'
#' @return a list from class \code{\link{psam_test}}, with values: \describe{
#'     \item{p_value}{a \code{numeric} scalar giving the p-value of the test}
#'     \item{mc_sample}{a \code{numeric} vector giving the test statistic for each of the Monte Carlo simulations}
#'     \item{alternative}{a \code{character} giving the alternative hypothesis}
#'     \item{alpha}{a \code{numeric} scalar giving the significance level}
#'     \item{rejects}{a \code{logical} scalar, TRUE if the null hypothesis is reject}
#'   }
#'
#' @export
#'
psam_mc <- function(obj_sp1, obj_sp2, n_sim = 499L,
                    unique_bbox = NULL, alpha = 0.05,
                    alternative = "two_sided", 
                    fixed = FALSE, hausdorff = F,
                    ...) {

    if((! "SpatialPolygons" %in% class(obj_sp1)) | (! "SpatialPolygons" %in% class(obj_sp2)))
        stop("obj_sp1 and obj_sp2 must be from class 'SpatialPolygons'")

    if(length(alternative) > 1)
        stop("Provide just one alternative.")

    if(! alternative %in% c('two_sided', 'attraction', 'repulsion'))
        stop("Alternative must be 'two_sided', 'attraction' or 'repulsion'.")

    if(length(alpha) > 1 | length(n_sim) > 1)
        stop('alpha and n_sim must be scalars.')

    if(!(alpha > 0 & alpha < 1))
        stop('alpha must lie between 0 and 1.')

    if((!is.null(unique_bbox)) & (!is.matrix(unique_bbox)))
        stop('unique_bbox must be NULL or matrix.')

    if(is.null(unique_bbox)) {
        bbox_1 <- obj_sp1@bbox
        bbox_2 <- obj_sp2@bbox
        unique_bbox <- matrix(c(min(c(bbox_1[1, 1], bbox_2[1, 1])),
                                max(c(bbox_1[1, 2], bbox_2[1, 2])),
                                min(c(bbox_1[2, 1], bbox_2[2, 1])),
                                max(c(bbox_1[2, 2], bbox_2[2, 2]))),
                              ncol = 2, byrow = T)
        rm(bbox_1, bbox_2)
    }

    obj1_shift <- poly_shift(obj_sp = obj_sp1, bbox_tot = unique_bbox)
    if(!fixed) {
        obj2_shift <- poly_shift(obj_sp = obj_sp2, bbox_tot = unique_bbox)
    }

    output <- vector(mode = "list", length = 5L)

    names(output) <- c("p_value", "rejects",
                       "mc_sample", "alternative",
                       "alpha")

    output$mc_sample <- vector(mode = "numeric", length = n_sim + 1L)

    output$alternative <- alternative
    output$alpha <- alpha
    output$rejects <- FALSE

    output$mc_sample[n_sim + 1L] <- psam(obj_sp1, obj_sp2,
                                         hausdorff = hausdorff)

    if(fixed) {
        for(i in seq_len(n_sim)) {
            obj2_rshift <- poly_rf2(obj_sp2)
            lm_sp <- limits_to_sp(obj2_rshift@bbox)
            lm_sp@proj4string <- obj2_rshift@proj4string
            obj1_aux <- rgeos::gIntersection(spgeom1 = obj1_shift,
                                             spgeom2 = lm_sp,
                                             byid = T)
            output$mc_sample[i] <- psam(obj1_aux, obj2_rshift,
                                        hausdorff = hausdorff,
                                        ...)
        }
    } else {
        index_1 <- seq.int(1L, n_sim, 2L)
        index_2 <- seq.int(2L, n_sim, 2L)

        for(i in index_1) {
            obj2_rshift <- poly_rf2(obj_sp2)
            lm_sp <- limits_to_sp(obj2_rshift@bbox)
            lm_sp@proj4string <- obj2_rshift@proj4string
            obj1_aux <- gIntersection(spgeom1 = obj1_shift,
                                      spgeom2 = lm_sp,
                                      byid = T)
            output$mc_sample[i] <- psam(obj1_aux, obj2_rshift,
                                        hausdorff = hausdorff,
                                        ...)
        }

        for(i in index_2) {
            obj1_rshift <- poly_rf2(obj_sp1)
            if(hausdorff) {
                obj2_aux <- poly_touch(obj2_shift, obj1_rshift@bbox)
            } else {
                lm_sp <- limits_to_sp(obj1_rshift@bbox)
                lm_sp@proj4string <- obj1_rshift@proj4string
                obj2_aux <- gIntersection(spgeom1 = obj2_shift,
                                          spgeom2 = lm_sp,
                                          byid = T)
            }
            output$mc_sample[i] <- psam(obj1_rshift, obj2_aux,
                                        hausdorff = hausdorff,
                                        ...)
        }
    }


    if(alternative == "two_sided") {
        output$p_value <- min(mean(output$mc_sample[n_sim + 1L] <= output$mc_sample),
                              mean(output$mc_sample[n_sim + 1L] >= output$mc_sample))
    }
    if(alternative == "attraction") {
        output$p_value <- mean(output$mc_sample[n_sim + 1L] >= output$mc_sample)
    }
    if(alternative == "repulsion") {
        output$p_value <- mean(output$mc_sample[n_sim + 1L] <= output$mc_sample)
    }

    if(output$p_value <= output$alpha) output$rejects <- TRUE

    class(output) <- mc_psa(output)
    class(output) <- psam_test(output)

    rm(list = ls()[!ls() %in% c("output")])

    return(output)
}

#' Polygons Spatial Association Test - Global Envelope
#'
#' @description A Monte Carlo test to verify if two sets of polygons are
#'    associated based in a global envelope of the functions
#'    \eqn{K_{12}(d)} and \eqn{L_{12}(d)} using different test statistics.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param n_sim an \code{integer} corresponding to the number of Monte Carlo simulations
#'  for the test
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains both sets
#' @param alpha a \code{numeric} indicating the confidence level.
#' @param H a \code{character} indicating the function to be used. Possible entries are:
#' \code{'K'} or \code{'L'}.
#' @param ts a \code{character} associated to a test statistic. Inputs acepted:
#' \code{c('IM', 'MAD', 'SIM', 'SMAD', 'IMDQ', 'MADDQ')}.
#' @param distances a \code{numeric vector} indicating the distances to evaluate \eqn{H(d)}. If
#' \code{NULL} then the range considered goes from 5% to 20% of the max distance that can be
#' observed inside the \code{unique_bbox}.
#' @param fixed a \code{boolean} indicating if the first pattern should be fixed on the toroidal
#' shift or the first will be fixed in half of iterations and then the other one. \code{TRUE} or
#' \code{FALSE}, respectively.
#' @param method a \code{character} specifying which kind of distance will be used
#' to evalueate \eqn{H}. Also, there is an option of using areas. Options available:
#' \code{c('hausdorff', 'euclidean', 'area')}.
#'
#' @importFrom rgeos gIntersection
#' @importFrom methods slot
#' @importFrom stats quantile
#' @import sp
#'
#' @return a list from class \code{\link{gof_test}}, with values: \describe{
#'     \item{p_value}{a \code{numeric} scalar giving the p-value of the test}
#'     \item{mc_sample}{a \code{numeric} vector giving the test statistic for each of the Monte Carlo simulations}
#'     \item{mc_funct}{a \code{matrix} where each line correspond to the function (\eqn{K} or \eqn{L}) estimated
#'     for the Monte Carlo simulations}
#'     \item{distances}{\code{numeric vector} containing the distances where mc_func were evaluated.}
#'     \item{alpha}{a \code{numeric} scalar giving the significance level}
#'     \item{rejects}{a \code{logical} scalar, TRUE if the null hypothesis is reject}
#'   }
#'
#' @export
#'
gof_mc <- function(obj_sp1, obj_sp2, n_sim = 499L,
                   unique_bbox = NULL, alpha = 0.01,
                   H = 'L', ts = 'SMAD', distances = NULL,
                   fixed = FALSE, method = 'hausdorff') {

    if((! "SpatialPolygons" %in% class(obj_sp1)) | (! "SpatialPolygons" %in% class(obj_sp2)))
        stop("obj_sp1 and obj_sp2 must be from class 'SpatialPolygons'")

    if(! H %in% c('K', 'L'))
        stop("H must be 'K' or 'L'.")

    if(! method %in% c('euclidean', 'hausdorff', 'area'))
        stop("method must be 'euclidean', 'hausdorff' or 'area'.")

    if(! ts %in% c('IM', 'MAD', 'SIM', 'SMAD', 'IMDQ', 'MADDQ'))
        stop("ts must be 'IM', 'MAD', 'SIM', 'SMAD', 'IMDQ' or 'MADDQ'.")

    if(length(alpha) > 1 | length(n_sim) > 1)
        stop('alpha and n_sim must be scalars.')

    if(length(ts) > 1 | length(H) > 1 | length(method) > 1)
        stop('ts, H, and method should have length == 1L.')

    if(!(alpha > 0 & alpha < 1))
        stop('alpha must lie between 0 and 1.')

    if((!is.null(unique_bbox)) & (!is.matrix(unique_bbox)))
        stop('unique_bbox must be NULL or matrix.')

    if(is.null(unique_bbox)) {
        bbox_1 <- obj_sp1@bbox
        bbox_2 <- obj_sp2@bbox
        unique_bbox <- matrix(c(min(c(bbox_1[1, 1], bbox_2[1, 1])),
                                max(c(bbox_1[1, 2], bbox_2[1, 2])),
                                min(c(bbox_1[2, 1], bbox_2[2, 1])),
                                max(c(bbox_1[2, 2], bbox_2[2, 2]))),
                              ncol = 2, byrow = T)
        rm(bbox_1, bbox_2)
    }

    if(is.null(distances)) {
        max_d <- sqrt(diff(unique_bbox[1, ])^2 + diff(unique_bbox[2, ])^2)
        distances <- seq(from = .05*max_d, to = .2*max_d,
                         length.out = 15L)
    }
    
    obj1_shift <- poly_shift(obj_sp = obj_sp1, bbox_tot = unique_bbox)
    if(!fixed) {
        obj2_shift <- poly_shift(obj_sp = obj_sp2, bbox_tot = unique_bbox)
    }

    output <- vector(mode = "list", length = 6L)

    names(output) <- c("p_value", "rejects",
                       "mc_sample", "mc_funct",
                       "distances", "alpha")

    output$mc_sample <- vector(mode = "numeric", length = n_sim + 1L)
    output$mc_funct <- matrix(ncol = length(distances), nrow = n_sim + 1L)

    output$alpha <- alpha
    output$rejects <- FALSE
    output$distances <- distances

    output$mc_funct[n_sim + 1L, ] <- h_func(obj_sp1, obj_sp2, unique_bbox,
                                            distances, method, H)

    if(fixed) {
        for(i in seq_len(n_sim)) {
            obj2_rshift <- poly_rf2(obj_sp2)
            lm_sp <- limits_to_sp(obj2_rshift@bbox)
            lm_sp@proj4string <- obj2_rshift@proj4string
            obj1_aux <- rgeos::gIntersection(spgeom1 = obj1_shift,
                                             spgeom2 = lm_sp,
                                             byid = T)
            output$mc_funct[i, ] <- h_func(obj1_aux, obj2_rshift,
                                           unique_bbox, distances,
                                           method, H)
        }
    } else {
        index_1 <- seq.int(1L, n_sim, 2L)
        index_2 <- seq.int(2L, n_sim, 2L)

        for(i in index_1) {
            obj2_rshift <- poly_rf2(obj_sp2)
            lm_sp <- limits_to_sp(obj2_rshift@bbox)
            lm_sp@proj4string <- obj2_rshift@proj4string
            obj1_aux <- gIntersection(spgeom1 = obj1_shift,
                                      spgeom2 = lm_sp,
                                      byid = T)
            
            output$mc_funct[i, ] <- h_func(obj1_aux, obj2_rshift,
                                           unique_bbox, distances,
                                           method, H)
        }

        for(i in index_2) {
            obj1_rshift <- poly_rf2(obj_sp1)
            lm_sp <- limits_to_sp(obj1_rshift@bbox)
            lm_sp@proj4string <- obj1_rshift@proj4string
            obj2_aux <- gIntersection(spgeom1 = obj2_shift,
                                      spgeom2 = lm_sp,
                                      byid = T)
            output$mc_funct[i, ] <- h_func(obj1_rshift, obj2_aux,
                                           unique_bbox, distances,
                                           method, H)
        }
    }

    if(length(distances > 1)) {
        h <- distances[length(distances)] - distances[length(distances) - 1L]
    } else {
        h <- 1
    }

    output$mc_sample <- switch(ts,
                               'IM'    = {im(x = output$mc_funct, h = h)},
                               'MAD'   = {mad(x = output$mc_funct)},
                               'SIM'   = {s_im(x = output$mc_funct, h = h)},
                               'SMAD'  = {s_mad(x = output$mc_funct)},
                               'IMDQ'  = {im_ac(x = output$mc_funct, h = h)},
                               'MADDQ' = {mad_ac(x = output$mc_funct)})

    output$p_value <- mean(output$mc_sample[n_sim + 1L] <= output$mc_sample)

    if(output$p_value <= output$alpha) output$rejects <- TRUE

    class(output) <- mc_psa(output)
    class(output) <- gof_test(output)

    rm(list = ls()[!ls() %in% c("output")])

    return(output)
}

#' @useDynLib tpsa
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL
