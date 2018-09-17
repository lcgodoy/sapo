#' Polygons Spatial Association Test
#'
#' @description A Monte Carlo test to verify if two sets of polygons are
#'    associated. The difference bewtween this function and the function
#'    \code{psat_mc} is that here you can choose the test statistic.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param n_sim an \code{integer} corresponding to the number of Monte Carlo simulations
#'  for the test
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains the both sets
#' @param alpha a \code{numeric} indicating the confidence level
#' @param ts a \code{character} indicating the test statistic used in the test, the options are
#'  \code{c('psam', 'pk_dist12', 'pk_area12', 'pf12')}.
#' @param alternative a \code{character} indicating the alternative hypothesis, it can be: "two_sided",
#' "repulsion", or "attraction" if you interest is  only check if the sets are independent
#'           or not, if the two sets repulses each other, or if the two sets attracts each other,
#'           respectively.
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param fixed a \code{boolean} indicating if the first pattern should be fixed on the toroidal
#' shift or the first will be fixed in half of iterations and then the other one. \code{TRUE} or
#' \code{FALSE}, respectively.
#' @param r_min min distance to calculate \eqn{K_{1,2}} functions. Used only if ts != 'psam'
#' @param r_max max distance to calculate \eqn{K_{1,2}} functions. Used only if ts != 'psam'
#' @param by controls how many values between \code{r_min} and \code{r_max} will by used
#' to calculate \eqn{K_{1,2}} functions. Used only if ts != 'psam'
#' @param ... parameters for test statistics functions
#'
#' @importFrom rgeos gIntersection
#' @importFrom methods slot
#' @importFrom stats quantile
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
psat_mc <- function(obj_sp1, obj_sp2, n_sim = 500L,
                    unique_bbox = NULL,
                    alpha = 0.01, ts = 'psam',
                    alternative = "two_sided",
                    correction = 'none', fixed = FALSE,
                    r_min = NULL, r_max = NULL,
                    by = NULL,
                    ...) {

  if((! "SpatialPolygons" %in% class(obj_sp1)) | (! "SpatialPolygons" %in% class(obj_sp2)))
    stop("obj_sp1 and obj_sp2 must be from class 'SpatialPolygons'")

  if(length(alternative) > 1)
    stop("Provide just one alternative.")

  if(length(ts) > 1)
    stop("Provide just one ts")

  if(length(correction) > 1)
    stop("Provide just one correction")

  if(! alternative %in% c('two_sided', 'attraction', 'repulsion'))
    stop("Alternative must be 'two_sided', 'attraction' or 'repulsion'.")

  if(! ts %in% c('psam', 'pk_dist12', 'pk_area12'))
    stop("alternative must be 'psam', 'pk_dist12' or 'pk_area12'.")

  if(ts == 'psam' & correction == 'adjust') {
    warning('This edge correction is not available for psam, none will be used.')
    correction <- 'none'
  }

  if(! correction %in% c('none', 'guard', 'torus', 'adjust'))
    stop("correction must be 'none', 'guard', 'torus' or 'adjust'.")

  if(length(alpha) > 1 | length(n_sim) > 1)
    stop('alpha and n_sim must be scalars.')

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

  if(ts != 'psam') {
    if(is.null(r_max) | is.null(r_min)) {
      r_x <- unique_bbox[1, 2] - unique_bbox[1, 1]
      r_y <- unique_bbox[2, 2] - unique_bbox[2, 1]
    }

    if(is.null(r_max)) {
      r_max <- .15*max(r_x, r_y)
    }

    if(is.null(r_min)) {
      r_min <- .05*max(r_x, r_y)
    }

    if(is.null(by)) {
      by <- (r_max - r_min)/sqrt(12)
    }
  }

  if(correction == 'torus') {
    obj1_shift <- torus_corr(objsp = obj_sp1, bbox_tot = unique_bbox)
    if(!fixed) {
      obj2_shift <- torus_corr(objsp = obj_sp2, bbox_tot = unique_bbox)
    }
  } else {
    obj1_shift <- poly_shift(obj_sp = obj_sp1, bbox_tot = unique_bbox)
    if(!fixed) {
      obj2_shift <- poly_shift(obj_sp = obj_sp2, bbox_tot = unique_bbox)
    }
  }

  output <- vector(mode = "list", length = 6)

  names(output) <- c("p_value", "rejects",
                     "sample_ts", "mc_ts",
                     "alternative", "alpha")

  # mc_values <- vector(mode = "numeric", length = n_sim)

  output$alternative <- alternative
  output$alpha <- alpha
  # output$rejects <- FALSE

  if(ts == 'psam') {

    output$rejects <- FALSE

    output$sample_ts <- psam(obj_sp1, obj_sp2, correction = correction)

    output$mc_ts <- vector(mode = 'numeric')

    if(fixed) {
      output$mc_ts <- mc_iterations(obj1_shift, obj_sp2,
                                    niter = n_sim, ts = ts,
                                    correction = correction, ...)
    } else {
      output$mc_ts <- c(mc_iterations(obj1_shift, obj_sp2,
                                      niter = length(1:(n_sim/2)), ts = ts,
                                      correction = correction, ...),
                        mc_iterations(obj2_shift, obj_sp1,
                                      niter = length(trunc(n_sim/2 + 1):n_sim), ts = ts,
                                      correction = correction, ...))
    }

    if(alternative == "two_sided") {
      output$p_value <- min(mean(output$sample_ts <= output$mc_ts), mean(output$sample_ts >= output$mc_ts))
    }
    if(alternative == "attraction") {
      output$p_value <- mean(output$sample_ts >= output$mc_ts)
    }
    if(alternative == "repulsion") {
      output$p_value <- mean(output$sample_ts <= output$mc_ts)
    }

    if(output$p_value <= output$alpha) output$rejects <- TRUE

    class(output) <- psa_psam(output)
  }

  if(ts == 'pk_dist12') {

    output$sample_ts <- pk_dist12(obj_sp1, obj_sp2, bbox = obj_sp1@bbox,
                                  correction = correction, r_min = r_min,
                                  r_max = r_max, by = by, ...)

    if(fixed) {
      mc_aux <-mc_iterations(obj1_shift, obj_sp2,
                             niter = n_sim, ts = ts,
                             correction = correction,
                             r_min = r_min, r_max = r_max,
                             by = by, ...)
    } else {
      mc_aux <- data.table::rbindlist(list(mc_iterations(obj1_shift, obj_sp2,
                                                         niter = length(1:(n_sim/2)), ts = ts,
                                                         correction = correction,
                                                         r_min = r_min, r_max = r_max,
                                                         by = by, ...),
                                           mc_iterations(obj2_shift, obj_sp1,
                                                         niter = length(trunc(n_sim/2 + 1):n_sim), ts = ts,
                                                         correction = correction,
                                                         r_min = r_min, r_max = r_max,
                                                         by = by, ...)),
                                      use.names = F, fill = F, idcol = NULL)
    }

    p_value <- vector(mode = 'numeric', length = nrow(output$sample_ts))

    if(alternative == "two_sided") {
      output$mc_ts <- data.frame(r = unique(mc_aux$r),
                                 k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = (alpha/2)),
                                 k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - (alpha/2))
      )
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- mc_aux[mc_aux$r == output$sample_ts[i, 1], ]

        p_value[i] <- min(mean(aux$pk12 >= output$sample_ts[i, 2]), mean(aux$pk12 <= output$sample_ts[i, 2]))

        rm(aux)
      }
      output$p_value <- p_value
    }

    if(alternative == "attraction") {
      output$mc_ts <- data.frame(r = unique(mc_aux$r),
                                 k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = alpha),
                                 k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - alpha)
      )
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- mc_aux[mc_aux$r == output$sample_ts[i, 1], ]

        p_value[i] <- mean(aux$pk12 >= output$sample_ts[i, 2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    if(alternative == "repulsion") {
      output$mc_ts <- data.frame(r = unique(mc_aux$r),
                                 k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = alpha),
                                 k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - alpha)
      )
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- mc_aux[mc_aux$r == output$sample_ts[i, 1], ]

        p_value[i] <- mean(aux$pk12 <= output$sample_ts[i, 2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    output$rejects <- (min(output$p_value, na.rm = T) < output$alpha)

    class(output) <- psa_pk12(output)
  }

  if(ts == 'pk_area12') {

    output$sample_ts <- pk_area12(obj_sp1, obj_sp2, bbox = obj_sp1@bbox,
                                  correction = correction,
                                  r_min = r_min, r_max = r_max,
                                  by = by, ...)

    if(fixed) {
      mc_aux <-mc_iterations(obj1_shift, obj_sp2,
                             niter = n_sim, ts = ts,
                             correction = correction,
                             r_min = r_min, r_max = r_max,
                             by = by, ...)
    } else {
      mc_aux <- data.table::rbindlist(list(mc_iterations(obj1_shift, obj_sp2,
                                                         niter = length(1:(n_sim/2)), ts = ts,
                                                         correction = correction,
                                                         r_min = r_min, r_max = r_max,
                                                         by = by, ...),
                                           mc_iterations(obj2_shift, obj_sp1,
                                                         niter = length(trunc(n_sim/2 + 1):n_sim), ts = ts,
                                                         correction = correction,
                                                         r_min = r_min, r_max = r_max,
                                                         by = by, ...)),
                                      use.names = F, fill = F, idcol = NULL)
    }

    p_value <- vector(mode = 'numeric', length = nrow(output$sample_ts))

    if(alternative == "two_sided") {
      output$mc_ts <- data.frame(r = unique(mc_aux$r),
                                 k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = (alpha/2)),
                                 k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - (alpha/2))
      )
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- mc_aux[mc_aux$r == output$sample_ts[i, 1], ]

        p_value[i] <- min(mean(aux$pk12 >= output$sample_ts[i, 2]), mean(aux$pk12 <= output$sample_ts[i, 2]))

        rm(aux)
      }
      output$p_value <- p_value
    }

    if(alternative == "attraction") {
      output$mc_ts <- data.frame(r = unique(mc_aux$r),
                                 k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = alpha),
                                 k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - alpha)
      )
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- mc_aux[mc_aux$r == output$sample_ts[i, 1], ]

        p_value[i] <- mean(aux$pk12 >= output$sample_ts[i, 2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    if(alternative == "repulsion") {
      output$mc_ts <- data.frame(r = unique(mc_aux$r),
                                 k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = alpha),
                                 k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - alpha)
      )
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- mc_aux[mc_aux$r == output$sample_ts[i, 1], ]

        p_value[i] <- mean(aux$pk12 <= output$sample_ts[i, 2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    output$rejects <- min(output$p_value, na.rm = T) < output$alpha

    class(output) <- psa_pk12(output)
  }

  class(output) <- psa_test(output)

  rm(list = ls()[!ls() %in% c("output")])

  return(output)
}

#' @useDynLib tpsa
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL
