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
#' @param same_bbox a \code{boolean} - (desnecessario)
#' @param bbox_1 a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  of the first spatial object
#' @param bbox_2 a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  of the second spatial object
#' @param alpha a \code{numeric} indicating the confidence level
#' @param ts a \code{character} indicating the test statistic used in the test, the options are
#'  \code{c('psam', 'pk_dist12', 'pk_area12', 'pf12')}.
#' @param alternative a \code{character} indicating the alternative hypothesis, it can be: "independece",
#' "repulsion", or "attraction" if you interest is  only check if the sets are independent
#'           or not, if the two sets repulses each other, or if the two sets attracts each other,
#'           respectively.
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
psat_mc <- function(obj_sp1, obj_sp2, n_sim = 500L, unique_bbox = NULL,
                     same_bbox = T, bbox_1 = NULL, bbox_2 = NULL,
                     alpha = 0.01, ts = 'psam',
                     alternative = "two_sided") {

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
    k <- which(rgeos::gWithin(obj_sp1, limits_to_sp(unique_bbox), byid = T))
    obj_sp1 <- obj_sp1[k, ]
    rm(k)
  }

  if(class(obj_sp2) %in% 'SpatialPolygons') {
    obj_sp2 <- gIntersection(obj_sp2, limits_to_sp(unique_bbox), byid = T,
                             id = suppressWarnings(names(obj_sp2)))
  } else {
    k <- which(rgeos::gWithin(obj_sp2, limits_to_sp(unique_bbox), byid = T))
    obj_sp2 <- obj_sp2[k, ]
    rm(k)
  }

  # if(ts == 'pk_area12') {
  #   if(sp::is.projected(obj_sp1)) obj_sp1 <- rgeos::gBuffer(obj_sp1, byid = T, width = 0)
  #   if(sp::is.projected(obj_sp2)) obj_sp2 <- rgeos::gBuffer(obj_sp2, byid = T, width = 0)
  # }

  attr(obj_sp1, "bbox") <- unique_bbox

  attr(obj_sp2, "bbox") <- unique_bbox

  obj1_shift <- poly_shift(obj_sp = obj_sp1, bbox_tot = unique_bbox)
  obj2_shift <- poly_shift(obj_sp = obj_sp2, bbox_tot = unique_bbox)

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

    output$sample_ts <- psam(obj_sp1, obj_sp2)

    output$mc_ts <- vector(mode = 'numeric')

    output$mc_ts <- c(mc_iterations(obj1_shift, obj_sp2,
                                    niter = round((n_sim/2) + .5), ts = ts),
                      mc_iterations(obj2_shift, obj_sp1,
                                    niter = round((n_sim/2)), ts = ts))

    if(alternative == "two_sided") {
      output$p_value <- min(mean(output$sample_ts < output$mc_ts), mean(output$sample_ts > output$mc_ts))
    }
    if(alternative == "attraction") {
      output$p_value <- mean(output$sample_ts > output$mc_ts)
    }
    if(alternative == "repulsion") {
      output$p_value <- mean(output$sample_ts < output$mc_ts)
    }

    if(output$p_value <= output$alpha) output$rejects <- TRUE

    class(output) <- psa_psam(output)
  }

  if(ts == 'pf12') {
    output$sample_ts <- as.data.frame(pf12(obj_sp1, obj_sp2))

    # mc_values <- vector(mode = 'numeric')
    r_x <- unique_bbox[1,2] - unique_bbox[1,1]
    r_y <- unique_bbox[2,2] - unique_bbox[2,1]
    rmax <- .4*max(r_x, r_y)
    rm(r_x, r_y)

    mc_aux <- rbind(mc_iterations(obj1_shift, obj_sp2,
                                  niter = round((n_sim/2) + .5), ts = ts,
                                  args = list(r_min = min(output$sample_ts$r), r_max = rmax)),
                    mc_iterations(obj2_shift, obj_sp1,
                                  niter = round((n_sim/2)), ts = ts,
                                  args = list(r_min = min(output$sample_ts$r), r_max = rmax)))

    output$mc_ts <- data.frame(r = unique(mc_aux$r),
                               f12_inf = tapply(mc_aux$pf12, mc_aux$r, quantile, p = (alpha/2)),
                               f12_up = tapply(mc_aux$pf12, mc_aux$r, quantile, p = 1 - (alpha/2)),
                               f21_inf = tapply(mc_aux$pf21, mc_aux$r, quantile, p = (alpha/2)),
                               f21_up = tapply(mc_aux$pf21, mc_aux$r, quantile, p = 1 - (alpha/2)))

    p12 <- vector(mode = 'numeric', length = nrow(output$sample_ts))
    p21 <- vector(mode = 'numeric', length = nrow(output$sample_ts))

    if(alternative == "two_sided") {
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p12[i] <- min(mean(aux$pf12 > output$sample_ts[i,2]), mean(aux$pf12 < output$sample_ts[i,2]))

        p21[i] <- min(mean(aux$pf12 > output$sample_ts[i,3]), mean(aux$pf12 < output$sample_ts[i,3]))

        rm(aux)
      }
      output$p_value <- data.frame(p12 = p12, p21 = p21)
    }

    if(alternative == "attraction") {
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p12[i] <- mean(aux$pf12 < output$sample_ts[i,2])

        p21[i] <- mean(aux$pf12 < output$sample_ts[i,3])

        rm(aux)
      }

      output$p_value <- data.frame(p12 = p12, p21 = p21)
    }

    if(alternative == "repulsion") {
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p12[i] <- mean(aux$pf12 > output$sample_ts[i,2])

        p21[i] <- mean(aux$pf12 > output$sample_ts[i,3])

        rm(aux)
      }

      output$p_value <- data.frame(p12 = p12, p21 = p21)
    }

    output$rejects <- vector(mode = 'logical', length = 2)

    output$rejects[1] <- ifelse(any(output$p_value$p12 <= output$alpha), TRUE, FALSE)
    output$rejects[2] <- ifelse(any(output$p_value$p21 <= output$alpha), TRUE, FALSE)

    names(output$rejects) <- c('f12', 'f21')

    class(output) <- psa_pf12(output)
  }

  if(ts == 'pk_dist12') {
    r_x <- unique_bbox[1,2] - unique_bbox[1,1]
    r_y <- unique_bbox[2,2] - unique_bbox[2,1]
    rmax <- .4*max(r_x, r_y)
    rm(r_x, r_y)

    output$sample_ts <- pk_dist12(obj_sp1, obj_sp2, r_max = rmax, bbox = obj_sp1@bbox)

    mc_aux <- rbind(mc_iterations(obj1_shift, obj_sp2,
                                  niter = round((n_sim/2) + .5), ts = ts,
                                  args = list(r_min = min(output$sample_ts$r), r_max = rmax)),
                    mc_iterations(obj2_shift, obj_sp1,
                                  niter = round((n_sim/2)), ts = ts,
                                  args = list(r_min = min(output$sample_ts$r), r_max = rmax)))

    output$mc_ts <- data.frame(r = unique(mc_aux$r),
                               k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = (alpha/2)),
                               k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - (alpha/2))
                               )

    p_value <- vector(mode = 'numeric', length = nrow(output$sample_ts))

    if(alternative == "two_sided") {
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p_value[i] <- min(mean(aux$pk12 > output$sample_ts[i,2]), mean(aux$pk12 < output$sample_ts[i,2]))

        rm(aux)
      }
      output$p_value <- p_value
    }

    if(alternative == "attraction") {
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p_value[i] <- mean(aux$pk12 < output$sample_ts[i,2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    if(alternative == "repulsion") {
      for(i in seq_len(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p_value[i] <- mean(aux$pk12 > output$sample_ts[i,2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    output$rejects <- ifelse(any(output$p_value <= output$alpha), TRUE, FALSE)

    class(output) <- psa_pk12(output)
  }

  if(ts == 'pk_area12') {
    r_x <- unique_bbox[1,2] - unique_bbox[1,1]
    r_y <- unique_bbox[2,2] - unique_bbox[2,1]
    rmax <- .4*max(r_x, r_y)
    rm(r_x, r_y)

    output$sample_ts <- pk_area12(obj_sp1, obj_sp2, r_max = rmax, bbox = obj_sp1@bbox)

    mc_aux <- rbind(mc_iterations(obj1_shift, obj_sp2,
                                  niter = round((n_sim/2) + .5), ts = ts,
                                  args = list(r_min = min(output$sample_ts$r), r_max = rmax)),
                    mc_iterations(obj2_shift, obj_sp1,
                                  niter = round((n_sim/2)), ts = ts,
                                  args = list(r_min = min(output$sample_ts$r), r_max = rmax)))

    output$mc_ts <- data.frame(r = unique(mc_aux$r),
                               k12_inf = tapply(mc_aux$pk12, mc_aux$r, quantile, p = (alpha/2)),
                               k12_up = tapply(mc_aux$pk12, mc_aux$r, quantile, p = 1 - (alpha/2))
    )

    p_value <- vector(mode = 'numeric', length = nrow(output$sample_ts))

    if(alternative == "two_sided") {
      for(i in seq_along(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p_value[i] <- min(mean(aux$pk12 > output$sample_ts[i,2]), mean(aux$pk12 < output$sample_ts[i,2]))

        rm(aux)
      }
      output$p_value <- p_value
    }

    if(alternative == "attraction") {
      for(i in seq_along(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p_value[i] <- mean(aux$pk12 < output$sample_ts[i,2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    if(alternative == "repulsion") {
      for(i in seq_along(nrow(output$sample_ts))) {
        aux <- subset(mc_aux, mc_aux$r == output$sample_ts[i,1])

        p_value[i] <- mean(aux$pk12 > output$sample_ts[i,2])

        rm(aux)
      }

      output$p_value <- p_value
    }

    output$rejects <- ifelse(any(output$p_value <= output$alpha), TRUE, FALSE)

    class(output) <- psa_pk12(output)
  }

  class(output) <- psa_test(output)

  rm(list = ls()[!ls() %in% c("output")])

  return(output)
}
