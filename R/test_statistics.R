#' Polygons' Spatial Association Measure
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard')}.
#' @param ... parameters associated with edge correction.
#'
#' @return a \code{numeric} \code{scalar} with the polygons' spatial
#'  association measure for the corresponding objects.
#' @export
#'
psam <- function(obj_sp1, obj_sp2, correction = 'none', ...) {

  switch (correction,
          'none' = {
            m_dist <- sp_ID_dist(obj_sp1, obj_sp2, ...)

            min_row <- apply(m_dist, 1, min)
            min_col <- apply(m_dist, 2, min)
          },
          'torus' = {
            obj_sp1_t <- torus_corr(obj_sp1)
            obj_sp2_t <- torus_corr(obj_sp2)
            min_row <- apply(sp_ID_dist(obj_sp1, obj_sp2_t, ...), 1, min)
            min_col <- apply(sp_ID_dist(obj_sp1_t, obj_sp2, ...), 2, min)
          },
          'guard' = {
            new_bbox <- make_guard(obj_sp1@bbox, ...)
            obj_sp1_ng <- rgeos::gIntersection(obj_sp1, limits_to_sp(new_bbox), byid = T)
            obj_sp2_ng <- rgeos::gIntersection(obj_sp2, limits_to_sp(new_bbox), byid = T)
            min_row <- apply(sp_ID_dist(obj_sp1_ng, obj_sp2, ...), 1, min, na.rm = T)
            min_col <- apply(sp_ID_dist(obj_sp1, obj_sp2_ng, ...), 2, min, na.rm = T)
          }
  )

  M <- (sum(min_row) + sum(min_col))/(length(min_col) + length(min_row))

  rm(list = ls()[ls() != "M"])

  return(M)

}

#' F function adaptation for polygons
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param r_min a \code{numeric} representing the mininmun distance to calculate
#' the function, \code{default = 0}.
#' @param r_max a \code{numeric} representing the maximum distance to calculate
#' the function \code{default = max(polygonds distance matrix)}.
#' @param by a \code{numeric} value that represents the how many units
#' between each \code{r} value.
#'
#' @return a \code{numeric} \code{matrix} corresponding to the estimated
#'  \eqn{F_{1,2}} and \eqn{F_{2,1}} function in the interval \eqn{[r_{min}, r_{max}]}.
#' @export
#'
pf12 <- function(obj_sp1, obj_sp2, r_min = 0, r_max = NULL, by = NULL) {
  m_dist <- sp_ID_dist(obj_sp1, obj_sp2)
  bbox <- obj_sp1@bbox

  if(is.null(r_max)) {
    r_x <- bbox[1,2] - bbox[1,1]
    r_y <- bbox[2,2] - bbox[2,1]
    r_max <- .4*max(r_x, r_y)
    rm(r_x, r_y)
  }

  if(is.null(by)) {
    by <- .5*(r_max - r_min)/sqrt(12)
  }

  r <- seq(from = r_min, to = r_max, by = by)

  output <- matrix(ncol = 3, nrow = length(r))
  colnames(output) <- c('r', 'pf12', 'pf21')
  output[ ,1] <- r

  ab <- apply(m_dist, 1, min, na.rm = T)
  ba <- apply(m_dist, 2, min, na.rm = T)

  for(i in seq_along(r)) {
    output[i, 2] <- mean(ab <= r[i])
    output[i, 3] <- mean(ba <= r[i])
  }

  rm(list = ls()[ls() != "output"])
  class(output) <- pF(output)

  return(output)
}

#' \eqn{K_{1,2}} - area-based adaptation for polygons
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param r_min a \code{numeric} representing the mininmun distance to calculate
#' the function, \code{default = 0}.
#' @param r_max a \code{numeric} representing the maximum distance to calculate
#' the function \code{default = max(polygonds distance matrix)}.
#' @param by a \code{numeric} value that represents the how many units
#' between each \code{r} value.
#' @param bbox a \code{matrix} giving the region of study.
#'
#' @return a \code{data.frame} corresponding to the estimated
#'  \eqn{K_{1,2}} function in the interval \eqn{[r_{min}, r_{max}]}.
#'
#' @export
#'
pk_area12m <- function(obj_sp1, obj_sp2, r_min = NULL, r_max = NULL, by = NULL, bbox) {
  # mat_dist <- sp_ID_dist(obj_sp1, obj_sp2)

  # if(sp::is.projected(obj_sp1)) obj_sp1 <- rgeos::gBuffer(obj_sp1, byid = T, width = 0)
  # if(sp::is.projected(obj_sp2)) obj_sp2 <- rgeos::gBuffer(obj_sp2, byid = T, width = 0)

  if(is.null(r_max)) {
    r_x <- bbox[1,2] - bbox[1,1]
    r_y <- bbox[2,2] - bbox[2,1]
    r_max <- .4*max(r_x, r_y)
    rm(r_x, r_y)
  }

  if(is.null(r_min)) {
    r_min <- 0.0001
  }

  if(is.null(by)) {
    by <- .5*(r_max - r_min)/sqrt(12)
  }

  # this calculations can be done outside the function
  N <- (bbox[1,2] - bbox[1,1])*(bbox[2,2] - bbox[2,1])
  tot_1 <- rgeos::gArea(obj_sp1)
  tot_2 <- rgeos::gArea(obj_sp2)
  l_1 <- tot_1/N
  l_2 <- tot_2/N

  r <- seq(from = r_min, to = r_max, by = by)

  output <- data.frame(r = rep(NA, length(r)), pk12 = rep(NA, length(r)))
  output$r <- r

  # areas_1 <- vector(mode = 'numeric', length = length(obj_sp1))
  # areas_2 <- vector(mode = 'numeric', length = length(obj_sp2))

  for(i in seq_along(r)) {

    buff_poly1 <- rgeos::gBuffer(obj_sp1, width = r[i])
    aux <- rgeos::gIntersection(buff_poly1, obj_sp2)
    areas_1 <- ifelse(is.null(aux), 0, rgeos::gArea(aux))
    rm(aux)

    buff_poly2 <- rgeos::gBuffer(obj_sp2, width = r[i])
    aux <- rgeos::gIntersection(buff_poly2, obj_sp1)
    areas_2 <- ifelse(is.null(aux), 0, rgeos::gArea(aux))
    rm(aux)

    k12 <- (l_2^-1)*areas_2/rgeos::gArea(rgeos::gIntersection(buff_poly1, limits_to_sp(obj_sp2@bbox)))
    k21 <- (l_1^-1)*areas_1/rgeos::gArea(rgeos::gIntersection(buff_poly2, limits_to_sp(obj_sp1@bbox)))
    output$pk12[i] <- (l_1*k21 + l_2*k12)/(l_1 + l_2)
  }

  rm(list = ls()[ls() != 'output'])
  return(output)

}

#' \eqn{K_{1,2}} - area-based adaptation for polygons - Modified
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param r_min a \code{numeric} representing the mininmun distance to calculate
#' the function, \code{default = 0}.
#' @param r_max a \code{numeric} representing the maximum distance to calculate
#' the function \code{default = max(polygonds distance matrix)}.
#' @param by a \code{numeric} value that represents the how many units
#' between each \code{r} value.
#' @param bbox a \code{matrix} giving the region of study.
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param ... parameters associated with edge correction.
#'
#' @return a \code{data.frame} corresponding to the estimated
#'  \eqn{K_{1,2}} function in the interval \eqn{[r_{min}, r_{max}]}.
#'
#' @export
#'
pk_area12 <- function(obj_sp1, obj_sp2, r_min = NULL, r_max = NULL, by = NULL, bbox,
                      correction = 'none', ...) {
  # mat_dist <- sp_ID_dist(obj_sp1, obj_sp2)

  # if(sp::is.projected(obj_sp1)) obj_sp1 <- rgeos::gBuffer(obj_sp1, byid = T, width = 0)
  # if(sp::is.projected(obj_sp2)) obj_sp2 <- rgeos::gBuffer(obj_sp2, byid = T, width = 0)

  if(is.null(r_max)) {
    r_x <- bbox[1,2] - bbox[1,1]
    r_y <- bbox[2,2] - bbox[2,1]
    r_max <- .15*max(r_x, r_y)
    rm(r_x, r_y)
  }

  if(is.null(r_min)) {
    r_min <- 0.0001
  }

  if(is.null(by)) {
    by <- (r_max - r_min)/sqrt(12)
  }

  # this calculations can be done outside the function
  N <- (bbox[1,2] - bbox[1,1])*(bbox[2,2] - bbox[2,1])
  # l_1 <- tot_1/N
  # l_2 <- tot_2/N

  r <- seq(from = r_min, to = r_max, by = by)/2

  output <- data.frame(r = rep(NA, length(r)), pk12 = rep(NA, length(r)))
  output$r <- r

  # areas_1 <- vector(mode = 'numeric', length = length(obj_sp1))
  # areas_2 <- vector(mode = 'numeric', length = length(obj_sp2))

  switch (correction,
          'none' = {
            tot_1 <- rgeos::gArea(obj_sp1)
            tot_2 <- rgeos::gArea(obj_sp2)
            for(i in seq_along(r)) {
              aux <- rgeos::gBuffer(obj_sp2, width = r[i])
              aux <- rgeos::gIntersection(aux, obj_sp1, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_1 <- rgeos::gArea(aux)
              }
              rm(aux)

              aux <- rgeos::gBuffer(obj_sp1, width = r[i])
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_2 <- rgeos::gArea(aux)
              }
              rm(aux)

              k12 <- tot_1*areas_2
              k21 <- tot_2*areas_1
              output$pk12[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'torus' = {
            tot_1 <- rgeos::gArea(obj_sp1)
            tot_2 <- rgeos::gArea(obj_sp2)
            for(i in seq_along(r)) {
              obj_sp1_t <- torus_corr(obj_sp1, bbox)
              obj_sp2_t <- torus_corr(obj_sp2, bbox)
              aux <- rgeos::gBuffer(obj_sp1_t, width = r[i], byid = T)
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_1 <- rgeos::gArea(aux)
              }
              rm(aux)

              aux <- rgeos::gBuffer(obj_sp2_t, width = r[i], T)
              aux <- rgeos::gIntersection(aux, obj_sp1, T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_2 <- rgeos::gArea(aux)
              }
              rm(aux)

              k12 <- tot_1*areas_2
              k21 <- tot_2*areas_1
              output$pk12[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'guard' = {
            new_bbox <- make_guard(obj_sp1@bbox, ...)
            obj_sp1_ng <- rgeos::gIntersection(obj_sp1, limits_to_sp(new_bbox), byid = T)
            obj_sp2_ng <- rgeos::gIntersection(obj_sp2, limits_to_sp(new_bbox), byid = T)
            tot_1 <- rgeos::gArea(obj_sp1_ng)
            tot_2 <- rgeos::gArea(obj_sp2_ng)
            for(i in seq_along(r)) {
              aux <- rgeos::gBuffer(obj_sp1_ng, width = r[i], byid = T)
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_1 <- rgeos::gArea(aux)
              }
              rm(aux)

              aux <- rgeos::gBuffer(obj_sp2_ng, width = r[i], byid = T)
              aux <- rgeos::gIntersection(aux, obj_sp1, byid = T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_2 <- rgeos::gArea(aux)
              }
              rm(aux)

              k12 <- tot_1*areas_2
              k21 <- tot_2*areas_1
              output$pk12[i] <- tryCatch(((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2)),
                                         error = function(cond) {
                                           return(NA)
                                         })
            }
          },
          'adjust' = {
            tot_1 <- rgeos::gArea(obj_sp1)
            tot_2 <- rgeos::gArea(obj_sp2)
            for(i in seq_along(r)) {
              aux <- rgeos::gBuffer(obj_sp1, width = r[i], byid = T)
              prop_inside_1 <- rgeos::gArea(rgeos::gIntersection(aux, limits_to_sp(obj_sp1@bbox)))/rgeos::gArea(aux)
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_1 <- rgeos::gArea(aux)
              }
              rm(aux)

              aux <- rgeos::gBuffer(obj_sp2, width = r[i], byid = T)
              prop_inside_2 <- rgeos::gArea(rgeos::gIntersection(aux, limits_to_sp(obj_sp2@bbox)))/rgeos::gArea(aux)
              aux <- rgeos::gIntersection(aux, obj_sp1, byid = T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_2 <- rgeos::gArea(aux)
              }
              rm(aux)
              k12 <- (tot_1/prop_inside_1)*areas_2
              k21 <- (tot_2/prop_inside_2)*areas_1
              output$pk12[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          })


  rm(list = ls()[ls() != 'output'])
  return(output)

}


#' \eqn{K_{1,2}} - distance-based adaptation for polygons
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param r_min a \code{numeric} representing the mininmun distance to calculate
#' the function, \code{default = 0}.
#' @param r_max a \code{numeric} representing the maximum distance to calculate
#' the function \code{default = max(polygonds distance matrix)}.
#' @param by a \code{numeric} value that represents the how many units
#' between each \code{r} value.
#' @param bbox a \code{matrix} giving the region of study
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param ... parameters associated with edge correction.
#'
#' @return a \code{data.frame} corresponding to the estimated
#'  \eqn{K_{1,2}} function in the interval \eqn{[r_{min}, r_{max}]}.
#'
#' @export
#'
pk_dist12 <- function(obj_sp1, obj_sp2, r_min = NULL, r_max = NULL, by = NULL, bbox,
                      correction = 'none', ...) {

  if(is.null(r_max)) {
    r_x <- bbox[1,2] - bbox[1,1]
    r_y <- bbox[2,2] - bbox[2,1]
    r_max <- .15*max(r_x, r_y)
    rm(r_x, r_y)
  }

  if(is.null(r_min)) {
    r_min <- 0.0001
  }

  if(is.null(by)) {
    by <- (r_max - r_min)/sqrt(12)
  }

  # this calculations can be done outside the function
  N <- (bbox[1,2] - bbox[1,1])*(bbox[2,2] - bbox[2,1])

  r <- seq(from = r_min, to = r_max, by = by)

  output <- data.frame(r = rep(NA, length(r)), pk12 = rep(NA, length(r)))
  output$r <- r

  switch (correction,
          'none' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            mat_dist <- sp_ID_dist(obj_sp1, obj_sp2, ...)
            for(i in seq_along(r)) {
              # k12 <- (l_2^(-1)) * sum(mat_dist < r[i])/N
              # k21 <- (l_1^(-1)) * sum(mat_dist < r[i])/N
              output$pk12[i] <- (N/(tot_1*tot_2))*sum(mat_dist < r[i], na.rm = T)
            }
          },
          'torus' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            # l_1 <- tot_1/N
            # l_2 <- tot_2/N
            obj_sp1_t <- torus_corr(obj_sp1, bbox)
            obj_sp2_t <- torus_corr(obj_sp2, bbox)
            mat_dist1 <- sp_ID_dist(obj_sp1, obj_sp2_t, ...)
            mat_dist2 <- sp_ID_dist(obj_sp1_t, obj_sp2, ...)
            for(i in seq_along(r)) {
              # k12 <- (tot_2^(-1)) * sum(mat_dist1 < r[i])
              # k21 <- (tot_1^(-1)) * sum(mat_dist2 < r[i])
              k12 <- sum(mat_dist1 < r[i], na.rm = T)*tot_2
              k21 <- sum(mat_dist2 < r[i], na.rm = T)*tot_1
              output$pk12[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'guard' = {
            new_bbox <- make_guard(obj_sp1@bbox, ...)
            obj_sp1_ng <- rgeos::gIntersection(obj_sp1, limits_to_sp(new_bbox), byid = T)
            obj_sp2_ng <- rgeos::gIntersection(obj_sp2, limits_to_sp(new_bbox), byid = T)
            tot_1 <- length(obj_sp1_ng)
            tot_1 <- ifelse(tot_1 == 0, 1e-16, tot_1)
            tot_2 <- length(obj_sp2_ng)
            tot_1 <- ifelse(tot_2 == 0, 1e-16, tot_2)
            mat_dist1 <- sp_ID_dist(obj_sp1_ng, obj_sp2, ...)
            mat_dist2 <- sp_ID_dist(obj_sp1, obj_sp2_ng, ...)
            for(i in seq_along(r)) {
              # k12 <- (l_2^(-1)) * sum(mat_dist1 < r[i])/N
              # k21 <- (l_1^(-1)) * sum(mat_dist2 < r[i])/N
              # output$pk12[i] <- ((tot_1*tot_2)^(-1))*N*(tot_2*k12 + tot_1*k21)
              k12 <- sum(mat_dist1 < r[i], na.rm = T)*tot_2
              k21 <- sum(mat_dist2 < r[i], na.rm = T)*tot_1
              output$pk12[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'adjust' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            for(i in seq_along(r)) {
              obj_sp1_bf <- rgeos::gBuffer(obj_sp1, width = r[i], byid = T)
              obj_sp2_bf <- rgeos::gBuffer(obj_sp2, width = r[i], byid = T)
              obj_sp1_bf <- rgeos::gArea(rgeos::gIntersection(obj_sp1_bf,
                                                              limits_to_sp(obj_sp1@bbox),byid = T),
                                         byid = T)/rgeos::gArea(obj_sp1_bf, byid = T)
              obj_sp2_bf <- rgeos::gArea(rgeos::gIntersection(obj_sp2_bf,
                                                              limits_to_sp(obj_sp2@bbox),
                                                              byid = T), byid = T)/rgeos::gArea(obj_sp2_bf, byid = T)
              mat_dist <- sp_ID_dist(obj_sp1, obj_sp2, ...)
              w <- (tot_1*obj_sp1_bf + tot_2*obj_sp2_bf)/(tot_1 + tot_2)
              output$pk12[i] <- (sum((mat_dist < r[i])/w[col(mat_dist)]))**(N/(tot_1*tot_2))
            }
          }
  )

  rm(list = ls()[ls() != "output"])

  return(output)
}
