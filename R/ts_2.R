#' \eqn{K_{1,2}} - area-based adaptation for polygons - 2
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
#' @param sim \code{boolean}
#' @param ... parameters associated with edge correction.
#'
#' @return a \code{data.frame} corresponding to the estimated
#'  \eqn{K_{1,2}} function in the interval \eqn{[r_{min}, r_{max}]}.
#'
#' @export
#'
karea <- function(obj_sp1, obj_sp2, r_min = NULL, r_max = NULL, by = NULL, bbox, sim = F, ...) {
  if(is.null(r_max) | is.null(r_min)) {
    r_x <- bbox[1,2] - bbox[1,1]
    r_y <- bbox[2,2] - bbox[2,1]
  }

  if(is.null(r_max)) {
    r_max <- .25*max(r_x, r_y)
  }

  if(is.null(r_min)) {
    r_min <- 0
  }

  if(is.null(by)) {
    by <- (r_max - r_min)/sqrt(12)
  }

  if(r_min == 0) {
    r_min <- by
  }

  # this calculations can be done outside the function
  N <- (bbox[1,2] - bbox[1,1])*(bbox[2,2] - bbox[2,1])

  r <- seq(from = r_min, to = r_max, by = by)/2

  # if(r[length(r)] != r_max/2) {
  #   r <- c(r, r_max/2)
  # }

  output <- data.frame(r = r, pk12 = rep(NA, length(r)))

  # areas_1 <- vector(mode = 'numeric', length = length(obj_sp1))
  # areas_2 <- vector(mode = 'numeric', length = length(obj_sp2))

  tot_1 <- rgeos::gArea(obj_sp1)
  if(sim) {
    tot_2 <- rgeos::gArea(obj_sp2)/4
  } else {
    tot_2 <- rgeos::gArea(obj_sp2)
  }
  for(i in seq_along(r)) {
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
    output$pk12[i] <- k12*(N/(tot_1*tot_2))
  }

  rm(list = ls()[ls() != 'output'])
  return(output)

}

#' \eqn{K_{1,2}} - distance-based adaptation for polygons 2
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
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param sim \code{boolean}
#' @param ... parameters associated with edge correction.
#'
#' @return a \code{data.frame} corresponding to the estimated
#'  \eqn{K_{1,2}} function in the interval \eqn{[r_{min}, r_{max}]}.
#'
#' @export
#'
kdist <- function(obj_sp1, obj_sp2, r_min = NULL, r_max = NULL, by = NULL, bbox,
                  sim = F, ...) {
  if(is.null(r_max) | is.null(r_min)) {
    r_x <- bbox[1,2] - bbox[1,1]
    r_y <- bbox[2,2] - bbox[2,1]
  }

  if(is.null(r_max)) {
    r_max <- .25*max(r_x, r_y)
  }

  if(is.null(r_min)) {
    r_min <- 0
  }

  if(is.null(by)) {
    by <- (r_max - r_min)/sqrt(12)
  }

  if(r_min == 0) {
    r_min <- by
  }

  N <- (bbox[1,2] - bbox[1,1])*(bbox[2,2] - bbox[2,1])

  r <- seq(from = r_min, to = r_max, by = by)

  # if(r[length(r)] != r_max) {
  #   r <- c(r, r_max)
  # }

  output <- data.frame(r = r, pk12 = rep(NA, length(r)))

  tot_1 <- length(obj_sp1)
  if(sim) {
    tot_2 <- length(obj_sp2)/4
  } else {
    tot_2 <- length(obj_sp2)
  }
  mat_dist <- sp_ID_dist(obj_sp1, obj_sp2, ...)
  for(i in seq_along(r)) {
    output$pk12[i] <- (N/(tot_1*tot_2))*sum(mat_dist <= r[i], na.rm = T)
  }

  rm(list = ls()[ls() != "output"])

  return(output)
}

