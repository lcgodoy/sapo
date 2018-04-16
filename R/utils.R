#' Limits to Spatial
#'
#' @description Internal function.
#'
#' @param b_box Boundary box from class \code{matrix}
#'
#' @return an object from class \code{SpatialPolygons}
#' @importFrom sp SpatialPolygons Polygons
#'
limits_to_sp <- function(b_box) {
  coords <- matrix(c(b_box[1,1], b_box[2,1],
                     b_box[1,1], b_box[2,2],
                     b_box[1,2], b_box[2,2],
                     b_box[1,2], b_box[2,1]
  ),
  ncol = 2, byrow = TRUE)

  out <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), ID = " ")))
  return(out)
}

#' Auxiliar function to \code{poly_shift}
#'
#' @param obj_sp a spatial object
#' @param obj_sp2 a spatial object
#' @param obj_sp3 a spatial object
#' @param obj_sp4 a spatial object
#'
#' @return a spatial object
shift_aux_point <- function(obj_sp, obj_sp2, obj_sp3, obj_sp4) {
  output <- rbind(obj_sp, obj_sp2,
                  obj_sp3, obj_sp4)

  return(output)
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
poly_shift <- function(obj_sp, bbox_tot = NULL) {
  if(!is.null(bbox_tot)) attr(obj_sp, "bbox") <- bbox_tot

  range_x <- (max(obj_sp@bbox[1,]) - min(obj_sp@bbox[1,]))
  range_y <- (max(obj_sp@bbox[2,]) - min(obj_sp@bbox[2,]))

  n_poly <- length(obj_sp)

  if(class(obj_sp) %in% "SpatialPolygons") {
    # Obj 1 = objsp

    # Obj 2 = copy 1 - translation in x

    obj_sp2 <- obj_sp

    aux <- obj_sp2@bbox

    # X Shift

    aux[1,] <- aux[1,] + range_x

    # Defining shifted bbox

    attr(obj_sp2, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    #index_aux <- as.numeric(row.names(obj2sp))

    for(i in 1:length(obj_sp2)) {
      aux <- obj_sp2@polygons[[i]]@Polygons[[1]]@coords
      aux[,1] <- aux[,1] + range_x
      attr(obj_sp2@polygons[[i]]@Polygons[[1]], "coords") <- aux
    }

    rm(aux)

    # Obj 3 = copy 2 - translation in y

    obj_sp3 <- obj_sp

    aux <- obj_sp3@bbox

    # Y Shift

    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp3, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    for(i in 1:length(obj_sp3)) {
      aux <- obj_sp3@polygons[[i]]@Polygons[[1]]@coords
      aux[,2] <- aux[,2] + range_y
      attr(obj_sp3@polygons[[i]]@Polygons[[1]], "coords") <- aux
    }

    rm(aux)

    # Obj 4 = copy 3 - translation in y

    obj_sp4 <- obj_sp

    aux <- obj_sp4@bbox

    # X and Y Shift

    aux[1,] <- aux[1,] + range_x
    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp4, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    for(i in 1:length(obj_sp4)) {
      aux <- obj_sp4@polygons[[i]]@Polygons[[1]]@coords
      aux[,1] <- aux[,1] + range_x
      aux[,2] <- aux[,2] + range_y
      attr(obj_sp4@polygons[[i]]@Polygons[[1]], "coords") <- aux
    }

    output <- shift_aux(obj_sp, obj_sp2, obj_sp3, obj_sp4,
                        length(obj_sp), range_x, range_y)
  } else {

    # Obj 1 = objsp

    # Obj 2 = copy 1 - translation in x

    obj_sp2 <- obj_sp

    aux <- obj_sp2@bbox

    # X Shift

    aux[1,] <- aux[1,] + range_x

    # Defining shifted bbox

    attr(obj_sp2, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    #index_aux <- as.numeric(row.names(obj2sp))

    aux <- obj_sp2@coords
    aux[, 1] <- aux[, 1] + range_x
    attr(obj_sp2, "coords") <- aux


    rm(aux)

    # Obj 3 = copy 2 - translation in y

    obj_sp3 <- obj_sp

    aux <- obj_sp3@bbox

    # Y Shift

    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp3, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    aux <- obj_sp3@coords
    aux[, 2] <- aux[, 2] + range_y
    attr(obj_sp3, 'coords') <- aux

    rm(aux)

    # Obj 4 = copy 3 - translation in y

    obj_sp4 <- obj_sp

    aux <- obj_sp4@bbox

    # X and Y Shift

    aux[1,] <- aux[1,] + range_x
    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp4, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    aux <- obj_sp4@coords
    aux[,1] <- aux[,1] + range_x
    aux[,2] <- aux[,2] + range_y
    attr(obj_sp4, "coords") <- aux

    output <- shift_aux_point(obj_sp, obj_sp2, obj_sp3, obj_sp4)
  }

  rm(list = ls()[ls() != "output"])

  return(output)
}

#' Polygons' Random Shift
#'
#' @param obj_sp object from class \code{SpatialPolygons}
#' @param bbox_max Boundary box from class \code{matrix}
#'
#' @importFrom stats runif
#'
#' @return an object from class \code{SpatialPolygons} randomly translated
#'
poly_rf <- function(obj_sp, bbox_max) {
  range_x <- (max(obj_sp@bbox[1,]) - min(obj_sp@bbox[1,]))
  range_y <- (max(obj_sp@bbox[2,]) - min(obj_sp@bbox[2,]))

  bbox_max2 <- bbox_max

  bbox_max2[1,1] <- bbox_max[1,1] + range_x/2
  bbox_max2[1,2] <- bbox_max[1,2] - range_x/2
  bbox_max2[2,1] <- bbox_max[2,1] + range_y/2
  bbox_max2[2,2] <- bbox_max[2,2] - range_y/2

  n_x <- runif(1, bbox_max2[1,1], bbox_max2[1,2])
  n_y <- runif(1, bbox_max2[2,1], bbox_max2[2,2])

  bbox2 <- matrix(c(xmin = n_x - range_x/2,
                    xmax = n_x + range_x/2,
                    ymin = n_y - range_y/2,
                    ymax = n_y + range_y/2),
                  nrow = 2, byrow = T, dimnames = list(c("x", "y"),
                                                       c("min", "max")))

  jump_x <- bbox2[1,1] - sp::bbox(obj_sp)[1,1]
  jump_y <- bbox2[2,1] - sp::bbox(obj_sp)[2,1]

  attr(obj_sp, "bbox") <- bbox2

  if(class(obj_sp) %in% 'SpatialPolygons') {
    for(i in 1:length(obj_sp)) {
      aux <- obj_sp@polygons[[i]]@Polygons[[1]]@coords
      aux[,1] <- aux[,1] + jump_x
      aux[,2] <- aux[,2] + jump_y
      attr(obj_sp@polygons[[i]]@Polygons[[1]], "coords") <- aux
    }
  } else {
    aux <- obj_sp@coords
    aux[,1] <- aux[,1] + jump_x
    aux[,2] <- aux[,2] + jump_y
    attr(obj_sp, "coords") <- aux
  }

  return(obj_sp)

}

#' Distance Matrix - ID
#'
#' @description Calculate the distance between polygons (or points)
#'              from different types considering their ID's.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#'
#' @importFrom rgeos gDistance
#'
#' @return a distance matrix
#'
sp_ID_dist <- function(obj_sp1, obj_sp2) {
  nameaux1 <- names(obj_sp1)

  n_poly1 <- length(unique(nameaux1))

  nameaux2 <- names(obj_sp2)

  n_poly2 <- length(unique(nameaux2))

  matrix_dist <- matrix(data = rep(0, n_poly1*n_poly2), nrow = n_poly2,
                        ncol = n_poly1, dimnames = list(unique(nameaux2),
                                                        unique(nameaux1)))

  matrix_aux <- gDistance(obj_sp1, obj_sp2, byid = T)

  if(any(duplicated(colnames(matrix_aux)))) {
    for(i in 1:ncol(matrix_dist)) {
      aux <- matrix_aux[, which(colnames(matrix_aux) == colnames(matrix_dist)[i]), drop = F]
      matrix_dist[, i] <- apply(aux, 1, min)
      rm(aux)
    }
  }
  if(any(duplicated(rownames(matrix_aux)))) {
    for(i in 1:nrow(matrix_dist)) {
      aux <- matrix_aux[which(rownames(matrix_aux) == rownames(matrix_dist)[i]), , drop = F]
      matrix_dist[i, ] <- apply(aux, 2, min)
      rm(aux)
    }
  }
  if(!(any(duplicated(colnames(matrix_aux)))) & (!any(duplicated(rownames(matrix_aux))))) {
    matrix_dist <- matrix_aux
  }

  rm(list = ls()[ls() != "matrix_dist"])

  return(matrix_dist)

}

#' Polygons' Spatial Association Measure
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#'
#' @return a \code{numeric} \code{scalar} with the polygons' spatial
#'  association measure for the corresponding objects.
#' @export
#'
psam <- function(obj_sp1, obj_sp2) {
  m_dist <- sp_ID_dist(obj_sp1, obj_sp2)

  min_row <- apply(m_dist, 1, min)
  min_col <- apply(m_dist, 2, min)

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

  ab <- apply(m_dist, 1, min)
  ba <- apply(m_dist, 2, min)

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
pk_area12 <- function(obj_sp1, obj_sp2, r_min = NULL, r_max = NULL, by = NULL, bbox) {
  # mat_dist <- sp_ID_dist(obj_sp1, obj_sp2)

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

  areas_1 <- vector(mode = 'numeric', length = length(obj_sp1))
  areas_2 <- vector(mode = 'numeric', length = length(obj_sp2))

  for(i in seq_along(r)) {

    for(j in seq_along(obj_sp1)) {
      aux <- rgeos::gBuffer(obj_sp1[j], width = r[i])
      aux <- rgeos::gIntersection(aux, obj_sp2)
      areas_1[j] <- ifelse(is.null(aux), 0, rgeos::gArea(aux))
      rm(aux)
    }

    for(j in seq_along(obj_sp2)) {
      aux <- rgeos::gBuffer(obj_sp2[j], width = r[i])
      aux <- rgeos::gIntersection(aux, obj_sp1)
      areas_2[j] <- ifelse(is.null(aux), 0, rgeos::gArea(aux))
      rm(aux)
    }

    k12 <- (l_2^(-1)) * sum(areas_2)/N
    k21 <- (l_1^(-1)) * sum(areas_1)/N
    output$pk12[i] <- (tot_1*k21 + tot_2*k12)/(tot_1 + tot_2)
  }

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
#'
#' @return a \code{data.frame} corresponding to the estimated
#'  \eqn{K_{1,2}} function in the interval \eqn{[r_{min}, r_{max}]}.
#'
#' @export
#'
pk_dist12 <- function(obj_sp1, obj_sp2, r_min = NULL, r_max = NULL, by = NULL, bbox) {

  if(sp::is.projected(obj_sp1)) obj_sp1 <- rgeos::gBuffer(obj_sp1, byid = T, width = 0)
  if(sp::is.projected(obj_sp2)) obj_sp2 <- rgeos::gBuffer(obj_sp2, byid = T, width = 0)

  mat_dist <- sp_ID_dist(obj_sp1, obj_sp2)

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

  r <- seq(from = r_min, to = r_max, by = by)

  output <- data.frame(r = rep(NA, length(r)), pk12 = rep(NA, length(r)))
  output$r <- r

  tot_1 <- length(obj_sp1)
  tot_2 <- length(obj_sp2)
  l_1 <- tot_1/N
  l_2 <- tot_2/N

  for(i in seq_along(r)) {
    k12 <- (l_2^(-1)) * sum(mat_dist < r[i])/N
    k21 <- (l_1^(-1)) * sum(mat_dist < r[i])/N
    output$pk12[i] <- (tot_1*k21 + tot_2*k12)/(tot_1 + tot_2)
  }

  rm(list = ls()[ls() != "output"])

  return(output)
}
