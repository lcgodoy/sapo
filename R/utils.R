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
      npolyaux <- length(obj_sp@polygons[[i]]@Polygons)
      if(!npolyaux > 1) {
        aux <- obj_sp@polygons[[i]]@Polygons[[1]]@coords
        aux[,1] <- aux[,1] + jump_x
        aux[,2] <- aux[,2] + jump_y
        attr(obj_sp@polygons[[i]]@Polygons[[1]], "coords") <- aux
      } else {
        for(k in seq_along(npolyaux)) {
          aux <- obj_sp@polygons[[i]]@Polygons[[k]]@coords
          aux[,1] <- aux[,1] + jump_x
          aux[,2] <- aux[,2] + jump_y
          attr(obj_sp@polygons[[i]]@Polygons[[k]], "coords") <- aux
        }
      }

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
#' @param m a \code{character} equal to \code{c('haus1', 'haus2')}.
#' If not supplied, euclidean distance is performed.
#'
#' @importFrom rgeos gDistance
#'
#' @return a distance matrix
#'
sp_ID_dist <- function(obj_sp1, obj_sp2, m = NULL) {
  if(is.null(m)) {
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
  if(m == 'haus1') {
    return(sp_ID_haus(obj_sp1 = obj_sp1,
                      obj_sp2 = obj_sp2))
  }
}

#' Hausdorff Distance Matrix - ID (Method 1)
#'
#' @description Calculate the hausdorff distance between polygons (or points)
#'              from different types considering their ID's.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#'
#' @importFrom rgeos gDistance
#'
#' @return a distance matrix
#'
sp_ID_haus <- function(obj_sp1, obj_sp2) {
  nameaux1 <- names(obj_sp1)

  n_poly1 <- length(unique(nameaux1))

  nameaux2 <- names(obj_sp2)

  n_poly2 <- length(unique(nameaux2))

  matrix_dist <- matrix(data = rep(0, n_poly1*n_poly2), nrow = n_poly2,
                        ncol = n_poly1, dimnames = list(unique(nameaux2),
                                                        unique(nameaux1)))

  matrix_aux <- gDistance(obj_sp1, obj_sp2, byid = T, hausdorff = T)

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


#' Make Guard
#'
#' Internal use
#'
#' @param bbox \code{numeric matrix}
#' @param guard_perc \code{numeric} between 0 and 1
#'
#' @return \code{numeric matrix}
#'
make_guard <- function(bbox, guard_perc = 0.125) {
  range_x <- bbox[1,2] - bbox[1,1]
  range_y <- bbox[2,2] - bbox[2,1]

  out <- bbox
  out[1, 1] <- out[1, 1] + guard_perc*range_x
  out[1, 2] <- out[1, 2] - guard_perc*range_x
  out[2, 1] <- out[2, 1] + guard_perc*range_y
  out[2, 2] <- out[2, 2] - guard_perc*range_y
  return(out)
}

#' H(d) - Euclidean Distance
#'
#' @description Auxiliar function to \code{\link{h_func}}.
#'
#' @details Internal use.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param distances a \code{numeric vector} indicating the distances to evaluate \eqn{H(d)}. If
#' \code{NULL} then the range considered goes from 5% to 20% of the max distance that can be
#' observed inside the \code{unique_bbox}.
#' @param H a \code{character} indicating the function to be used. Possible entries are:
#' \code{'K'} or \code{'L'}.
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains both sets
#' @param ... parameters for test statistics functions
#'
#' @return \code{numeric vector}.
h_euc <- function(obj_sp1, obj_sp2, correction, distances, H, unique_bbox, ...) {
  N <- (unique_bbox[1,2] - unique_bbox[1,1])*(unique_bbox[2,2] - unique_bbox[2,1])
  output <- vector(mode = 'numeric', length = length(distances))
  switch (correction,
          'none' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            mat_dist <- sp_ID_dist(obj_sp1, obj_sp2)
            for(i in seq_along(distances)) {
              output[i] <- (N/(tot_1*tot_2))*sum(mat_dist < distances[i], na.rm = T)
            }
          },
          'torus' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            obj_sp1_t <- torus_corr(obj_sp1, unique_bbox)
            obj_sp2_t <- torus_corr(obj_sp2, unique_bbox)
            mat_dist1 <- sp_ID_dist(obj_sp1, obj_sp2_t)
            mat_dist2 <- sp_ID_dist(obj_sp1_t, obj_sp2)
            for(i in seq_along(distances)) {
              k12 <- sum(mat_dist1 < distances[i], na.rm = T)*tot_2
              k21 <- sum(mat_dist2 < distances[i], na.rm = T)*tot_1
              output[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'guard' = {
            new_bbox <- make_guard(obj_sp1@bbox, ...)
            obj_sp1_ng <- rgeos::gIntersection(obj_sp1, limits_to_sp(new_bbox), byid = T)
            obj_sp2_ng <- rgeos::gIntersection(obj_sp2, limits_to_sp(new_bbox), byid = T)
            tot_1 <- length(obj_sp1_ng)
            tot_1 <- ifelse(tot_1 == 0, 1e-10, tot_1)
            tot_2 <- length(obj_sp2_ng)
            tot_1 <- ifelse(tot_2 == 0, 1e-10, tot_2)
            mat_dist1 <- sp_ID_dist(obj_sp1_ng, obj_sp2)
            mat_dist2 <- sp_ID_dist(obj_sp1, obj_sp2_ng)
            for(i in seq_along(distances)) {
              k12 <- sum(mat_dist1 < distances[i], na.rm = T)*tot_2
              k21 <- sum(mat_dist2 < distances[i], na.rm = T)*tot_1
              output[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'adjust' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            nm_1 <- row.names(obj_sp1)
            nm_2 <- row.names(obj_sp2)
            nm_aux_1 <- as.character(1:length(obj_sp1))
            nm_aux_2 <- as.character(1:length(obj_sp2))
            row.names(obj_sp1) <- nm_aux_1
            row.names(obj_sp2) <- nm_aux_2
            for(i in seq_along(distances)) {
              obj_sp1_bf <- rgeos::gBuffer(obj_sp1,
                                           width = distances[i], byid = T)
              obj_sp2_bf <- rgeos::gBuffer(obj_sp2,
                                           width = distances[i], byid = T)
              row.names(obj_sp1_bf) <- nm_1
              row.names(obj_sp2_bf) <- nm_2
              obj_sp1_bf <- rgeos::gArea(rgeos::gIntersection(obj_sp1_bf,
                                                              limits_to_sp(obj_sp1@bbox),byid = T),
                                         byid = T)/rgeos::gArea(obj_sp1_bf, byid = T)
              obj_sp2_bf <- rgeos::gArea(rgeos::gIntersection(obj_sp2_bf,
                                                              limits_to_sp(obj_sp2@bbox),
                                                              byid = T), byid = T)/rgeos::gArea(obj_sp2_bf, byid = T)

              row.names(obj_sp1) <- nm_1
              row.names(obj_sp2) <- nm_2
              mat_dist <- sp_ID_dist(obj_sp1, obj_sp2)
              row.names(obj_sp1) <- nm_aux_1
              row.names(obj_sp2) <- nm_aux_2
              w <- (tot_1*obj_sp1_bf + tot_2*obj_sp2_bf)/(tot_1 + tot_2)
              output[i] <- (sum((mat_dist < distances[i])/w[col(mat_dist)]))*(N/(tot_1*tot_2))
            }
          }
  )

  if(H == 'L') {
    output <- sqrt(output/pi)
  }

  return(output)
}

#' H(d) - Hausdorff Distance
#'
#' @description Auxiliar function to \code{\link{h_func}}.
#'
#' @details Internal use.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param distances a \code{numeric vector} indicating the distances to evaluate \eqn{H(d)}. If
#' \code{NULL} then the range considered goes from 5% to 20% of the max distance that can be
#' observed inside the \code{unique_bbox}.
#' @param H a \code{character} indicating the function to be used. Possible entries are:
#' \code{'K'} or \code{'L'}.
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains both sets
#' @param ... parameters for test statistics functions
#'
#' @return \code{numeric vector}.
h_haus <- function(obj_sp1, obj_sp2, correction, distances, H, unique_bbox, ...) {
  N <- (unique_bbox[1,2] - unique_bbox[1,1])*(unique_bbox[2,2] - unique_bbox[2,1])
  output <- vector(mode = 'numeric', length = length(distances))
  switch (correction,
          'none' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            mat_dist <- sp_ID_haus(obj_sp1, obj_sp2)
            for(i in seq_along(distances)) {
              output[i] <- (N/(tot_1*tot_2))*sum(mat_dist < distances[i], na.rm = T)
            }
          },
          'torus' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            obj_sp1_t <- torus_corr(obj_sp1, unique_bbox)
            obj_sp2_t <- torus_corr(obj_sp2, unique_bbox)
            mat_dist1 <- sp_ID_haus(obj_sp1, obj_sp2_t)
            mat_dist2 <- sp_ID_haus(obj_sp1_t, obj_sp2)
            for(i in seq_along(distances)) {
              k12 <- sum(mat_dist1 < distances[i], na.rm = T)*tot_2
              k21 <- sum(mat_dist2 < distances[i], na.rm = T)*tot_1
              output[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'guard' = {
            new_bbox <- make_guard(obj_sp1@bbox, ...)
            obj_sp1_ng <- rgeos::gIntersection(obj_sp1, limits_to_sp(new_bbox), byid = T)
            obj_sp2_ng <- rgeos::gIntersection(obj_sp2, limits_to_sp(new_bbox), byid = T)
            tot_1 <- length(obj_sp1_ng)
            tot_1 <- ifelse(tot_1 == 0, 1e-10, tot_1)
            tot_2 <- length(obj_sp2_ng)
            tot_1 <- ifelse(tot_2 == 0, 1e-10, tot_2)
            mat_dist1 <- sp_ID_haus(obj_sp1_ng, obj_sp2)
            mat_dist2 <- sp_ID_haus(obj_sp1, obj_sp2_ng)
            for(i in seq_along(distances)) {
              k12 <- sum(mat_dist1 < distances[i], na.rm = T)*tot_2
              k21 <- sum(mat_dist2 < distances[i], na.rm = T)*tot_1
              output[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'adjust' = {
            tot_1 <- length(obj_sp1)
            tot_2 <- length(obj_sp2)
            nm_1 <- row.names(obj_sp1)
            nm_2 <- row.names(obj_sp2)
            nm_aux_1 <- as.character(1:length(obj_sp1))
            nm_aux_2 <- as.character(1:length(obj_sp2))
            row.names(obj_sp1) <- nm_aux_1
            row.names(obj_sp2) <- nm_aux_2
            for(i in seq_along(distances)) {
              obj_sp1_bf <- rgeos::gBuffer(obj_sp1,
                                           width = distances[i], byid = T)
              obj_sp2_bf <- rgeos::gBuffer(obj_sp2,
                                           width = distances[i], byid = T)
              row.names(obj_sp1_bf) <- nm_1
              row.names(obj_sp2_bf) <- nm_2
              obj_sp1_bf <- rgeos::gArea(rgeos::gIntersection(obj_sp1_bf,
                                                              limits_to_sp(obj_sp1@bbox),byid = T),
                                         byid = T)/rgeos::gArea(obj_sp1_bf, byid = T)
              obj_sp2_bf <- rgeos::gArea(rgeos::gIntersection(obj_sp2_bf,
                                                              limits_to_sp(obj_sp2@bbox),
                                                              byid = T), byid = T)/rgeos::gArea(obj_sp2_bf, byid = T)

              row.names(obj_sp1) <- nm_1
              row.names(obj_sp2) <- nm_2
              mat_dist <- sp_ID_haus(obj_sp1, obj_sp2)
              row.names(obj_sp1) <- nm_aux_1
              row.names(obj_sp2) <- nm_aux_2
              w <- (tot_1*obj_sp1_bf + tot_2*obj_sp2_bf)/(tot_1 + tot_2)
              output[i] <- (sum((mat_dist < distances[i])/w[col(mat_dist)]))*(N/(tot_1*tot_2))
            }
          }
  )

  if(H == 'L') {
    output <- sqrt(output/pi)
  }

  return(output)
}

#' H(d) - Area-based
#'
#' @description Auxiliar function to \code{\link{h_func}}.
#'
#' @details Internal use.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param distances a \code{numeric vector} indicating the distances to evaluate \eqn{H(d)}. If
#' \code{NULL} then the range considered goes from 5% to 20% of the max distance that can be
#' observed inside the \code{unique_bbox}.
#' @param H a \code{character} indicating the function to be used. Possible entries are:
#' \code{'K'} or \code{'L'}.
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains both sets
#' @param ... parameters for test statistics functions
#'
#' @return \code{numeric vector}.
h_area <- function(obj_sp1, obj_sp2, correction, distances, H, unique_bbox, ...) {
  distances <- distances/2
  output <- vector(mode = 'numeric', length = length(distances))
  N <- (unique_bbox[1,2] - unique_bbox[1,1])*(unique_bbox[2,2] - unique_bbox[2,1])

  switch (correction,
          'none' = {
            tot_1 <- rgeos::gArea(obj_sp1)
            tot_2 <- rgeos::gArea(obj_sp2)
            for(i in seq_along(distances)) {
              aux <- rgeos::gBuffer(obj_sp2, width = distances[i])
              aux <- rgeos::gIntersection(aux, obj_sp1, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                areas_1 <- rgeos::gArea(aux)
              }
              rm(aux)

              aux <- rgeos::gBuffer(obj_sp1, width = distances[i])
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                areas_2 <- rgeos::gArea(aux)
              }

              k12 <- tot_1*areas_2
              k21 <- tot_2*areas_1
              output[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'torus' = {
            tot_1 <- rgeos::gArea(obj_sp1)
            tot_2 <- rgeos::gArea(obj_sp2)
            for(i in seq_along(distances)) {
              obj_sp1_t <- torus_corr(obj_sp1, unique_bbox)
              obj_sp2_t <- torus_corr(obj_sp2, unique_bbox)
              aux <- rgeos::gBuffer(obj_sp1_t, width = distances[i], byid = F)
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                # row.names(aux) <- stringr::str_extract(string = row.names(aux), pattern = '^[^\\s]+')
                areas_1 <- rgeos::gArea(aux)
              }

              aux <- rgeos::gBuffer(obj_sp2_t, width = distances[i], F)
              aux <- rgeos::gIntersection(aux, obj_sp1, T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                areas_2 <- rgeos::gArea(aux)
              }

              k12 <- tot_1*areas_2
              k21 <- tot_2*areas_1
              output[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          },
          'guard' = {
            new_bbox <- make_guard(obj_sp1@bbox, ...)
            obj_sp1_ng <- rgeos::gIntersection(obj_sp1, limits_to_sp(new_bbox), byid = T)
            obj_sp2_ng <- rgeos::gIntersection(obj_sp2, limits_to_sp(new_bbox), byid = T)
            tot_1 <- rgeos::gArea(obj_sp1_ng)
            tot_2 <- rgeos::gArea(obj_sp2_ng)
            for(i in seq_along(distances)) {
              aux <- rgeos::gBuffer(obj_sp1_ng, width = distances[i], byid = T)
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                areas_1 <- rgeos::gArea(aux)
              }
              rm(aux)

              aux <- rgeos::gBuffer(obj_sp2_ng, width = distances[i], byid = T)
              aux <- rgeos::gIntersection(aux, obj_sp1, byid = T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                areas_2 <- rgeos::gArea(aux)
              }

              k12 <- tot_1*areas_2
              k21 <- tot_2*areas_1
              output[i] <- tryCatch(((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2)),
                                    error = function(cond) {
                                      return(NA)
                                    })
            }
          },
          'adjust' = {
            tot_1 <- rgeos::gArea(obj_sp1)
            tot_2 <- rgeos::gArea(obj_sp2)
            for(i in seq_along(distances)) {
              aux <- rgeos::gBuffer(obj_sp1, width = distances[i], byid = T)
              prop_inside_1 <- rgeos::gArea(rgeos::gIntersection(aux, limits_to_sp(obj_sp1@bbox)))/rgeos::gArea(aux)
              aux <- rgeos::gIntersection(aux, obj_sp2, byid = T)
              if(is.null(aux)) {
                areas_1 <- 0
              } else {
                areas_1 <- rgeos::gArea(aux)
              }
              rm(aux)

              aux <- rgeos::gBuffer(obj_sp2, width = distances[i], byid = T)
              prop_inside_2 <- rgeos::gArea(rgeos::gIntersection(aux, limits_to_sp(obj_sp2@bbox)))/rgeos::gArea(aux)
              aux <- rgeos::gIntersection(aux, obj_sp1, byid = T)
              if(is.null(aux)) {
                areas_2 <- 0
              } else {
                areas_2 <- rgeos::gArea(aux)
              }
              rm(aux)
              k12 <- (tot_1/prop_inside_1)*areas_2
              k21 <- (tot_2/prop_inside_2)*areas_1
              output[i] <- ((k12 + k21)/(tot_1 + tot_2))*(N/(tot_1*tot_2))
            }
          })

  if(H == 'L') {
    output <- sqrt(output/pi)
  }

  return(output)
}

#' \eqn{H_{12}(d)}
#'
#' @description Evaluate a given function in the inputted distances.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#' @param distances a \code{numeric vector} indicating the distances to evaluate \eqn{H(d)}. If
#' \code{NULL} then the range considered goes from 5% to 20% of the max distance that can be
#' observed inside the \code{unique_bbox}.
#' @param H a \code{character} indicating the function to be used. Possible entries are:
#' \code{'K'} or \code{'L'}.
#' @param unique_bbox a \code{matrix} \eqn{2 \times 2} corresponding to the boundary box
#'  that contains both sets
#' @param method a \code{character} specifying which kind of distance will be used
#' to evalueate \eqn{H}. Also, there is an option of using areas. Options available:
#' \code{c('hausdorff', 'euclidean', 'area')}.
#' @param ... parameters for test statistics functions
#'
#' @export
#'
#' @return \code{numeric vector}.
h_func <- function(obj_sp1, obj_sp2, unique_bbox = NULL, distances = NULL,
                   method = 'hausdorff', H = 'L', correction = 'none', ...) {
  if((! "SpatialPolygons" %in% class(obj_sp1)) | (! "SpatialPolygons" %in% class(obj_sp2)))
    stop("obj_sp1 and obj_sp2 must be from class 'SpatialPolygons'")

  if(length(H) > 1 | length(correction) > 1 | length(method) > 1)
    stop('correction, method, and H should have length == 1L.')

  if(! correction %in% c('none', 'adjust', 'guard', 'torus'))
    stop("correction must be 'none', 'adjust', 'guard' or 'torus'.")

  if(! method %in% c('hausdorff', 'euclidean', 'area'))
    stop("method must be 'hausdorff', 'euclidean' or 'area'.")

  if(! H %in% c('K', 'L'))
    stop("correction must be 'K' or 'L'.")

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
    max_d <- sqrt((unique_bbox[1, 2] - unique_bbox[1, 1])^2 + (unique_bbox[2, 2] - unique_bbox[2, 1])^2)
    distances <- seq(from = .05*max_d, to = .2*max_d, length.out = 16L)
  }

  return(switch(method,
                'hausdorff' = {
                  h_haus(obj_sp1, obj_sp2, correction, distances, H, unique_bbox, ...)
                },
                'euclidean' = {
                  h_euc(obj_sp1, obj_sp2, correction, distances, H, unique_bbox, ...)
                },
                'area' = {
                  h_area(obj_sp1, obj_sp2, correction, distances, H, unique_bbox, ...)
                }))
}
