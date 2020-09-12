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

#' b(P1, P2, r) function
#'
#' Internal use
#'
#' @param r \code{numeric scalar} a distance
#' @param poly_1 \code{SpatialPolygons} - Polygon pattern 1
#' @param poly_2 \code{SpatialPolygons} - Polygon pattern 2
#'
#' @return \code{numeric scaler}
b_12 <- function(r, poly_1, poly_2) {
    x <- sf::st_as_sfc(poly_1)
    y <- sf::st_as_sfc(poly_2)

    xx <- vector(mode = "list",
                 length = length(r))

    xy_bool <- vector(mode = "list",
                      length = length(r))

    intsection <- vector(mode = "list",
                         length = length(r))

    b12 <- vector(mode = "numeric", length = length(r))

    for(i in seq_along(r)) {
        xx[[i]] <- sf::st_buffer(x = x, dist = r[i],
                                 nQuadSegs = 5L)
        xy_bool[[i]] <- sf::st_intersects(x = xx[[i]],
                                          y = y,
                                          sparse = FALSE)
        if(any(xy_bool[[i]])) {
            b12[i] <- sum(sf::st_area(sf::st_intersection(xx[[i]], y)))
        } else {
            b12[i] <- 0
        }
    }
    return(b12)
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
h_euc <- function(obj_sp1, obj_sp2, distances, H, unique_bbox, ...) {
    N <- diff(unique_bbox[1, ]) * diff(unique_bbox[2, ])

    output <- vector(mode = 'numeric', length = length(distances))

    mat_dist <- rgeos::gDistance(obj_sp1, obj_sp2, byid = TRUE)
    output <- calc_h(x = mat_dist, dists = distances)
    output <- output*N
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
h_haus <- function(obj_sp1, obj_sp2, distances, H, unique_bbox, ...) {
    N <- diff(unique_bbox[1, ]) * diff(unique_bbox[2, ])
    output <- vector(mode = 'numeric', length = length(distances))
    mat_dist <- rgeos::gDistance(obj_sp1, obj_sp2,
                                 byid = TRUE,
                                 hausdorff = TRUE)
    output <- calc_h(x = mat_dist, dists = distances)
    output <- output * N
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
h_area <- function(obj_sp1, obj_sp2, distances, H, unique_bbox, ...) {
    x <- sf::st_sfc(obj_sp1)
    y <- sf::st_sfc(obj_sp2)

    if(is.null(unique_bbox))
        bb_xy <- sf::st_bbox(x)
    else
        bb_xy <- as.vector(unique_bbox)

    ## areas
    a_d <- diff(bb_xy[c(1, 3)])*diff(bb_xy[c(2, 4)])
    a_x <- sf::st_area(sf::st_union(x))
    a_y <- sf::st_area(sf::st_union(y))
    ## intensities
    l_x <- a_x/a_d
    l_y <- a_y/a_d
    ## number of elements
    n_x <- length(x)
    n_y <- length(y)

    k12 <- vector(mode = "numeric", length = length(distances))
    k21 <- vector(mode = "numeric", length = length(distances))

    k12 <- ((l_y*n_x)^(-1))*b_12(x, y, distances)
    k21 <- ((l_x*n_y)^(-1))*b_12(x, y, distances)

    out <- (n_x*k12 + n_y*k12)/(n_x + n_y)
    
    if(H == 'L') {
        out <- sqrt(out/pi)
    }

    return(out)
}


#' \eqn{H_{12}(d)}
#'
#' @description Evaluate a given function in the inputted distances.
#'
#' @param obj_sp1 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
#' @param obj_sp2 an object from class \code{SpatialPolygons} or \code{SpatialPointsDataFrame}
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
h_func <- function(obj_sp1, obj_sp2,
                   unique_bbox = NULL,
                   distances = NULL,
                   method = 'hausdorff', H = 'L', ...) {
    if((! "SpatialPolygons" %in% class(obj_sp1)) | (! "SpatialPolygons" %in% class(obj_sp2)))
        stop("obj_sp1 and obj_sp2 must be from class 'SpatialPolygons'")

    if(length(H) > 1 | length(method) > 1)
        stop(' method and H should have length == 1L.')

    if(! method %in% c('hausdorff', 'euclidean', 'area'))
        stop("method must be 'hausdorff', 'euclidean' or 'area'.")

    if(! H %in% c('K', 'L'))
        stop("H must be 'K' or 'L'.")

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
        max_d <- sqrt(diff(unique_bbox[1, ])^2 +
                      diff(unique_bbox[2, ])^2)
        distances <- seq(from = .05*max_d, to = .2*max_d,
                         length.out = 15L)
    }

    return(switch(method,
                  "hausdorff" = {
                      h_haus(obj_sp1, obj_sp2, distances, H, unique_bbox, ...)
                  },
                  "euclidean" = {
                      h_euc(obj_sp1, obj_sp2, distances, H, unique_bbox, ...)
                  },
                  "area" = {
                      h_area(obj_sp1, obj_sp2, distances, H, unique_bbox, ...)
                  }))
}

#' Fix distance matrix containing broken polygons
#'
#' @description fix a polygons' distance matrix based on a given method.
#'
#' @param x distance matrix
#' @param method method used to fix. The options are "min", "max", "mean", "rnd_poly", "rnd_dist", "min_norm", "max_norm", "hybrid", "hyb_center", "hybrid_nc", "old_min"
#' @return a distance matrix
#' @importFrom stats median
fix_dist <- function(x, method = "rnd_poly") {
    ## cleaning row names
    rownames(x) <- trimws(rownames(x))
    row_nms <- gsub("t", "", rownames(x))
    ## verifying if there are broken polygons
    if(method == "nothing") {
        x_out <- x
    } else {
        
    if(any(duplicated(row_nms))) {
        ## extracting duplicated rows (broken polys)
        extr <- unique(row_nms[duplicated(row_nms)])
        ids  <- which(row_nms %in% extr)

        ## auxiliar matrix without broken polygons
        x2 <- x[-ids, ]

        ## auxiliar matrix that will store the "right" distances
        x_extr <- matrix(nrow = length(extr), ncol = ncol(x2))
        rownames(x_extr) <- extr
        colnames(x_extr) <- colnames(x)

        ## auxiliar vector to filter the rows of the matrix
        ids_aux <- vector(mode = "list", length = length(extr))

        if(method == "min") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, min, na.rm = TRUE)
                }
            } 
        } else if(method == "old_min") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, min, na.rm = TRUE)
            } 
        } else if(method == "max") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, max, na.rm = TRUE)
                }
            }
        } else if(method == "mean") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, mean, na.rm = TRUE)
                }
            }
        } else if(method == "rnd_poly") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    x_extr[i, ] <- x[sample(ids_aux[[i]], size = 1), ]
                }
            }
        } else if(method == "rnd_dist") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, sample, size = 1)
                }
            }
        } else if(method == "min_norm") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    norm_aux <- apply(x[ids_aux[[i]], ], 1, function(x) sqrt(sum(x^2)))
                    x_extr[i, ] <- x[ids_aux[[i]][which.min(norm_aux)], ]
                }
            }
        } else if(method == "max_norm") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    norm_aux <- apply(x[ids_aux[[i]], ], 1, function(x) sqrt(sum(x^2)))
                    x_extr[i, ] <- x[ids_aux[[i]][which.max(norm_aux)], ]
                }
            }
        } else if(method == "hybrid") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    if(nrow(x[ids_aux[[i]], ]) > 2) {
                        norm_aux <- apply(x[ids_aux[[i]], ], 1, function(x) sqrt(sum(x^2)))
                        x_extr[i, ] <- x[ids_aux[[i]][which.min(norm_aux)], ]
                    } else {
                        x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, min, na.rm = TRUE)
                    }
                }
            }
        } else if(method == "hyb_center") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    if(nrow(x[ids_aux[[i]], ]) > 2) {
                        x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, median, na.rm = TRUE)
                    } else {
                        x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, mean, na.rm = TRUE)
                    }
                }
            }
        } else if(method == "hybrid_nc") {
            for(i in seq_along(extr)) {
                ids_aux[[i]] <- which(row_nms %in% extr[i])
                if(any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
                    ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
                    x_extr[i, ] <- x[ids_aux[[i]], ]
                } else {
                    if(nrow(x[ids_aux[[i]], ]) > 2) {
                        norm_aux <- apply(x[ids_aux[[i]], ], 1, function(x) sqrt(sum(x^2)))
                        x_extr[i, ] <- x[ids_aux[[i]][
                            which.min(abs(norm_aux - median(norm_aux)))],
                            ]
                    } else {
                        x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, mean, na.rm = TRUE)
                    }
                }
            }
        } else {
            stop("Provide a valid method!")
        }
        x_out <- rbind(x2, x_extr)
    } else {
        x_out <- x
    }
    }
    
    return(x_out)
}

#' Calculates the psam
#'
#' @description Computes the PSAM based on a distance matrix based on a method
#' to deal with broken polygons (represented by duplicated rows)
#'
#' @param x a distance matrix
#' @param method method to deal with broken polygons
#' @return scalar psam
calc_psam <- function(x, method = "rnd_poly") {
    x <- fix_dist(x, method = method)
    min_row <- apply(x, 1, min)
    min_col <- apply(x, 2, min)
    (sum(min_row) + sum(min_col))/(length(min_col) + length(min_row))
}

#' \eqn{h_{12}(t)}
#' 
#' @description Computes the \eqn{h_{12}} (K or L) based on a distance matrix based on a method
#'
#' @param x distance matrix
#' @param method method to deal with broken polygons
#' @param L logical scalar indicating if the L function should be used instead
#' @param dists vector of distances to compute \eqn{h_{12}(t)}.
#' @return a numeric vector
calc_h <- function(x, method = "min",
                   L = FALSE,
                   dists = NULL) {
    x <- fix_dist(x, method = method)
    n <- nrow(x)
    m <- ncol(x)

    if(is.null(dists)) dists <- seq(from = .05*sqrt(2),
                                    to   = .2*sqrt(2),
                                    length.out = 15L)
    
    h12 <- vector(mode = "numeric", length = length(dists))

    for(i in seq_along(h12)) {
        h12[i] <- sum(x <= dists[i])/(n*m)
    }

    if(L)
        return(sqrt(h12/pi))
    else
        return(h12)
}

