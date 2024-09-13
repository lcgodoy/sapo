##' Fix distance matrix containing broken polygons
##'
##' @description fix a polygons' distance matrix based on a given method. This
##'   function assumes the polygon that has been broken is represented by the
##'   rows of the distance matrix.
##'
##' @param x distance matrix
##' @param method method used to fix. The options are "min", "max", "mean",
##'   "rnd_poly", "rnd_dist", "min_norm", "max_norm", "hybrid", "hyb_center",
##'   "hybrid_nc", "old_min"
##'
##' @return a distance matrix
fix_dist <- function(x, method = "rnd_poly") {
  ## cleaning row names
  rownames(x) <- trimws(rownames(x))
  row_nms <- gsub("t", "", rownames(x))
  ## verifying if there are broken polygons
  if (method == "nothing") {
    x_out <- x
  } else {
    if (any(duplicated(row_nms))) {
      ## extracting duplicated rows (broken polys)
      extr <- unique(row_nms[duplicated(row_nms)])
      ids <- which(row_nms %in% extr)
      ## auxiliar matrix without broken polygons
      x2 <- x[-ids, ]
      ## auxiliar matrix that will store the "right" distances
      x_extr <- matrix(nrow = length(extr), ncol = ncol(x2))
      rownames(x_extr) <- extr
      colnames(x_extr) <- colnames(x)
      ## auxiliar vector to filter the rows of the matrix
      ids_aux <- vector(mode = "list", length = length(extr))
      if (method == "min") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, min, na.rm = TRUE)
          }
        }
      } else if (method == "old_min") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, min, na.rm = TRUE)
        }
      } else if (method == "max") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, max, na.rm = TRUE)
          }
        }
      } else if (method == "mean") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, mean, na.rm = TRUE)
          }
        }
      } else if (method == "rnd_poly") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            x_extr[i, ] <- x[sample(ids_aux[[i]], size = 1), ]
          }
        }
      } else if (method == "rnd_dist") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, sample, size = 1)
          }
        }
      } else if (method == "min_norm") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            norm_aux <- apply(x[ids_aux[[i]], ], 1, function(x) sqrt(sum(x^2)))
            x_extr[i, ] <- x[ids_aux[[i]][which.min(norm_aux)], ]
          }
        }
      } else if (method == "max_norm") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            norm_aux <- apply(x[ids_aux[[i]], ], 1, function(x) sqrt(sum(x^2)))
            x_extr[i, ] <- x[ids_aux[[i]][which.max(norm_aux)], ]
          }
        }
      } else if (method == "hybrid") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            if (nrow(x[ids_aux[[i]], ]) > 2) {
              norm_aux <- apply(
                x[ids_aux[[i]], ], 1,
                function(x) sqrt(sum(x^2))
              )
              x_extr[i, ] <- x[ids_aux[[i]][which.min(norm_aux)], ]
            } else {
              x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, min, na.rm = TRUE)
            }
          }
        }
      } else if (method == "hyb_center") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            if (nrow(x[ids_aux[[i]], ]) > 2) {
              x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, median, na.rm = TRUE)
            } else {
              x_extr[i, ] <- apply(x[ids_aux[[i]], ], 2, mean, na.rm = TRUE)
            }
          }
        }
      } else if (method == "hybrid_nc") {
        for (i in seq_along(extr)) {
          ids_aux[[i]] <- which(row_nms %in% extr[i])
          if (any(grepl("t", rownames(x)[ids_aux[[i]]]))) {
            ids_aux[[i]] <- ids_aux[[i]][grepl("t", rownames(x)[ids_aux[[i]]])]
            x_extr[i, ] <- x[ids_aux[[i]], ]
          } else {
            if (nrow(x[ids_aux[[i]], ]) > 2) {
              norm_aux <- apply(
                x[ids_aux[[i]], ], 1,
                function(x) sqrt(sum(x^2))
              )
              x_extr[i, ] <-
                x[ids_aux[[i]][which.min(abs(norm_aux - median(norm_aux)))], ]
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

##' Distance between polygons accounting for toroidal shift.
##'
##' @title ID aware distance matrix
##' @param p1 a \code{sf} object containing one column specifying the objects
##'   id.
##' @param p2 a \code{sf} object containing one column specifying the objects
##'   id.
##' @param hausdorff \code{logical} scalar indicating whether the Hausdorff
##'   distance should be used.
##' @param method method for "fixing" the distance matrix.
##' @return a distance matrix.
##' @author Lucas Godoy
iadist <- function(p1, p2, hausdorff = TRUE,
                   method = "rnd_poly") {
  stopifnot(sf::st_crs(p1) == sf::st_crs(p2))
  stopifnot("id" %in% intersect(names(p1), names(p2)))
  ll <- isTRUE(sf::st_is_longlat(p1))
  dtype <- ifelse(hausdorff, "Hausdorff",
    ifelse(ll, "Great Circle",
      "Eucliean"
    )
  )
  out <- sf::st_distance(p1, p2,
    which = dtype
  )
  rownames(out) <- p1[["id"]]
  colnames(out) <- p2[["id"]]
  out <- fix_dist(out, method)
  return(out)
}

##' @title \eqn{h_{12}(t)} from matrix
##'
##' @description Computes the \eqn{h_{12}} (K or L) based on a distance matrix
##'   based on a method
##'
##' @param x distance matrix
##' @param var_st logical scalar indicating if the L function should be used
##'   instead
##' @param dists vector of distances to compute \eqn{h_{12}(t)}.
##' @return a numeric vector
calc_h <- function(x,
                   var_st = FALSE,
                   dists = NULL) {
  n <- nrow(x)
  m <- ncol(x)
  if (is.null(dists)) {
    dd <- max(x)
    dists <- seq(
      from = .05 * dd,
      to = .2 * dd,
      length.out = 15L
    )
  }
  h12 <- vector(mode = "numeric", length = length(dists))
  for (i in seq_along(h12)) {
    h12[i] <- sum(x <= dists[i]) / (n * m)
  }
  if (var_st) {
    return(sqrt(h12 / pi))
  } else {
    return(h12)
  }
}

##' @title \eqn{h_{12}(t)} from polygons
##'
##' @description Computes the \eqn{h_{12}} (K or L) based on a distance matrix
##'   based on a method
##'
##' @param p1 sf object
##' @param p2 sf object
##' @param x a list with two sf objects.
##' @param hausdorff logical parameter indicating whether the Hausdorff distance
##'   should be used
##' @param method method to deal with broken polygons
##' @param var_st logical scalar indicating if the L function should be used
##'   instead
##' @param dists vector of distances to compute \eqn{h_{12}(t)}.
##' @name hfun
##' @return a numeric vector
h_func <- function(p1, p2,
                           hausdorff = TRUE,
                           method = "rnd_poly",
                           var_st = FALSE,
                           dists = NULL) {
  output <- iadist(p1, p2, hausdorff, method) |>
    calc_h(var_st, dists)
  return(output)
}

##' @name hfun
h_func.list <- function(x, ...) {
  output <- h_func(x[[1]], x[[2]], ...)
  return(output)
}
