#' Polygons' Spatial Association Measure
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard')}.
#' @param hausdorff a \code{boolean}. If \code{TRUE}, then the Hausdorff distance is used.
#' Otherwise the Euclidean distance is used. Default is \code{FALSE}.
#'
#' @return a \code{numeric} \code{scalar} with the polygons' spatial
#'  association measure for the corresponding objects.
#' @export
#'
psam <- function(obj_sp1, obj_sp2, correction = 'none', hausdorff = FALSE) {
    mat_dist <- rgeos::gDistance(obj_sp1, obj_sp2,
                                 byid = TRUE,
                                 hausdorff = TRUE)

    calc_psam(mat_dist, method = "rnd_poly")
}

#' Maximum Absolute Deviation
#'
#' @param x \code{numeric matrix}
#'
#' @return \code{numeric vector}
mad <- function(x) {
    out <- vector(mode = 'numeric', length = nrow(x))
    mu <- apply(x, 2, mean_vec)
    for(i in seq_len(nrow(x))) {
        out[i] <- max(abs(x[i, ] - mu[i, ]))
    }
    return(out)
}

#' Studentized Maximum Absolute Deviation
#'
#' @param x \code{numeric matrix}
#'
#' @return \code{numeric vector}
s_mad <- function(x) {
    out <- vector(mode = 'numeric', length = nrow(x))
    mu <- apply(x, 2, mean_vec)
    sd_mad <- apply(x[-nrow(x),], 2, stats::sd)
    for(i in seq_len(nrow(x))) {
        out[i] <- max(abs(x[i, ] - mu[i, ])/sd_mad)
    }
    return(out)
}

#' Maximum Absolute Deviation with Assimetry Correction
#'
#' @param x \code{numeric matrix}
#'
#' @return \code{numeric vector}
mad_ac <- function(x) {
  out <- vector(mode = 'numeric', length = nrow(x))
  mu <- apply(x, 2, mean_vec)
  quant <- apply(x[-nrow(x),], 2, quantile, c(.025, .975))
  for(i in seq_len(nrow(x))) {
    Ind <- as.numeric(x[i, ] >= mu[i, ])
    out[i] <- max(Ind*abs(x[i, ] - mu[i, ])/abs(quant[1, ] - mu[i, ]) +
                    (1 - Ind)*abs(x[i, ] - mu[i, ])/abs(quant[2, ] - mu[i, ]))
  }
  return(out)
}

#' Integram Measure
#'
#' @param x \code{numeric matrix}
#' @param h \code{numeric}
#'
#' @return \code{numeric vector}
im <- function(x, h = 1) {
  out <- vector(mode = 'numeric', length = nrow(x))
  if(ncol(x) %% 2 != 0) {
    x <- x[, -ncol(x)]
  }
  aux <- apply(x, 2, mean_vec)
  for(i in seq_len(nrow(x))) {
    out[i] <- (h/3)*((x[i, 1] - aux[i, 1])^2 +
                       4*sum((x[i, seq(2, ncol(x) - 2, 2)] - aux[i, seq(2, ncol(x) - 2, 2)])^2) +
                       2*sum((x[i, seq(3, ncol(x) - 1, 2)] - x[i, seq(3, ncol(x) - 1, 2)])^2) +
                       (x[i, ncol(x)] - aux[i, ncol(x)])^2)
  }
  return(out)
}

#' Studentized Integram Measure
#'
#' @param x \code{numeric matrix}
#' @param h \code{numeric}
#'
#' @return \code{numeric vector}
s_im <- function(x, h = 1) {
  out <- vector(mode = 'numeric', length = nrow(x))
  aux <- apply(x, 2, mean_vec)
  sd_gof <- apply(x[-nrow(x),], 2, stats::sd)
  for(i in seq_len(nrow(x))) {
    out[i] <- sum(((x[i, ] - aux[i, ])^2)/sd_gof)*h
  }
  return(out)
}

#' Integram Measure with Assimetry Correction
#'
#' @param x \code{numeric matrix}
#' @param h \code{numeric}
#'
#' @return \code{numeric vector}
im_ac <- function(x, h = 1) {
  out <- vector(mode = 'numeric', length = nrow(x))
  mu <- apply(x, 2, mean_vec)
  quant <- apply(x[-nrow(x),], 2, quantile, c(.025, .975))
  for(i in seq_len(nrow(x))) {
    Ind <- as.numeric(x[i, ] >= mu[i, ])
    out[i] <- sum(Ind*((x[i, ] - mu[i, ])^2)*h/abs(quant[1, ] - mu[i, ]) +
                    (1 - Ind)*(((x[i, ] - mu[i, ])^2))*h/abs(quant[2, ] - mu[i, ]))
  }
  return(out)
}
