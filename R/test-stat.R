##' aux function to calculate the mean of a vector when removing each of its
##' elements one by one.
##'
##' @title auxiliary mean
##' @param x a numeric vector
##' @return a numeric vector
##' @author Lucas Godoy
mean_aux <- function(x) {
  sapply(seq_along(x), function(id, x) mean(x[-id]),
         x = x)
}


##' Maximum Absolute Deviation
##'
##' @param x \code{numeric matrix}
##'
##' @return \code{numeric vector}
mad <- function(x) {
  out <- vector(mode = "numeric", length = nrow(x))
  mu <- apply(x, 2, mean_aux)
  for (i in seq_len(nrow(x))) {
    out[i] <- max(abs(x[i, ] - mu[i, ]))
  }
  return(out)
}

##' Studentized Maximum Absolute Deviation
##'
##' @param x \code{numeric matrix}
##'
##' @return \code{numeric vector}
s_mad <- function(x) {
  out <- vector(mode = "numeric", length = nrow(x))
  mu <- apply(x, 2, mean_aux)
  sd_mad <- apply(x[-nrow(x), ], 2, stats::sd)
  for (i in seq_len(nrow(x))) {
    out[i] <- max(abs(x[i, ] - mu[i, ]) / sd_mad)
  }
  return(out)
}

##' Maximum Absolute Deviation with Assimetry Correction
##'
##' @param x \code{numeric matrix}
##'
##' @return \code{numeric vector}
mad_ac <- function(x) {
  out <- vector(mode = "numeric", length = nrow(x))
  mu <- apply(x, 2, mean_aux)
  quant <- apply(x[-nrow(x), ], 2, quantile, c(.025, .975))
  for (i in seq_len(nrow(x))) {
    Ind <- as.numeric(x[i, ] >= mu[i, ])
    out[i] <- max(Ind * abs(x[i, ] - mu[i, ]) / abs(quant[1, ] - mu[i, ]) +
      (1 - Ind) * abs(x[i, ] - mu[i, ]) / abs(quant[2, ] - mu[i, ]))
  }
  return(out)
}

##' Integram Measure
##'
##' @param x \code{numeric matrix}
##' @param h \code{numeric}
##'
##' @return \code{numeric vector}
im <- function(x, h = 1) {
  out <- vector(mode = "numeric", length = nrow(x))
  if (ncol(x) %% 2 != 0) {
    x <- x[, -ncol(x)]
  }
  aux <- apply(x, 2, mean_aux)
  for (i in seq_len(nrow(x))) {
    out[i] <-
      (h / 3) * ((x[i, 1] - aux[i, 1])^2 +
        4 * sum((x[i, seq(2, ncol(x) - 2, 2)] -
          aux[i, seq(2, ncol(x) - 2, 2)])^2) +
        2 * sum((x[i, seq(3, ncol(x) - 1, 2)] -
          x[i, seq(3, ncol(x) - 1, 2)])^2) +
        (x[i, ncol(x)] - aux[i, ncol(x)])^2)
  }
  return(out)
}

##' Studentized Integram Measure
##'
##' @param x \code{numeric matrix}
##' @param h \code{numeric}
##'
##' @return \code{numeric vector}
s_im <- function(x, h = 1) {
  out <- vector(mode = "numeric", length = nrow(x))
  aux <- apply(x, 2, mean_aux)
  sd_gof <- apply(x[-nrow(x), ], 2, stats::sd)
  for (i in seq_len(nrow(x))) {
    out[i] <- sum(((x[i, ] - aux[i, ])^2) / sd_gof) * h
  }
  return(out)
}

##' Integram Measure with Assimetry Correction
##'
##' @param x \code{numeric matrix}
##' @param h \code{numeric}
##'
##' @return \code{numeric vector}
im_ac <- function(x, h = 1) {
  out <- vector(mode = "numeric", length = nrow(x))
  mu <- apply(x, 2, mean_aux)
  quant <- apply(x[-nrow(x), ], 2, quantile, c(.025, .975))
  for (i in seq_len(nrow(x))) {
    Ind <- as.numeric(x[i, ] >= mu[i, ])
    out[i] <- sum(Ind * ((x[i, ] - mu[i, ])^2) *
      h / abs(quant[1, ] - mu[i, ]) +
      (1 - Ind) * (((x[i, ] - mu[i, ])^2)) *
        h / abs(quant[2, ] - mu[i, ]))
  }
  return(out)
}
