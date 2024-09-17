##' Polygons Spatial Association Test - Global Envelope
##'
##' @description A Monte Carlo test to verify if two sets of polygons are
##'    associated based in a global envelope of the functions
##'    \eqn{K_{12}(d)} and \eqn{L_{12}(d)} using different test statistics.
##'
##' @param p1 a \code{sf} object containing one column specifying the objects
##'   id.
##' @param p2 a \code{sf} object containing one column specifying the objects
##'   id.
##' @param id_col a \code{character} or \code{integer} indicating the column of
##'   \code{p1} storing the unique identifier for the polygons/sample units.
##' @param n_sim an \code{integer} corresponding to the number of Monte Carlo
##'   simulations for the test
##' @param alpha a \code{numeric} indicating the confidence level.
##' @param var_st use the variance stabilizing funciton?
##' @param ts a \code{character} associated to a test statistic. Inputs acepted:
##'   \code{c('IM', 'MAD', 'SIM', 'SMAD', 'IMDQ', 'MADDQ')}.
##' @param hausdorff a \code{logical} scalar indicating whether the Hausdorff
##'   distance should be used (default is TRUE).
##' @param distances a \code{numeric vector} indicating the distances to
##'   evaluate \eqn{H(d)}. If \code{NULL} then the range considered goes from 5%
##'   to 20% of the max distance that can be observed inside the study region.
##' @param method a \code{character} specifying which kind of distance will be
##'   used to evalueate \eqn{H}. Also, there is an option of using
##'   areas. Options available: \code{c("hausdorff", "euclidean")}.
##'
##'
##' @return a \code{list} with values: \describe{
##'     \item{p_value}{a \code{numeric} scalar giving the p-value of the test}
##'     \item{mc_sample}{a \code{numeric} vector giving the test statistic for each of the Monte Carlo simulations}
##'     \item{mc_funct}{a \code{matrix} where each line correspond to the function (\eqn{K} or \eqn{L}) estimated
##'     for the Monte Carlo simulations}
##'     \item{distances}{\code{numeric vector} containing the distances where mc_func were evaluated.}
##'     \item{alpha}{a \code{numeric} scalar giving the significance level}
##'     \item{rejects}{a \code{logical} scalar, TRUE if the null hypothesis is reject}
##'   }
##'
##' @export
##'
gof_mc <- function(p1, p2,
                   id_col = NULL,
                   n_sim = 499L,
                   alpha = 0.01,
                   var_st = TRUE,
                   ts = "SMAD",
                   distances = NULL,
                   hausdorff = TRUE,
                   method = "rnd_poly") {
  stopifnot("sf" %in% union(class(p1), class(p2)))
  stopifnot(length(n_sim) == 1)
  stopifnot(length(alpha) == 1)
  stopifnot(length(var_st) == 1)
  stopifnot(length(ts) == 1)
  stopifnot(length(hausdorff) == 1)
  stopifnot(length(method) == 1)
  stopifnot(method %in% c(
    "min", "max", "mean",
    "rnd_poly", "rnd_dist", "min_norm",
    "max_norm", "hybrid", "hyb_center",
    "hybrid_nc", "old_min"
  ))
  if (is.null(distances)) {
    max_d1 <- sf::st_bbox(p1)
    max_d1 <- diff(max_d1[c(1, 3)])^2 +
      diff(max_d1[c(2, 4)])^2
    max_d2 <- sf::st_bbox(p2)
    max_d2 <- diff(max_d2[c(1, 3)])^2 +
      diff(max_d2[c(2, 4)])^2
    max_d <- sqrt(max(c(max_d1, max_d2)))
    distances <- seq(
      from = .05 * max_d, to = .2 * max_d,
      length.out = 15L
    )
  }
  unique_bbox <- sf::st_intersection(
    sf::st_as_sfc(sf::st_bbox(p1)),
    sf::st_as_sfc(sf::st_bbox(p2))
    ) |>
    sf::st_bbox()
  output <- vector(mode = "list", length = 6L)
  names(output) <- c(
    "p_value", "rejects",
    "mc_sample", "mc_funct",
    "distances", "alpha"
  )
  output$mc_sample <- vector(mode = "numeric", length = n_sim + 1L)
  output$mc_funct <- matrix(ncol = length(distances), nrow = n_sim + 1L)
  output$alpha <- alpha
  output$rejects <- FALSE
  output$distances <- distances
  output$mc_funct[n_sim + 1L, ] <- h_func(
      p1, p2, hausdorff, method,
      var_st, dists = distances
  )
  p1_shifted <- pre_ts(p1, bb = unique_bbox, id_col = id_col)
  for (i in seq_len(n_sim))
    output$mc_funct[i, ] <- toroidal_shift(p1_shifted, p2,
                                           shifted = TRUE,
                                           unique_bb = unique_bbox) |>
      h_func.list(hausdorff, method, var_st, distances)
  if (length(distances) > 1) {
    h <- distances[length(distances)] - distances[length(distances) - 1L]
  } else {
    h <- 1
  }
  output$mc_sample <- switch(ts,
    "IM" = {
      im(x = output$mc_funct, h = h)
    },
    "MAD" = {
      mad(x = output$mc_funct)
    },
    "SIM" = {
      s_im(x = output$mc_funct, h = h)
    },
    "SMAD" = {
      s_mad(x = output$mc_funct)
    },
    "IMDQ" = {
      im_ac(x = output$mc_funct, h = h)
    },
    "MADDQ" = {
      mad_ac(x = output$mc_funct)
    }
  )
  output$p_value <- mean(output$mc_sample[n_sim + 1L] <= output$mc_sample)
  if (output$p_value <= output$alpha) output$rejects <- TRUE
  ## class(output) <- gof_test(output)
  return(output)
}
