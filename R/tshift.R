##' @title Pre-TS
##'
##' @description Create rigid copies of a polygon. This function an auxilliary
##'   function for the Toroidal Shift method
##'
##' @param poly an object of class \code{sf} or \code{sfc}.
##' @param bb (optional) a unique bounding box.
##' @param id_col a \code{character} indicating the id column in \code{poly}.
##'
##' @return an \code{sf} with 8 additional rigid copies of \code{poly}.
##'
##' @author Lucas Godoy
pre_ts <- function(poly,
                   bb = NULL,
                   id_col = NULL) {
  if (methods::is(poly, "sfc")) {
    poly <- sf::st_as_sf(poly)
  }
  if (is.null(id_col)) {
    poly$id <- seq_len(NROW(poly))
  }
  if (is.null(bb)) {
    bb <- sf::st_bbox(poly)
  } else {
    stopifnot(methods::is(bb, "bbox"))
  }
  range_x <- diff(bb[c(1, 3)]) |>
    unname()
  range_y <- diff(bb[c(2, 4)]) |>
    unname()
  translations <-
    expand.grid(
      x = c(0, range_x, -range_x),
      y = c(0, range_y, -range_y),
      KEEP.OUT.ATTRS = FALSE
    )[-1, ] |>
    apply(1, sf::st_point, simplify = FALSE) |>
    unname() |>
    lapply(sf::st_sfc)
  shifted <- lapply(translations,
    translate_by_pt,
    poly = poly
  )
  shifted <- lapply(shifted, \(x) {
    x$orig <- 0
    return(x)
  })
  poly$orig <- 1
  shifted <-
    do.call(rbind, shifted) |>
    rbind(poly)
  shifted$old_id <- shifted$id
  shifted$id <- with(shifted, ifelse(orig == 1, paste0("t", id), id))
  return(shifted)
}

##' @title Create jumps for random movements
##'
##' @details This is an internal function.
##'
##' @param unique_bb a bbox shared between both "Polygon Patterns"
##'
##' @return a `sfc` object representing a random jump or shift.
##'
##' @author Lucas Godoy
create_jump <- function(unique_bb) {
  jump_size <- c(
    "x" = stats::runif(1,
      min = unique_bb[1],
      max = unique_bb[3]
    ),
    "y" = stats::runif(1,
      min = unique_bb[2],
      max = unique_bb[4]
    )
  )
  possible_jumps <-
    expand.grid(
      x = c(0, -jump_size["x"], jump_size["x"]),
      y = c(0, -jump_size["y"], jump_size["y"])
    )[-1, ]
  jump <-
    possible_jumps[sample(seq_len(NROW(possible_jumps)), 1), ] |>
    as.numeric() |>
    sf::st_point() |>
    sf::st_sfc()
  return(jump)
}

##' @title Translate an `sf` object by a "point"
##' @param pt `sfc` representing a shift.
##' @param poly `sfc` of `sf` to be shifted
##'
##' @return a `sf` or `sfc` representing `poly` shifted by `pt`
##'
##' @author Lucas Godoy
translate_by_pt <- function(pt, poly) {
  sf::st_agr(poly) <- "constant"
  sf::st_agr(pt) <- "constant"
  sf::st_geometry(poly) <-
    sf::st_geometry(poly) + pt
  return(poly)
}

##' @title Toroidal Shift
##'
##' @param x a `sf` or `sfc` object. Its `geometry` may contain POLYGONS and/or
##'   POINTS.
##' @param y a `sf` or `sfc` object. Its `geometry` may contain POLYGONS and/or
##'   POINTS.
##' @param shifted `logical` indicating whether `x` has been "shifted". This
##'   parameter is mainly for internal use and testing.
##' @param unique_bb a bbox shared between both "Polygon Patterns"
##'
##' @return a list
##'
##' @author Lucas Godoy
toroidal_shift <- function(x, y,
                           shifted = FALSE,
                           unique_bb = NULL) {
  stopifnot(
    methods::is(x, "sf") || methods::is(x, "sfc"),
    methods::is(y, "sf") || methods::is(y, "sfc"),
    sf::st_crs(x) == sf::st_crs(y)
  )
  if (is.null(unique_bb)) {
    if (shifted) {
      stop(shifted, "Unique bbox must be provided when `shifted = TRUE`")
    }
    sf::st_agr(x) <- "constant"
    sf::st_agr(y) <- "constant"
    unique_bb <-
      sf::st_intersection(
        sf::st_as_sfc(sf::st_bbox(x)),
        sf::st_as_sfc(sf::st_bbox(y))
        ) |>
      sf::st_bbox()
  }
  if (!shifted) {
    x <- pre_ts(x)
  }
  jump <- create_jump(unique_bb)
  y <- translate_by_pt(
    pt = jump,
    poly = y
  )
  new_bb <- sf::st_bbox(y)
  x <- sf::st_intersection(x, sf::st_as_sfc(new_bb))
  return(list(
    "rigid" = x,
    "random" = y
  ))
}
