#' Polygons' Random Shift - 2
#'
#' @param obj_sp object from class \code{SpatialPolygons}
#' @param bbox_max Boundary box from class \code{matrix}
#'
#' @importFrom stats runif
#'
#' @return an object from class \code{SpatialPolygons} randomly translated
#'
poly_rf2 <- function(obj_sp, bbox_max) {
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

  for(i in 1:length(obj_sp)) {
    npolyaux <- obj_sp@polygons[[i]]@Polygons %>% length
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

  return(obj_sp)

}

#' Points' Random Shift
#'
#' @param obj_sp object from class \code{SpatialPoints}
#' @param bbox_max Boundary box from class \code{matrix}
#'
#' @importFrom stats runif
#'
#' @return an object from class \code{SpatialPoints} randomly translated
#'
pt_rf <- function(obj_sp, bbox_max) {
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

  aux <- obj_sp@coords
  aux[,1] <- aux[,1] + jump_x
  aux[,2] <- aux[,2] + jump_y
  attr(obj_sp, "coords") <- aux

  return(obj_sp)

}

#' psat_mc2 - auxiliar function
#'
#' @description Internal use.
#'
#' @param obj1_shift object from class \code{SpatialPolygons} or \code{SpatialPoints}
#' @param obj_sp2 object from class \code{SpatialPolygons} or \code{SpatialPoints}
#' @param niter \code{numeric} value indicating number of iterations
#' @param ts \code{string} giving the test statistic, the options are \code{c('psam', 'pf12', 'pk12')}
#' @param args \code{list} only necessary when the test statistic is different from psam
#'
#' @return a \code{numeric vector}.
#'
mc_iterations <- function(obj1_shift, obj_sp2, niter, ts, args = NULL) {

  if(ts == 'psam') {
    output <- vector(mode = 'numeric', length = niter)

    if('SpatialPolygons' %in% class(obj1_shift)) {
      if('SpatialPolygons' %in% class(obj_sp2)) {
        for(i in seq_len(niter)) {
          obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          obj1_aux <- gIntersection(obj1_shift, limits_to_sp(obj2_rshift@bbox), byid = T,
                                    id = suppressWarnings(names(obj1_aux)))
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[i] <- psam(obj1_aux, obj2_rshift)
        }
      } else {
        for(i in seq_len(niter)) {
          obj2_rshift <- pt_rf(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          obj1_aux <- gIntersection(obj1_shift, limits_to_sp(obj2_rshift@bbox), byid = T,
                                    id = suppressWarnings(names(obj1_aux)))
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[i] <- psam(obj1_aux, obj2_rshift)
        }
      }
    } else {
      if('SpatialPolygons' %in% class(obj_sp2)) {
        for(i in seq_len(niter)) {
          obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          k <- which(rgeos::gWithin(obj1_aux, limits_to_sp(obj2_rshift@bbox), byid = T))
          obj1_aux <- obj1_aux[k, ]
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[i] <- psam(obj1_aux, obj2_rshift)
        }
      } else {
        for(i in seq_len(niter)) {
          obj2_rshift <- pt_rf(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          k <- which(rgeos::gWithin(obj1_aux, limits_to_sp(obj2_rshift@bbox), byid = T))
          obj1_aux <- obj1_aux[k, ]
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[i] <- psam(obj1_aux, obj2_rshift)
        }
      }
    }
  }

  if(ts == 'pf12') {

    output <- vector(mode = 'list', length = niter)

    if('SpatialPolygons' %in% class(obj1_shift)) {
      if('SpatialPolygons' %in% class(obj_sp2)) {
        for(i in seq_len(niter)) {
          obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          obj1_aux <- gIntersection(obj1_shift, limits_to_sp(obj2_rshift@bbox), byid = T,
                                    id = suppressWarnings(names(obj1_aux)))
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[[i]] <- pf12(obj1_aux, obj2_rshift,
                              r_min = args$r_min, r_max = args$r_max)
        }
      } else {
        for(i in seq_len(niter)) {
          obj2_rshift <- pt_rf(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          obj1_aux <- gIntersection(obj1_shift, limits_to_sp(obj2_rshift@bbox), byid = T,
                                    id = suppressWarnings(names(obj1_aux)))
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[[i]] <- pf12(obj1_aux, obj2_rshift,
                              r_min = args$r_min, r_max = args$r_max)
        }
      }
    } else {
      if('SpatialPolygons' %in% class(obj_sp2)) {
        for(i in seq_len(niter)) {
          obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          k <- which(rgeos::gWithin(obj1_aux, limits_to_sp(obj2_rshift@bbox), byid = T))
          obj1_aux <- obj1_aux[k, ]
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[[i]] <- pf12(obj1_aux, obj2_rshift,
                              r_min = args$r_min, r_max = args$r_max)
        }
      } else {
        for(i in seq_len(niter)) {
          obj2_rshift <- pt_rf(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          k <- which(rgeos::gWithin(obj1_aux, limits_to_sp(obj2_rshift@bbox), byid = T))
          obj1_aux <- obj1_aux[k, ]
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[[i]] <- pf12(obj1_aux, obj2_rshift,
                              r_min = args$r_min, r_max = args$r_max)
        }
      }
    }

    output <- do.call('rbind', output) %>% as.data.frame

  }

  if(ts == 'pk_dist12') {

    output <- vector(mode = 'list', length = niter)

    if('SpatialPolygons' %in% class(obj1_shift)) {
      if('SpatialPolygons' %in% class(obj_sp2)) {
        for(i in seq_len(niter)) {
          obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          obj1_aux <- gIntersection(obj1_shift, limits_to_sp(obj2_rshift@bbox), byid = T,
                                    id = suppressWarnings(names(obj1_aux)))
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[[i]] <- pk_dist12(obj1_aux, obj2_rshift,
                                   r_min = args$r_min, r_max = args$r_max,
                                   bbox = obj1_aux@bbox)
        }
      }
    }

    output <- do.call('rbind', output) %>% as.data.frame

  }

  if(ts == 'pk_area12') {

    output <- vector(mode = 'list', length = niter)

    if('SpatialPolygons' %in% class(obj1_shift)) {
      if('SpatialPolygons' %in% class(obj_sp2)) {
        for(i in seq_len(niter)) {
          obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
          obj1_aux <- obj1_shift
          obj1_aux <- gIntersection(obj1_shift, limits_to_sp(obj2_rshift@bbox), byid = T,
                                    id = suppressWarnings(names(obj1_aux)))
          attr(obj1_aux, "bbox") <- obj2_rshift@bbox
          output[[i]] <- pk_area12(obj1_aux, obj2_rshift,
                                   r_min = args$r_min, r_max = args$r_max,
                                   bbox = obj1_aux@bbox)
        }
      }
    }

    output <- do.call('rbind', output) %>% as.data.frame

  }

  rm(list = ls()[ls() != 'output'])

  return(output)
}
