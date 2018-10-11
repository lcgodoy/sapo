#' psat_mc - auxiliar function
#'
#' @description Internal use.
#'
#' @param obj1_shift object from class \code{SpatialPolygons} or \code{SpatialPoints}
#' @param obj_sp2 object from class \code{SpatialPolygons} or \code{SpatialPoints}
#' @param niter \code{numeric} value indicating number of iterations
#' @param ts \code{string} giving the test statistic, the options are \code{c('psam', 'pf12', 'pk12')}
#' @param correction a \code{character} giving the edge correction to be used. Possible
#' entries are \code{c('none', 'torus', 'guard', 'adjust')}.
#'
#' @param ... parameters for test statistics functions
#' @return a \code{numeric vector}.
#'
mc_iterations <- function(obj1_shift, obj_sp2, niter, ts, correction, ...) {

  if(ts == 'psam') {
    output <- vector(mode = 'numeric', length = niter)

    for(i in seq_len(niter)) {
      obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
      if(correction == 'torus') {
        output[i] <- psam(obj1_shift, obj2_rshift, correction = 'none', ...)
      } else {
        obj1_aux <- obj1_shift
        lm_sp <- limits_to_sp(obj2_rshift@bbox)
        lm_sp@proj4string <- obj2_rshift@proj4string
        obj1_aux <- gIntersection(obj1_shift, lm_sp, byid = T,
                                  id = suppressWarnings(names(obj1_aux)))
        attr(obj1_aux, "bbox") <- obj2_rshift@bbox
        output[i] <- psam(obj1_aux, obj2_rshift, correction = correction, ...)
      }
    }
  }

  if(ts == 'pk_dist12') {

    # output <- matrix(nrow = niter, ncol = 2)
    output <- vector(mode = 'list', length = niter)

    for(i in seq_len(niter)) {
      obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
      if(correction == 'torus') {
        output[[i]] <- pk_dist12(obj1_shift, obj2_rshift,
                                 bbox = obj1_shift@bbox,
                                 correction = 'none',
                                 ...)
      } else {
        obj1_aux <- obj1_shift
        lm_sp <- limits_to_sp(obj2_rshift@bbox)
        lm_sp@proj4string <- obj2_rshift@proj4string
        obj1_aux <- gIntersection(obj1_shift, lm_sp, byid = T,
                                  id = suppressWarnings(names(obj1_aux)))
        attr(obj1_aux, "bbox") <- obj2_rshift@bbox
        output[[i]] <- pk_dist12(obj1_aux, obj2_rshift,
                                 bbox = obj1_aux@bbox,
                                 correction = correction,
                                 ...)
      }
    }

    output <- data.table::rbindlist(output) %>% as.data.frame

  }

  if(ts == 'pk_area12') {

    output <- vector(mode = 'list', length = niter)

    for(i in seq_len(niter)) {
      obj2_rshift <- poly_rf2(obj_sp2, obj1_shift@bbox)
      if(correction == 'torus') {
        output[[i]] <- pk_area12(obj1_shift, obj2_rshift,
                                 bbox = obj1_shift@bbox,
                                 correction = 'none',
                                 ...)
      } else {
        obj1_aux <- obj1_shift
        lm_sp <- limits_to_sp(obj2_rshift@bbox)
        lm_sp@proj4string <- obj2_rshift@proj4string
        obj1_aux <- gIntersection(obj1_shift, lm_sp, byid = T,
                                  id = suppressWarnings(names(obj1_aux)))
        attr(obj1_aux, "bbox") <- obj2_rshift@bbox
        output[[i]] <- pk_area12(obj1_aux, obj2_rshift,
                                 bbox = obj1_aux@bbox,
                                 correction = correction,
                                 ...)
      }
    }

    output <- data.table::rbindlist(output) %>% as.data.frame

  }

  rm(list = ls()[ls() != 'output'])

  return(output)
}
