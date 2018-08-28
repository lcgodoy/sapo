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

# Create copies of a set of polygons
#
# @param obj_sp object from class \code{SpatialPolygons}
# @param bbox_tot Boundary box from class \code{matrix}
#
# @description Auxiliar function.
#
# @importFrom stats runif
# @return an object from class \code{SpatialPolygons}
#
# poly_shift <- function(obj_sp, bbox_tot = NULL) {
#   if(!is.null(bbox_tot)) attr(obj_sp, "bbox") <- bbox_tot
#
#   range_x <- (max(obj_sp@bbox[1,]) - min(obj_sp@bbox[1,]))
#   range_y <- (max(obj_sp@bbox[2,]) - min(obj_sp@bbox[2,]))
#
#   n_poly <- length(obj_sp)
#
#   if(class(obj_sp) %in% "SpatialPolygons") {
#     # Obj 1 = objsp
#
#     # Obj 2 = copy 1 - translation in x
#
#     obj_sp2 <- obj_sp
#
#     aux <- obj_sp2@bbox
#
#     # X Shift
#
#     aux[1,] <- aux[1,] + range_x
#
#     # Defining shifted bbox
#
#     attr(obj_sp2, "bbox") <- aux
#
#     rm(aux)
#
#     # Shifting each polygon
#
#     #index_aux <- as.numeric(row.names(obj2sp))
#
#     for(i in 1:length(obj_sp2)) {
#       npolyaux <- obj_sp2@polygons[[i]]@Polygons %>% length
#       if(!npolyaux > 1) {
#         aux <- obj_sp2@polygons[[i]]@Polygons[[1]]@coords
#         aux[,1] <- aux[,1] + range_x
#         attr(obj_sp2@polygons[[i]]@Polygons[[1]], "coords") <- aux
#       } else {
#         for(k in seq_len(npolyaux)) {
#           aux <- obj_sp2@polygons[[i]]@Polygons[[k]]@coords
#           aux[,1] <- aux[,1] + range_x
#           attr(obj_sp2@polygons[[i]]@Polygons[[k]], "coords") <- aux
#
#         }
#       }
#     }
#
#     rm(aux)
#
#     # Obj 3 = copy 2 - translation in y
#
#     obj_sp3 <- obj_sp
#
#     aux <- obj_sp3@bbox
#
#     # Y Shift
#
#     aux[2,] <- aux[2,] + range_y
#
#     # Defining shifted bbox
#
#     attr(obj_sp3, "bbox") <- aux
#
#     rm(aux)
#
#     # Shifting each polygon
#
#     for(i in 1:length(obj_sp3)) {
#       npolyaux <- obj_sp3@polygons[[i]]@Polygons %>% length
#       if(!npolyaux > 1) {
#         aux <- obj_sp3@polygons[[i]]@Polygons[[1]]@coords
#         aux[,2] <- aux[,2] + range_y
#         attr(obj_sp3@polygons[[i]]@Polygons[[1]], "coords") <- aux
#       } else {
#         for(k in seq_len(npolyaux)) {
#           aux <- obj_sp3@polygons[[i]]@Polygons[[k]]@coords
#           aux[,2] <- aux[,2] + range_y
#           attr(obj_sp3@polygons[[i]]@Polygons[[k]], "coords") <- aux
#         }
#       }
#     }
#
#     rm(aux)
#
#     # Obj 4 = copy 3 - translation in y
#
#     obj_sp4 <- obj_sp
#
#     aux <- obj_sp4@bbox
#
#     # X and Y Shift
#
#     aux[1,] <- aux[1,] + range_x
#     aux[2,] <- aux[2,] + range_y
#
#     # Defining shifted bbox
#
#     attr(obj_sp4, "bbox") <- aux
#
#     rm(aux)
#
#     # Shifting each polygon
#
#     for(i in 1:length(obj_sp4)) {
#       npolyaux <- obj_sp4@polygons[[i]]@Polygons %>% length
#       if(!npolyaux > 1) {
#         aux <- obj_sp4@polygons[[i]]@Polygons[[1]]@coords
#         aux[,1] <- aux[,1] + range_x
#         aux[,2] <- aux[,2] + range_y
#         attr(obj_sp4@polygons[[i]]@Polygons[[1]], "coords") <- aux
#       } else {
#         for(k in seq_len(npolyaux)) {
#           aux <- obj_sp4@polygons[[i]]@Polygons[[k]]@coords
#           aux[,1] <- aux[,1] + range_x
#           aux[,2] <- aux[,2] + range_y
#           attr(obj_sp4@polygons[[i]]@Polygons[[k]], "coords") <- aux
#         }
#       }
#     }
#
#     output <- shift_aux(obj_sp, obj_sp2, obj_sp3, obj_sp4,
#                         length(obj_sp), range_x, range_y)
#   } else {
#
#     # Obj 1 = objsp
#
#     # Obj 2 = copy 1 - translation in x
#
#     obj_sp2 <- obj_sp
#
#     aux <- obj_sp2@bbox
#
#     # X Shift
#
#     aux[1,] <- aux[1,] + range_x
#
#     # Defining shifted bbox
#
#     attr(obj_sp2, "bbox") <- aux
#
#     rm(aux)
#
#     # Shifting each polygon
#
#     #index_aux <- as.numeric(row.names(obj2sp))
#
#     aux <- obj_sp2@coords
#     aux[, 1] <- aux[, 1] + range_x
#     attr(obj_sp2, "coords") <- aux
#
#
#     rm(aux)
#
#     # Obj 3 = copy 2 - translation in y
#
#     obj_sp3 <- obj_sp
#
#     aux <- obj_sp3@bbox
#
#     # Y Shift
#
#     aux[2,] <- aux[2,] + range_y
#
#     # Defining shifted bbox
#
#     attr(obj_sp3, "bbox") <- aux
#
#     rm(aux)
#
#     # Shifting each polygon
#
#     aux <- obj_sp3@coords
#     aux[, 2] <- aux[, 2] + range_y
#     attr(obj_sp3, 'coords') <- aux
#
#     rm(aux)
#
#     # Obj 4 = copy 3 - translation in y
#
#     obj_sp4 <- obj_sp
#
#     aux <- obj_sp4@bbox
#
#     # X and Y Shift
#
#     aux[1,] <- aux[1,] + range_x
#     aux[2,] <- aux[2,] + range_y
#
#     # Defining shifted bbox
#
#     attr(obj_sp4, "bbox") <- aux
#
#     rm(aux)
#
#     # Shifting each polygon
#
#     aux <- obj_sp4@coords
#     aux[,1] <- aux[,1] + range_x
#     aux[,2] <- aux[,2] + range_y
#     attr(obj_sp4, "coords") <- aux
#
#     output <- shift_aux_point(obj_sp, obj_sp2, obj_sp3, obj_sp4)
#   }
#
#   rm(list = ls()[ls() != "output"])
#
#   return(output)
# }

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
  if(m == 'haus2') {
    return(sp_ID_haus2(obj_sp1 = obj_sp1,
                       obj_sp2 = obj_sp2))
  }


}

#' Hausdorff Distance Matrix - ID
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

  matrix_hauss <- matrix(data = rep(0, n_poly1*n_poly2), nrow = n_poly2,
                         ncol = n_poly1, dimnames = list(unique(nameaux2),
                                                         unique(nameaux1))
  )
  matrix_aux <- gDistance(obj_sp1, obj_sp2, byid = T)
  matrix_aux_hauss <- gDistance(obj_sp1, obj_sp2, byid = T, hausdorff = T)

  matrix_out <- data.table::melt(matrix_aux, value.name = 'euc') %>%
    cbind(hauss = data.table::melt(matrix_aux_hauss, value.name = 'hauss')$hauss)

  for(i in seq_len(nrow(matrix_hauss))) {
    for(j in seq_len(ncol(matrix_hauss))) {
      cond1 <- matrix_out$Var1 == rownames(matrix_hauss)[i] & matrix_out$Var2 == colnames(matrix_hauss)[j]
      k <- which(cond1 & matrix_out$euc == min(matrix_out$euc[cond1]))
      matrix_hauss[i, j] <- matrix_out$hauss[k]
    }
  }

  rm(list = ls()[ls() != "matrix_hauss"])

  return(matrix_hauss)

}

#' Hausdorff Distance Matrix - ID (Method 2)
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
sp_ID_haus2 <- function(obj_sp1, obj_sp2) {
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
