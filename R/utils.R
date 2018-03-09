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

#' Create copies of a set of polygons
#'
#' @param obj_sp object from class \code{SpatialPolygons}
#' @param bbox_tot Boundary box from class \code{matrix}
#'
#' @description Auxiliar function.
#'
#' @importFrom stats runif
#' @return an object from class \code{SpatialPolygons}
#'
poly_shift <- function(obj_sp, bbox_tot = NULL) {
  if(!is.null(bbox_tot)) attr(obj_sp, "bbox") <- bbox_tot

  range_x <- (max(obj_sp@bbox[1,]) - min(obj_sp@bbox[1,]))
  range_y <- (max(obj_sp@bbox[2,]) - min(obj_sp@bbox[2,]))

  n_poly <- length(obj_sp)

  if(class(obj_sp) %in% "SpatialPolygons") {
    # Obj 1 = objsp

    # Obj 2 = copy 1 - translation in x

    obj_sp2 <- obj_sp

    aux <- obj_sp2@bbox

    # X Shift

    aux[1,] <- aux[1,] + range_x

    # Defining shifted bbox

    attr(obj_sp2, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    #index_aux <- as.numeric(row.names(obj2sp))

    for(i in 1:length(obj_sp2)) {
      aux <- obj_sp2@polygons[[i]]@Polygons[[1]]@coords
      aux[,1] <- aux[,1] + range_x
      attr(obj_sp2@polygons[[i]]@Polygons[[1]], "coords") <- aux
    }

    rm(aux)

    # Obj 3 = copy 2 - translation in y

    obj_sp3 <- obj_sp

    aux <- obj_sp3@bbox

    # Y Shift

    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp3, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    for(i in 1:length(obj_sp3)) {
      aux <- obj_sp3@polygons[[i]]@Polygons[[1]]@coords
      aux[,2] <- aux[,2] + range_y
      attr(obj_sp3@polygons[[i]]@Polygons[[1]], "coords") <- aux
    }

    rm(aux)

    # Obj 4 = copy 3 - translation in y

    obj_sp4 <- obj_sp

    aux <- obj_sp4@bbox

    # X and Y Shift

    aux[1,] <- aux[1,] + range_x
    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp4, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    for(i in 1:length(obj_sp4)) {
      aux <- obj_sp4@polygons[[i]]@Polygons[[1]]@coords
      aux[,1] <- aux[,1] + range_x
      aux[,2] <- aux[,2] + range_y
      attr(obj_sp4@polygons[[i]]@Polygons[[1]], "coords") <- aux
    }

    output <- shift_aux(obj_sp, obj_sp2, obj_sp3, obj_sp4,
                        length(obj_sp), range_x, range_y)
  } else {

    # Obj 1 = objsp

    # Obj 2 = copy 1 - translation in x

    obj_sp2 <- obj_sp

    aux <- obj_sp2@bbox

    # X Shift

    aux[1,] <- aux[1,] + range_x

    # Defining shifted bbox

    attr(obj_sp2, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    #index_aux <- as.numeric(row.names(obj2sp))

    aux <- obj_sp2@coords
    aux[, 1] <- aux[, 1] + range_x
    attr(obj_sp2, "coords") <- aux


    rm(aux)

    # Obj 3 = copy 2 - translation in y

    obj_sp3 <- obj_sp

    aux <- obj_sp3@bbox

    # Y Shift

    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp3, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    aux <- obj_sp3@coords
    aux[, 2] <- aux[, 2] + range_y
    attr(obj_sp3, 'coords') <- aux

    rm(aux)

    # Obj 4 = copy 3 - translation in y

    obj_sp4 <- obj_sp

    aux <- obj_sp4@bbox

    # X and Y Shift

    aux[1,] <- aux[1,] + range_x
    aux[2,] <- aux[2,] + range_y

    # Defining shifted bbox

    attr(obj_sp4, "bbox") <- aux

    rm(aux)

    # Shifting each polygon

    aux <- obj_sp4@coords
    aux[,1] <- aux[,1] + range_x
    aux[,2] <- aux[,2] + range_y
    attr(obj_sp4, "coords") <- aux

    output <- shift_aux_point(obj_sp, obj_sp2, obj_sp3, obj_sp4)
  }

  rm(list = ls()[ls() != "output"])

  return(output)
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
      aux <- obj_sp@polygons[[i]]@Polygons[[1]]@coords
      aux[,1] <- aux[,1] + jump_x
      aux[,2] <- aux[,2] + jump_y
      attr(obj_sp@polygons[[i]]@Polygons[[1]], "coords") <- aux
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
#'
#' @importFrom rgeos gDistance
#'
#' @return a distance matrix
#'
sp_ID_dist <- function(obj_sp1, obj_sp2) {
  nameaux <- suppressWarnings(row.names(obj_sp1))

  n_poly1 <- length(table(nameaux))

  matrix_dist <- matrix(data = rep(0, n_poly1*length(obj_sp2)), nrow = length(obj_sp2),
                        ncol = n_poly1)

  colnames(matrix_dist) <- names(table(nameaux))
  rownames(matrix_dist) <- suppressWarnings(row.names(obj_sp2))

  matrix_aux <- gDistance(obj_sp1, obj_sp2, byid = T)

  if(any(duplicated(colnames(matrix_aux)))) {
    for(i in 1:ncol(matrix_dist)) {
      aux <- matrix_aux[, which(colnames(matrix_aux) == colnames(matrix_dist)[i]), drop = F]
      matrix_dist[, i] <- apply(aux, 1, min)
      rm(aux)
    }
  }
  if(any(duplicated(rownames(matrix_aux)))) {
    for(i in 1:ncol(matrix_dist)) {
      aux <- matrix_aux[which(rownames(matrix_aux) == rownames(matrix_dist)[i]), , drop = F]
      matrix_dist[i, ] <- apply(aux, 2, min)
      rm(aux)
    }
  } else {
    matrix_dist <- matrix_aux
  }

  rm(list = ls()[ls() != "matrix_dist"])

  return(matrix_dist)

}

#' Polygons' Association Measure
#'
#' @description test statistic for the MC test.
#'
#' @param obj_sp1 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#' @param obj_sp2 a \code{SpatialPolygons} or \code{SpatialPointsDataFrame} object
#'
#' @return a \code{scalar} \code{numeric}
#' @export
#'
pam <- function(obj_sp1, obj_sp2) {
  m_dist <- sp_ID_dist(obj_sp1, obj_sp2)

  min_row <- apply(m_dist, 1, min)
  min_col <- apply(m_dist, 2, min)

  M <- (sum(min_row) + sum(min_col))/(length(min_col) + length(min_row))

  rm(list = ls()[ls() != "M"])

  return(M)

}
