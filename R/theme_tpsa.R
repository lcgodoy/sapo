#' \code{ggplo2} theme for \code{tpsa}
#'
#' @param ... Named argument to modify the theme
#'
#' @export
theme_tpsa <- function(...) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = 'black'),
                     strip.text = ggplot2::element_text(color = 'white'),
                     panel.grid = ggplot2::element_line(colour = '#BEBEBE66',
                                               linetype = 'longdash'),
                     axis.ticks.length = ggplot2::unit(.05, units = 'cm'), ...)
  } else {
    stop('\n Install ggplot2 in order to use it. \n')
  }
}
