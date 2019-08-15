#' \code{\link{mc_psa}} class
#'
#' Internal use
#'
#' @param x Output from \code{\link{psam_mc}} or
#' \code{\link{gof_mc}}
#'
#'
mc_psa <- function(x) {
  x <- append(class(x), "mc_psa")
}

#' \code{\link{psam_mc}} class
#'
#' Internal use
#'
#' @param x Output from \code{\link{psam_mc}}
#'
#'
psam_test <- function(x) {
  x <- append(class(x), "psam_test")
}

#' \code{\link{gof_mc}} class
#'
#' Internal use
#'
#' @param x Output from \code{\link{gof_mc}}
#'
gof_test <- function(x) {
  x <- append(class(x), "gof_test")
}

#' Print method for \code{\link{mc_psa}}
#'
#' @param x A \code{\link{mc_psa}} object
#' @param ... inherits from \code{print}.
#'
#' @method print mc_psa
#' @export
#'
print.mc_psa <- function(x, ...) {
  conf <- round((1 - x$alpha)*100, digits = 4)

  if(F) print(x$mc_sample, ...)

  cat('Monte Carlo Polygons Spatial Association Test', '\n', '\n')
  cat('Null hypothesis:', "The polygons' sets are independent.", "\n")
  cat('p-value:', round(x$p_value, 8), "\n")
}

#' Summary method for \code{\link{mc_psa}}
#'
#' @param object a \code{\link{mc_psa}} object
#' @param ... inherits from \code{summary}.
#'
#' @method summary mc_psa
#' @export
summary.mc_psa <- function(object, ...) {
  cat("Monte Carlo Polygons Spatial Association Test", '\n \n')
  cat('Test Statistic Summary \n')
  cat('\t')
  summary(object$mc_sample, ...)
  cat('\n')
  cat('p-value:', round(object$p_value, 6), "\n")
}

#' Plot method for \code{\link{mc_psa}}
#'
#' @param x a \code{\link{psa_test}} object
#' @param ... inherits from \code{plot}
#'
#' @import graphics
#' @importFrom stats quantile
#' @aliases plot plot.psa_test
#' @method plot psa_test
#' @export
#'
# plot.psa_test <- function(x, ...) {
#   if('psam_test' %in% class(x)) {
#     if(x$alternative == 'two_sided') {
#       if (requireNamespace("ggplot2", quietly = TRUE)) {
#         ggplot2::ggplot() +
#           ggplot2::geom_line(ggplot2::aes(x = x$mc_sample), stat = 'density') +
#           ggplot2::geom_area(ggplot2::aes(x = x$mc_sample,
#                                           fill = {x$mc_sample >= quantile(x$mc_sample, 1 - (x$alpha/2))}
#           ), stat = 'density'
#           ) +
#           ggplot2::geom_area(ggplot2::aes(x$mc_sample,
#                                           fill = {x$mc_sample <= quantile(x$mc_sample, x$alpha/2)}
#           ), stat = 'density'
#           ) +
#           ggplot2::scale_fill_manual(values = c('FALSE' = 'transparent', 'TRUE' = '#ff000080')) +
#           ggplot2::geom_vline(xintercept = x$sample_ts,
#                               linetype = 'dashed', col = 'red') +
#           theme_tpsa() +
#           ggplot2::labs(x = 'Test Statistic', y = 'Kernel density') +
#           ggplot2::ggtitle(label = 'Monte Carlo Test',
#                            subtitle = 'Polygons Spatial Association Measure') +
#           ggplot2::guides(fill = F)
#       } else {
#         stop('install ggplot2 to visualize the plot.')
#       }
#     } else {
#       if(x$alternative == 'repulsion') {
#         if (requireNamespace("ggplot2", quietly = TRUE)) {
#           ggplot2::ggplot() +
#             ggplot2::geom_line(ggplot2::aes(x = x$mc_sample), stat = 'density') +
#             ggplot2::geom_area(ggplot2::aes(x = x$mc_sample,
#                                             fill = {x$mc_sample >= quantile(x$mc_sample, 1 - x$alpha)}
#             ), stat = 'density'
#             ) +
#             ggplot2::scale_fill_manual(values = c('FALSE' = 'transparent', 'TRUE' = '#ff000080')) +
#             ggplot2::geom_vline(xintercept = x$sample_ts,
#                                 linetype = 'dashed', col = 'red') +
#             theme_tpsa() +
#             ggplot2::labs(x = 'Test Statistic', y = 'Kernel density') +
#             ggplot2::ggtitle(label = 'Monte Carlo Test',
#                              subtitle = 'Polygons Spatial Association Measure') +
#             ggplot2::guides(fill = F)
#         } else {
#           stop('install ggplot2 to visualize the plot.')
#         }
#       } else {
#         if(x$alternative == 'attraction') {
#           if (requireNamespace("ggplot2", quietly = TRUE)) {
#             ggplot2::ggplot() +
#               ggplot2::geom_line(ggplot2::aes(x = x$mc_sample), stat = 'density') +
#               ggplot2::geom_area(ggplot2::aes(x = x$mc_sample,
#                                               fill = {x$mc_sample <= quantile(x$mc_sample, x$alpha)}
#               ), stat = 'density'
#               ) +
#               ggplot2::scale_fill_manual(values = c('FALSE' = 'transparent', 'TRUE' = '#ff000080')) +
#               ggplot2::geom_vline(xintercept = x$sample_ts,
#                                   linetype = 'dashed', col = 'red') +
#               theme_tpsa() +
#               ggplot2::labs(x = 'Test Statistic', y = 'Kernel density') +
#               ggplot2::ggtitle(label = 'Monte Carlo Test',
#                                subtitle = 'Polygons Spatial Association Measure') +
#               ggplot2::guides(fill = F)
#
#           } else {
#             stop('install ggplot2 to visualize the plot.')
#           }
#         }
#       }
#     }
#
#   } else {
#     if (requireNamespace("ggplot2", quietly = TRUE)) {
#       ggplot2::ggplot() +
#         ggplot2::geom_line(ggplot2::aes(x = x$distances, y = x$mc_sample)) +
#         ggplot2::geom_line(ggplot2::aes(x = r, y = k_inf),
#                            color = 'red', linetype = 2, inherit.aes = F) +
#         ggplot2::geom_line(ggplot2::aes(x = r, y = k_up),
#                            color = 'red', linetype = 2, inherit.aes = F) +
#         theme_tpsa() +
#         ggplot2::labs(x = 'Distance', y = expression(K[12](d))) +
#         ggplot2::ggtitle(label = 'Monte Carlo Test')
#     } else {
#       plot(x$sample_ts[, 1], x$sample_ts[, 2],
#            type = 'l',
#            xlab = 'Distance',
#            ylab = '',
#            main = expression(K['12'](d)),
#            bty = 'l', ...)
#       grid()
#       lines(x$mc_sample$r, x$mc_sample$k12_up, lty = 2, col = 'red')
#       lines(x$mc_sample$r, x$mc_sample$k12_inf, lty = 2, col = 'red')
#     }
#   }
# }
