#' Class of the output from \code{\link{psat_mc}}
#'
#' Internal use
#'
#' @param x Output from \code{\link{psat_mc}}
#'
#' @aliases psa_test
#'
#' @export
#'
psa_test <- function(x) {
  x <- append(class(x), "psa_test")
}

#' Class of the output from \code{\link{psat_mc}} when the test statistic is \code{\link{psam}}
#'
#' Internal use
#'
#' @param x Output from \code{\link{psat_mc}}
#'
#' @aliases psa_psam
#'
#' @export
#'
psa_psam <- function(x) {
  x <- append(class(x), "psa_psam")
}


#' Class of the output from \code{\link{psat_mc}} when the test statistic is \code{\link{pk_dist12}}
#' or \code{\link{pk_area12}}
#'
#' Internal use
#'
#' @param x Output from \code{\link{psat_mc}}
#'
#' @aliases psa_pk12
#'
#' @export
#'
psa_pk12 <- function(x) {
  x <- append(class(x), "psa_pk12")
}

#' Print method for \code{\link{psa_test}}
#'
#' @param x A \code{\link{psa_test}} object
#' @param ... inherits from \code{print} and \code{format}
#'
#' @aliases print print.psa_test
#' @method print psa_test
#' @export
#'
print.psa_test <- function(x, ...) {
  if(!"psa_test" %in% class(x)) stop('Invalid object')

  conf <- (1 - x$alpha)*100 %>%
    round(digits = 4)

  alt <- switch(x$alternative,
                "two_sided" = "The sets have some kind of association.",
                "attraction" = "The sets attracts each other.",
                "repulsion" = "The sets repulses each other.")

  if('psa_psam' %in% class(x)) {

    if(FALSE) invisible(a <- print(x, ...))

    cat(format('Monte Carlo Polygons Spatial Association Test', ...), '\n', '\n')
    cat(paste('Confidence level:', paste(conf, "%", sep = ""), "\n", sep = " "))
    cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
    cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
    cat(paste('p-value:', round(x$p_value, 8), "\n", sep = " "))
  }

  if('psa_pk12' %in% class(x)) {
    cat(format('Monte Carlo Polygons Spatial Association Test', ...), '\n', '\n')
    cat(paste('Confidence level:', paste(conf, "%", sep = ""), "\n", sep = " "))
    cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
    cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
    cat(paste('Decision:', "\n", sep = " "))
    cat(paste('\t', 'K12:', ifelse(x$rejects, 'Reject H0', 'Do not reject H0'), sep = " "))
  }

}

#' Summary method for \code{\link{psa_test}}
#'
#' @param object a \code{\link{psa_test}} object
#' @param ... inherits from \code{print}, \code{format}, and \code{summary}
#'
#' @aliases summary summary.psa_test
#' @method summary psa_test
#' @export
summary.psa_test <- function(object, ...) {
  if(!"psa_test" %in% class(object)) stop('Invalid object')

  if(FALSE) a <- summary(object, ...)

  if('psa_psam' %in% class(object)) {
    alt <- switch(object$alternative,
                  "two_sided" = "The sets have some kind of association.",
                  "attraction" = "The sets attracts each other.",
                  "repulsion" = "The sets repulses each other.")

    conf <- (1 - object$alpha)*100 %>%
      round(digits = 4)

    if(object$alternative == "two_sided") {
      dec <- ifelse(object$p_value < (object$alpha/2),
                    "Reject the null hypothesis",
                    "Do not reject the null hypothesis")
    } else {
      dec <- ifelse(object$p_value < object$alpha,
                    "Reject the null hypothesis",
                    "Do not reject the null hypothesis")
    }

    cat("Monte Carlo Polygons Spatial Association Test", '\n \n')
    cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
    cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
    cat(paste('Decision:', dec,
              "with confidence level of",
              paste(conf, "%", sep = ""),
              "\n", sep = " "))
    cat(paste('p-value:', round(object$p_value, 6), "\n", sep = " "))
  } else {
    cat('Summary not implemented for this test statistic yet.')
  }
}

#' PLot method for \code{\link{psa_test}}
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
plot.psa_test <- function(x, ...) {
  if(!"psa_test" %in% class(x)) stop('Invalid object')

  if('psa_psam' %in% class(x)) {
    if(x$alternative == 'two_sided') {
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(x = x$mc_ts), stat = 'density') +
          ggplot2::geom_area(ggplot2::aes(x = x$mc_ts,
                                 fill = {x$mc_ts >= quantile(x$mc_ts, 1 - (x$alpha/2))}
          ), stat = 'density'
          ) +
          ggplot2::geom_area(ggplot2::aes(x$mc_ts,
                                 fill = {x$mc_ts <= quantile(x$mc_ts, x$alpha/2)}
          ), stat = 'density'
          ) +
          ggplot2::scale_fill_manual(values = c('FALSE' = 'transparent', 'TRUE' = '#ff000080')) +
          ggplot2::geom_vline(xintercept = x$sample_ts,
                              linetype = 'dashed', col = 'red') +
          theme_tpsa() +
          ggplot2::labs(x = 'Test Statistic', y = 'Kernel density') +
          ggplot2::ggtitle(label = 'Monte Carlo Test',
                           subtitle = 'Polygons Spatial Association Measure') +
          ggplot2::guides(fill = F)
      } else {
        stop('install ggplot2 to visualize the plot.')
      }
    } else {
      if(x$alternative == 'repulsion') {
        if (requireNamespace("ggplot2", quietly = TRUE)) {
          ggplot2::ggplot() +
            ggplot2::geom_line(ggplot2::aes(x = x$mc_ts), stat = 'density') +
            ggplot2::geom_area(ggplot2::aes(x = x$mc_ts,
                                            fill = {x$mc_ts >= quantile(x$mc_ts, 1 - x$alpha)}
            ), stat = 'density'
            ) +
            ggplot2::scale_fill_manual(values = c('FALSE' = 'transparent', 'TRUE' = '#ff000080')) +
            ggplot2::geom_vline(xintercept = x$sample_ts,
                                linetype = 'dashed', col = 'red') +
            theme_tpsa() +
            ggplot2::labs(x = 'Test Statistic', y = 'Kernel density') +
            ggplot2::ggtitle(label = 'Monte Carlo Test',
                             subtitle = 'Polygons Spatial Association Measure') +
            ggplot2::guides(fill = F)
        } else {
          stop('install ggplot2 to visualize the plot.')
        }
      } else {
        if(x$alternative == 'attraction') {
          if (requireNamespace("ggplot2", quietly = TRUE)) {
            ggplot2::ggplot() +
              ggplot2::geom_line(ggplot2::aes(x = x$mc_ts), stat = 'density') +
              ggplot2::geom_area(ggplot2::aes(x = x$mc_ts,
                                              fill = {x$mc_ts <= quantile(x$mc_ts, x$alpha)}
              ), stat = 'density'
              ) +
              ggplot2::scale_fill_manual(values = c('FALSE' = 'transparent', 'TRUE' = '#ff000080')) +
              ggplot2::geom_vline(xintercept = x$sample_ts,
                                  linetype = 'dashed', col = 'red') +
              theme_tpsa() +
              ggplot2::labs(x = 'Test Statistic', y = 'Kernel density') +
              ggplot2::ggtitle(label = 'Monte Carlo Test',
                               subtitle = 'Polygons Spatial Association Measure') +
              ggplo2::guides(fill = F)

          } else {
            stop('install ggplot2 to visualize the plot.')
          }
        }
      }
    }

  } else {
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        df_gg <- x$mc_ts[,1:3]
        names(df_gg) <- c('r', 'k_inf', 'k_up')
        df_gg$obs <- x$sample_ts[,2]
        df_gg$func <- 'K[12](d)'

        ggplot2::ggplot(data = df_gg) +
          ggplot2::geom_line(ggplot2::aes(x = r, y = obs)) +
          ggplot2::geom_line(ggplot2::aes(x = r, y = k_inf),
                             color = 'red', linetype = 2, inherit.aes = F) +
          ggplot2::geom_line(ggplot2::aes(x = r, y = k_up),
                             color = 'red', linetype = 2, inherit.aes = F) +
          theme_tpsa() +
          ggplot2::labs(x = 'Distance', y = expression(K[12](d))) +
          ggplot2::ggtitle(label = 'Monte Carlo Test')
      } else {
        plot(x$sample_ts[, 1], x$sample_ts[, 2],
             type = 'l',
             xlab = 'Distance',
             ylab = '',
             main = expression(K['12'](d)),
             bty = 'l', ...)
        grid()
        lines(x$mc_ts$r, x$mc_ts$k12_up, lty = 2, col = 'red')
        lines(x$mc_ts$r, x$mc_ts$k12_inf, lty = 2, col = 'red')
      }
    }
}
