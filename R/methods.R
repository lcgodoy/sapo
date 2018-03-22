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

#' Class of the output from \code{\link{psat_mc}} when the test statistic is \code{\link{pf12}}
#'
#' Internal use
#'
#' @param x Output from \code{\link{psat_mc}}
#'
#' @aliases psa_pf12
#'
#' @export
#'
psa_pf12 <- function(x) {
  x <- append(class(x), "psa_pf12")
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

  if('psa_psam' %in% class(x)){

    if(FALSE) invisible(a <- print(x, ...))

    cat(format('Monte Carlo Polygons Spatial Association Test', ...), '\n', '\n')
    cat(paste('Confidence level:', paste(conf, "%", sep = ""), "\n", sep = " "))
    cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
    cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
    cat(paste('p-value:', round(x$p_value, 8), "\n", sep = " "))
  }

  if('psa_pf12' %in% class(x)){
    cat(format('Monte Carlo Polygons Spatial Association Test', ...), '\n', '\n')
    cat(paste('Confidence level:', paste(conf, "%", sep = ""), "\n", sep = " "))
    cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
    cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
    cat(paste('Decision:', "\n", sep = " "))
    cat(paste('\t', 'F12:', ifelse(x$rejects[1], 'Reject H0', 'Do not reject H0'), '\n', sep = " "))
    cat(paste('\t', 'F21:', ifelse(x$rejects[2], 'Reject H0', 'Do not reject H0'), sep = " "))
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

# plot.psa_test <- function(x) {
#   if(!class(x) %in% "psa_test") stop('Invalid object')
#
#
# }

#' Class of the output from \code{\link{psat_mc}}
#'
#' Internal use
#'
#' @param x Output from \code{\link{pf12}}
#'
#' @aliases pF
#'
#' @export
#'
pF <- function(x) {
  x <- append(class(x), "pF")
}

#' Plot method for \code{\link{pF}}
#'
#' @param x a \code{\link{pF}} object
#' @param ... inherits from \code{plot}.
#'
#' @aliases plot plot.pF
#' @method plot pF
#'
#' @importFrom graphics legend lines
#'
#' @export
plot.pF <- function(x, ...) {
  if(!"pF" %in% class(x)) stop('Invalid object')

  plot(x[,1], x[,2],
       type = 'l', xlab = 'd',
       ylab = '', lty = 2,
       col = 'red',
       ylim = c(min(x[,2], x[,3]),
                max(x[,2], x[,3])),
       ...)

  lines(x[,1], x[,2], lty = 3,
        col = 'blue',
        ...)

  legend('bottomright', col = c('red', 'blue'),
         lty = c(2, 3),
         legend = c((expression(F[12])), (expression(F[21]))),
         bg = 'n',
         bty = 'n')
}
