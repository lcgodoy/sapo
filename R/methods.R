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
  if(!class(x) %in% "psa_test") stop('Invalid object')

  conf <- (1 - x$alpha)*100 %>%
    round(digits = 4)

  if(FALSE) invisible(a <- print(x, ...))

  alt <- switch(x$alternative,
                "two_sided" = "The sets have some kind of association.",
                "attraction" = "The sets attracts each other.",
                "repulsion" = "The sets repulses each other.")

  cat(format('Monte Carlo Polygons Spatial Association Test', ...), '\n', '\n')
  cat(paste('Confidence level:', paste(conf, "%", sep = ""), "\n", sep = " "))
  cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
  cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
  cat(paste('p-value:', round(x$p_value, 8), "\n", sep = " "))
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
  if(!class(object) %in% "psa_test") stop('Invalid object')

  if(FALSE) a <- summary(object, ...)

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
}

# plot.psa_test <- function(x) {
#   if(!class(x) %in% "psa_test") stop('Invalid object')
#
#
# }
