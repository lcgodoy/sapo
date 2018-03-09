#' Class of the output from \code{\link{pat_mc}}
#'
#' Internal use
#'
#' @param x Output from \code{\link{pat_mc}}
#'
#' @aliases pat_test
#'
#' @export
#'
pat_test <- function(x) {
  x <- append(class(x), "pat_test")
}


#' Print method for \code{\link{pat_test}}
#'
#' @param x A \code{\link{pat_test}} object
#' @param ... inherits from \code{print} and \code{format}
#'
#' @aliases print print.pat_test
#' @method print pat_test
#' @export
#'
print.pat_test <- function(x, ...) {
  if(!class(x) %in% "pat_test") stop('Invalid object')

  conf <- (1 - x$alpha)*100 %>%
    round(digits = 4)

  if(FALSE) invisible(a <- print(x, ...))

  alt <- switch(x$alternative,
                "two_sided" = "The sets have some kind of association.",
                "attraction" = "The sets attracts each other.",
                "repulsion" = "The sets repulses each other.")

  cat(format('Monte Carlo Polygons Association Test', ...), '\n', '\n')
  cat(paste('Confidence level:', paste(conf, "%", sep = ""), "\n", sep = " "))
  cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
  cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
  cat(paste('p-value:', round(x$p_value, 8), "\n", sep = " "))
}

#' Summary method for \code{\link{pat_test}}
#'
#' @param object a \code{\link{pat_test}} object
#' @param ... inherits from \code{print}, \code{format}, and \code{summary}
#'
#' @aliases summary summary.pat_test
#' @method summary pat_test
#' @export
summary.pat_test <- function(object, ...) {
  if(!class(object) %in% "pat_test") stop('Invalid object')

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

  cat("Monte Carlo Polygons Association Test", '\n \n')
  cat(paste('Null hypothesis:', "The sets are independent.", "\n", sep = " "))
  cat(paste('Alternative hypothesis:', alt, "\n", sep = " "))
  cat(paste('Decision:', dec,
            "with confidence level of",
            paste(conf, "%", sep = ""),
            "\n", sep = " "))
  cat(paste('p-value:', round(object$p_value, 6), "\n", sep = " "))
}

# plot.pat_test <- function(x) {
#   if(!class(x) %in% "pat_test") stop('Invalid object')
#
#
# }
