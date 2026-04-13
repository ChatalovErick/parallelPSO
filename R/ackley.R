#' Ackley Function
#'
#' The Ackley function is a non-convex function used for testing optimization algorithms.
#' It has many local minima, making it challenging for optimization algorithms.
#' The global minimum is at x = (0, 0, ..., 0) with f(x) = 0.
#'
#' @param xx Numeric vector of parameters to evaluate.
#' @param a Constant parameter (default: 20).
#' @param b Constant parameter (default: 0.2).
#' @param c Constant parameter (default: 2*pi).
#'
#' @return The value of the Ackley function at the given point.
#'
#' @details
#' The Ackley function is defined as:
#' f(x) = -a * exp(-b * sqrt(sum(x^2) / d)) - exp(sum(cos(c * x)) / d) + a + exp(1)
#' where d is the dimension of x.
#'
#' @references
#' Surjanovic, S., & Bingham, D. (2013). Virtual Library of Simulation Experiments:
#' Test Functions and Datasets. Simon Fraser University.
#' http://www.sfu.ca/~ssurjano/ackley.html
#'
#' @examples
#' # Evaluate Ackley function at origin (should be 0)
#' ackley(rep(0, 5))
#'
#' # Evaluate at random point
#' ackley(runif(5, -5, 5))
#'
#' @export
ackley <- function(xx, a = 20, b = 0.2, c = 2 * pi) {

  d <- length(xx)

  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))

  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)

  y <- term1 + term2 + a + exp(1)
  return(y)
}
