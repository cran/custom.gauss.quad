#' Moments Computed in Multiple Precision Using the Package Rmpfr
#'
#' @description This module computes the \code{r}'th moment
#'              \deqn{\int_{-\infty}^{\infty} x^r f(x) dx,}
#'  where \eqn{f} is the weight function (specified by the list \code{which.f}), for any nonnegative integer \code{r} using \code{nbits} bits of precision for its computation, via the \code{R} package \code{Rmpfr}.
#'
#'
#' @param which.f a list specifying the nonnegative integrable
#' weight function \eqn{f}, with the following 3 components:
#'  (i)  name (in the form of a character string),
#'  (ii) support specified by a 2-vector of the endpoints of the interval,
#'  (iii)  parameter vector when \eqn{f} belongs to a family of
#'        weight functions and is specified by the value of this
#'        parameter vector (if \eqn{f} is already fully specified
#'        then the parameter vector is set to \code{NULL})
#'
#' @param r nonnegative integer, specifying that it is the \code{r}'th moment for the weight function \eqn{f} that is to be computed
#' @param nbits number of bits in the multiple precision numbers used by the \code{R} package \code{Rmpfr}
#'              to carry out the computation of the \code{r}'th moment
#'
#' @return The \code{r}'th moment with number of bits of precision
#'         \code{nbits} used in its computation, via the \code{R} package \code{Rmpfr}
#'
#' @details
#' Suppose, for example, that we wish to find the Gauss quadrature nodes and weights
#' for the weight function \eqn{f} that is the probability density function of a random
#' variable with the same distribution as \eqn{R/m^{1/2}} where \eqn{R} has a
#' \eqn{\chi_m} distribution (i.e. \eqn{R^2} has a \eqn{\chi_m^2} distribution).
#' In this case, the \code{r}'th moment is
#'     \deqn{\int_{-\infty}^{\infty} x^r f(x) dx
#'     = \left(\frac{2}{m} \right)^{r/2}
#'     \frac{\Gamma((r+m)/2)}{\Gamma(m/2)},}
#' which can be computed to an arbitrary number of bits of precision
#' \code{nbits} using the \code{R} package \code{Rmpfr}.
#' In this case, we specify this weight function \eqn{f} by first
#' assigning the value of \code{m} and then using the \code{R} command
#'
#'  \code{which.f <- list(name="scaled.chi.pdf", support=c(0, Inf), parameters=m)}
#'
#' The code within the function \code{moments} used to compute the
#' \eqn{r}'th moment, to an arbitrary number of bits of precision
#' \code{nbits} using the package \code{Rmpfr}, is listed in the
#' Examples section.
#'
#' @export
#'
#' @examples
#' # The code for the function moments must include a section
#' # that computes the r th moment to an arbitrary number of bits
#' # of precision nbits using the R package Rmpfr for the particular
#' # weight function f of interest.
#' # Suppose that the weight function f is the probability density
#' # function of a random variable with the same probability
#' # distribution as R divided by the square root of m, where R has a
#' # chi distribution with m degrees of freedom.
#' # The code for the function moments includes the following:
#' #
#' #    if (which.f$name == "scaled.chi.pdf"){
#' #    m <- which.f$parameters
#' #    if (r == 0){
#' #    return(mpfr(1, nbits))
#' #    }
#' #    mp.2 <- mpfr(2, nbits)
#' #    mp.r <- mpfr(r, nbits)
#' #    mp.m <- mpfr(m, nbits)
#' #    term1 <- (mp.r/ mp.2) * log(mp.2 / mp.m)
#' #    term2 <- lgamma((mp.r + mp.m) / mp.2)
#' #    term3 <- lgamma(mp.m / mp.2)
#' #    return(exp(term1 + term2 - term3))
#' #    }
#'
#' @seealso
#' \code{custom}
#'
#' @import Rmpfr
#'
moments <- function(which.f, r, nbits){
  # This module computes the r'th moment
  #
  #        infinity
  #   integral       x^r f(x) dx,
  #       -infinity
  #
  # for any nonnegative integer r.
  #
  # The nonnegative integrable function f

  #
  # Multiple precision is used with nbits.
  #
  # Inputs (which are not in multiple precision):
  # which.f: a list specifying f
  # r: nonnegative integer
  # nbits: number of bits in the multiple precision
  #        numbers
  #
  # Output:
  # The r'th moment in multiple precision.
  #
  # Written by P.Kabaila in Sept 2022

  if (which.f$name == "scaled.chi.pdf"){

    m <- which.f$parameters

    if (r == 0){
      return(mpfr(1, nbits))
    }

    mp.2 <- mpfr(2, nbits)
    mp.r <- mpfr(r, nbits)
    mp.m <- mpfr(m, nbits)

    term1 <- (mp.r/ mp.2) * log(mp.2 / mp.m)

    term2 <- lgamma((mp.r + mp.m) / mp.2)

    term3 <- lgamma(mp.m / mp.2)

    return(exp(term1 + term2 - term3))

  }

  if (which.f$name == "Hermite"){

    if (r == 0){
      pi.mp <- Const("pi", nbits)
      return(sqrt(pi.mp))
    }

    if (2 * as.integer(r/2) != r){
      return(mpfr(0, nbits))
    }

    num.mp <- mpfr(r + 1, nbits)
    denom.mp <- mpfr(2, nbits)
    return(gamma(num.mp / denom.mp))

  }

  if (which.f$name == "Generalized.Laguerre"){
    alpha.GGL <- which.f$parameters
    term.mp <- mpfr(r + alpha.GGL + 1, nbits)
    return(gamma(term.mp))

  }

  if (which.f$name == "Legendre"){
    if (2 * as.integer(r/2) != r){
      return(mpfr(0, nbits))
    }
    num.mp <- mpfr(2, nbits)
    denom.mp <- mpfr(r + 1, nbits)
    return(num.mp / denom.mp)

  }

  if (which.f$name == "chemistry.example"){
    mp.1 <- mpfr(1, nbits)
    mp.2 <- mpfr(2, nbits)
    mp.3 <- mpfr(3, nbits)
    mp.r <- mpfr(r, nbits)
    term1 <- mp.3^((mp.r - mp.2) / mp.3)
    term2 <- gamma((mp.r + mp.1) / mp.3)
    return(term1 * term2)

  }

  stop("Unknown which.f in moments.mp.R")

}
