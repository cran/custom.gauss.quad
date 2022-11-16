#' Custom-Made Gauss Quadrature Nodes and Weights
#'
#' For the nonnegative weight function specified by \code{which.f} and given
#' number \code{n} of nodes, the function \code{custom} computes the
#' Gauss quadrature nodes \code{asNumeric.nodes} and corresponding weights
#' \code{asNumeric.weights} which are double precision vectors.
#'
#' @param which.f a list specifying the nonnegative integrable
#' weight function \eqn{f}, with the following three components:
#'  (i)  name (in the form of a character string),
#'  (ii) support specified by a 2-vector of the endpoints of the interval,
#'  (iii) parameter vector when \eqn{f} belongs to a family of
#'        weight functions and is specified by the value of this
#'        parameter vector (if \eqn{f} is already fully specified
#'        then the parameter vector is set to \code{NULL})
#'
#' @param n number of Gauss quadrature nodes
#'
#' @return
#' A list with the following elements:
#' \code{which.f},
#' \code{n},
#' \code{nbits.vec},
#' \code{list.Gauss.nodes},
#' \code{list.Gauss.weights},
#' \code{mat.timings},
#' \code{max.abs.diffs.nodes},
#' \code{sum.abs.diffs.weights},
#' \code{L.nodes},
#' \code{L.weights},
#' \code{asNumeric.nodes},
#' \code{asNumeric.weights}.
#'
#' @details Suppose that we wish to evaluate
#' \deqn{\int_{-\infty}^{\infty} g(x) f(x) dx,}
#' where \eqn{f} is a specified nonnegative integrable
#' weight function. The Gauss quadrature approximation to this
#' integral has the form
#' \deqn{\sum_{i=1}^n \lambda_i \, g(\tau_i),}
#' where \eqn{\tau_1, \dots, \tau_n} are called the nodes
#' and \eqn{\lambda_1, \dots, \lambda_n} are called the
#' corresponding weights. This approximation is exact
#' whenever \eqn{g} is a polynomial of degree less
#' than or equal to \eqn{2n - 1}.
#'
#' If \eqn{f} takes a form that leads to Gauss quadrature
#' rules with nodes that are the roots of classical orthogonal
#' polynomials of a continuous variable then these rules
#' (such as Gauss Legendre, Gauss Hermite and Gauss Laguerre) are
#' readily accessible to statisticians via \code{R} packages
#' such as \code{statmod}.
#' If, however, \eqn{f} does not take one of these
#' particular forms then the Gauss quadrature nodes and weights
#' need to be custom-made.
#'
#' \code{custom} computes the Gauss quadrature nodes
#' and weights, for given \code{n}, using a user-supplied formula
#' for the \eqn{r}'th moment
#' \deqn{\int_{-\infty}^{\infty} x^r f(x) dx}
#' for all nonnegative integers \eqn{r}. This formula must be inserted
#' by the user into the code for the function
#' \code{moments} and must be able to be computed
#' to an arbitrary number of bits (\code{nbits}) of precision
#' using the \code{R} package \code{Rmpfr}.
#'
#' To obtain some assurance that the Gauss quadrature nodes
#' and weights are computed to sufficient precision for subsequent
#' double precision computations in \code{R}, these nodes and
#' weights are computed for the increasing numbers of bits of
#' precision given in the 5-vector \code{nbits.vec} and the results
#' compared. This comparison results in the criteria
#' \code{max.abs.diffs.nodes}, \code{sum.abs.diffs.weights},
#' \code{L.nodes} and \code{L.weights} described in detail by
#' Kabaila (2022). The execution times for various parts of the
#' code are stored in \code{mat.timings} whose components are
#' described by Kabaila (2022).
#'
#' \code{list.Gauss.nodes[[i]]} and \code{list.Gauss.weights[[i]]}
#' are the \code{n}-vectors of Gauss quadrature nodes and weights,
#' respectively, computed using \code{nbits.vec[i]} bits of precision
#' (\code{i=1,...,5}).
#'
#' The most accurate approximations to the Gauss quadrature nodes and weights are
#' \code{list.Gauss.nodes[[5]]} and \code{list.Gauss.weights[[5]]}. These are
#' converted to double precision by applying the \code{asNumeric}
#' function from the \code{R} package \code{Rmpfr}, resulting
#' in \code{asNumeric.nodes} and \code{asNumeric.weights}, respectively.
#'
#' @references
#' Kabaila, P. (2022) Custom-made Gauss quadrature for statisticians. arXiv:2211.04729
#'
#' @seealso
#' \code{moments}
#'
#' @export
#'
#' @examples
#' # Suppose that the weight function f is the probability density
#' # function of a random variable with the same probability
#' # distribution as R divided by the square root of m, where R has a
#' # chi distribution with m degrees of freedom.
#' # Also suppose that we wish to compute the Gauss quadrature nodes
#' # and weights, for number of nodes n = 5, when the parameter m = 160.
#' # The r th moment can be computed to an arbitrary number of bits of
#' # precision using the R package Rmpfr. We describe the weight function
#' # f using the following R commands:
#'
#' m <- 160
#' which.f <- list(name="scaled.chi.pdf", support=c(0, Inf),
#' parameters=m)
#'
#' # Here, "scaled.chi.pdf" is the name (a character string) that we
#' # have given to the weight function f. The R function moments includes
#' # the code needed to compute the r th moment to an arbitrary number
#' # of bits of precision using the R package Rmpfr.
#' # We compute the Gauss quadrature node and weight, for the toy example
#' # with number of nodes n=1, using the following R commands:
#'
#' n <- 1
#' gauss.list <- custom(which.f, n)
#' old <- options(digits = 17)
#' gauss.list$asNumeric.nodes
#' gauss.list$asNumeric.weights
#' options(old)
#'
#' # These commands take less than 1 second to run. The resulting
#' # of node and corresponding weight in double precision are:
#' # > gauss.list$asNumeric.nodes
#' # [1] 0.99843873022375829
#' # > gauss.list$asNumeric.weights
#' # [1] 1
#'
#' # The computation times for number of nodes n=5, 17 and 33 are roughly
#' # 160 seconds, 31 minutes and 5 hours,respectively.
#' #
#' # We compute the Gauss quadrature nodes and weights, for number of
#' # nodes n=5, using the following R commands:
#'\donttest{
#' n <- 5
#' gauss.list <- custom(which.f, n)
#' old <- options(digits = 17)
#' gauss.list$asNumeric.nodes
#' gauss.list$asNumeric.weights
#' options(old)
#' }
#' # These commands take roughly 3 minutes to run. The resulting vectors
#' # of nodes and corresponding weights in double precision are:
#' # > gauss.list$asNumeric.nodes
#' # [1] 0.84746499810651410 0.92785998378868118 1.00262691212158761
#' # [4] 1.07930375924992528 1.16628363226782716
#' # > gauss.list$asNumeric.weights
#' # [1] 0.0144433732487188448 0.2483585328946608384 0.5305446123744097520
#' # [4] 0.1977278905956056654 0.0089255908866048821
#'
#'
#' @import Rmpfr
#'
custom <- function(which.f, n){
    # This module implements the moment-based method via
    # moment determinants, where multiple precision
    # evaluation of moments through a formula is used
    # throughout, for the values of nbits specified in the
    # vector nbits.vec.
  nbits.first <- ceiling(60 + 6.5 * n)
  nbits.vec <- seq(nbits.first, nbits.first + 136, length.out = 5)
  len.nbits.vec <- length(nbits.vec)

  # cat(
  #   "nbits.vec = (",
  #   nbits.vec[1],
  #   ",",
  #   nbits.vec[2],
  #   ",",
  #   nbits.vec[3],
  #   ",",
  #   nbits.vec[4],
  #   ",",
  #   nbits.vec[5],
  #   ")",
  #   "\n"
  # )

  mat.timings <- matrix(0, nrow = 3, ncol = len.nbits.vec)

  if (n == 1) {
    list.Gauss.nodes <- NULL
    list.Gauss.weights <- NULL
    for (j in c(1:len.nbits.vec)) {
      nbits <- nbits.vec[j]
      mu0 <- moments(which.f, 0, nbits)
      mu1 <- moments(which.f, 1, nbits)
      Gauss.nodes <- mu1 / mu0
      list.Gauss.nodes <- append(list.Gauss.nodes, Gauss.nodes)
      Gauss.weights <- mu0
      list.Gauss.weights <-
        append(list.Gauss.weights, Gauss.weights)
    }
  } else{
    nbits <- nbits.vec[1]
    # times.taken[1] is the "user time"
    times.taken <- system.time(orthog.poly <-
                                 orthog_poly_mp_rev2(which.f, n, nbits))
    mat.timings[1, 1] <- times.taken[1]

    times.taken <- system.time(Gauss.nodes <-
          roots_orthog_poly_robust_mp_rev3(which.f, orthog.poly, n, nbits))
    list.Gauss.nodes <- list(Gauss.nodes)
    mat.timings[2, 1] <- times.taken[1]

    times.taken <- system.time(Gauss.weights <-
              weights_from_moments_mp_rev2(which.f, Gauss.nodes, n, nbits))
    list.Gauss.weights <- list(Gauss.weights)
    mat.timings[3, 1] <- times.taken[1]

    # cat("Computations for element number",
    #     1,
    #     "of nbits.vec completed",
    #     "\n")

    for (j in c(2:len.nbits.vec)) {
      nbits <- nbits.vec[j]

      times.taken <- system.time(orthog.poly <-
                                   orthog_poly_mp_rev2(which.f, n, nbits))
      mat.timings[1, j] <- times.taken[1]

      times.taken <- system.time(Gauss.nodes <-
          roots_orthog_poly_robust_mp_rev3(which.f, orthog.poly, n, nbits))
      list.Gauss.nodes <-
        append(list.Gauss.nodes, list(Gauss.nodes))
      mat.timings[2, j] <- times.taken[1]

      times.taken <- system.time(Gauss.weights <-
              weights_from_moments_mp_rev2(which.f, Gauss.nodes, n, nbits))
      list.Gauss.weights <-
        append(list.Gauss.weights, list(Gauss.weights))
      mat.timings[3, j] <- times.taken[1]

      # cat("Computations for element number",
      #     j,
      #     "of nbits.vec completed",
      #     "\n")

    }

  }

  max.abs.diffs.nodes <- rep(0, (len.nbits.vec - 1))
  for (j in c(2:len.nbits.vec)) {
    diff.nodes <- list.Gauss.nodes[[j]] - list.Gauss.nodes[[(j - 1)]]
    max.abs.diffs.nodes[(j - 1)] <-
      max(abs(asNumeric(diff.nodes)))
  }

  sum.abs.diffs.weights <- rep(0, (len.nbits.vec - 1))
  for (j in c(2:len.nbits.vec)) {
    diff.weights <-
      list.Gauss.weights[[j]] - list.Gauss.weights[[(j - 1)]]
    sum.abs.diffs.weights[(j - 1)] <-
      sum(abs(asNumeric(diff.weights)))
  }

  # Compute L.nodes
  mat.nodes.dble.prec <- matrix(0, nrow = n, ncol = len.nbits.vec)
  for (j in c(1:len.nbits.vec)) {
    mat.nodes.dble.prec[, j] <-
      asNumeric(list.Gauss.nodes[[j]])
  }

  diff.mat.nodes.dble.prec <-
    matrix(0, nrow = n, ncol = (len.nbits.vec - 1))
  for (j in c(2:len.nbits.vec)) {
    diff.mat.nodes.dble.prec[, (j - 1)] <-
      mat.nodes.dble.prec[, j] -
      mat.nodes.dble.prec[, (j - 1)]
  }

  temp.nodes.vec <- colSums(abs(diff.mat.nodes.dble.prec))

  if (temp.nodes.vec[(len.nbits.vec - 1)] != 0) {
    L.nodes <- NULL
  } else{
    j <- len.nbits.vec - 1
    while (j > 0 && temp.nodes.vec[j] == 0) {
      j <- j - 1
    }
    L.nodes <- j + 1
  }

  # Compute L.weights
  mat.weights.dble.prec <- matrix(0, nrow = n, ncol = len.nbits.vec)
  for (j in c(1:len.nbits.vec)) {
    mat.weights.dble.prec[, j] <-
      asNumeric(list.Gauss.weights[[j]])
  }

  diff.mat.weights.dble.prec <-
    matrix(0, nrow = n, ncol = (len.nbits.vec - 1))
  for (j in c(2:len.nbits.vec)) {
    diff.mat.weights.dble.prec[, (j - 1)] <-
      mat.weights.dble.prec[, j] -
      mat.weights.dble.prec[, (j - 1)]
  }

  temp.weights.vec <- colSums(abs(diff.mat.weights.dble.prec))

  if (temp.weights.vec[(len.nbits.vec - 1)] != 0) {
    L.weights <- NULL
  } else{
    j <- len.nbits.vec - 1
    while (j > 0 && temp.weights.vec[j] == 0) {
      j <- j - 1
    }
    L.weights <- j + 1
  }

  asNumeric.nodes <- mat.nodes.dble.prec[, len.nbits.vec]
  asNumeric.weights <- mat.weights.dble.prec[, len.nbits.vec]

  list(
    which.f = which.f,
    n = n,
    nbits.vec = nbits.vec,
    list.Gauss.nodes = list.Gauss.nodes,
    list.Gauss.weights = list.Gauss.weights,
    mat.timings = mat.timings,
    max.abs.diffs.nodes = max.abs.diffs.nodes,
    sum.abs.diffs.weights = sum.abs.diffs.weights,
    L.nodes = L.nodes,
    L.weights = L.weights,
    asNumeric.nodes = asNumeric.nodes,
    asNumeric.weights = asNumeric.weights
  )

}
