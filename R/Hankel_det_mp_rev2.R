Hankel_det_mp_rev2 <- function(which.f, n, nbits){
  # This module computes the Hankel determinant
  # of order n, where n is a nonnegative integer,
  # as defined by (2.1.1) on page 53 of
  #
  # Gautschi, W. (2004) Orthogonal Polynomials:
  # Computation and Approximation.
  #
  # The nonnegative integrable function f
  # is specified by the list which.f with
  # the following 3 components:
  #  (i)  name,
  #  (ii) support specified by a 2-vector
  #    of the endpoints of the interval,
  #  (iii) parameter vector
  #
  # Multiple precision is used with nbits.
  #
  # Inputs (which are not in multiple precision):
  # which.f: specifies f
  # n: order of the determinant, a nonnegative integer
  # nbits: number of bits in the multiple precision
  #        numbers
  #
  # Output:
  # Hankel determinant of order n
  # in multiple precision.
  #
  # Written by P.Kabaila in Sept 2022

  if (n == 0){
    return(mpfr(1, nbits))
  }

  Hankel.mat <- mpfr(rep(0, n^2), nbits)
  dim(Hankel.mat) <- c(n, n)

  for (i in c(1:n)){
    for (j in c(1:n)){
      r <- (j - 1) + (i - 1)
      Hankel.mat[i, j] <- moments(which.f, r, nbits)
    }
  }

  if (n <= 2){
    tmp <- determinant(Hankel.mat, logarithm = FALSE, asNumeric = FALSE)
    return(tmp$sign * tmp$modulus)
  }else{
    return(det_LU_Dolittle(Hankel.mat, nbits))
  }

}
