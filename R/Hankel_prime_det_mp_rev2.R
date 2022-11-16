Hankel_prime_det_mp_rev2 <- function(which.f, n, nbits){
  # This module computes the Hankel determinant
  # of order n, where n is a nonnegative integer,
  # as defined by (2.1.2) on page 53 of
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
  # Hankel prime determinant of order n
  # in multiple precision.
  #
  # Written by P.Kabaila in Sept 2022

  if (n == 0){
    return(mpfr(0, nbits))
  }

  if (n == 1){
    tmp <- moments(which.f, 1, nbits)
    return(tmp)
  }

  Hankel.prime.part1.mat <- mpfr(rep(0, n*(n-1)), nbits)
  dim(Hankel.prime.part1.mat) <- c(n, (n-1))

  for (i in c(1:n)){
    for (j in c(1:(n-1))){
      r <- (j - 1) + (i - 1)
      Hankel.prime.part1.mat[i, j] <- moments(which.f, r, nbits)
    }
  }

  Hankel.prime.last.col <- mpfr(rep(0, n), nbits)
  dim(Hankel.prime.last.col) <- c(n, 1)
  for (i in c(1:n)){
    r <- i + n - 1
    Hankel.prime.last.col[i] <- moments(which.f, r, nbits)
  }

  Hankel.prime.mat <- cbind(Hankel.prime.part1.mat,
                            Hankel.prime.last.col)

  colnames(Hankel.prime.mat) <- NULL

  if (n <= 2){
    tmp <- determinant(Hankel.prime.mat, logarithm = FALSE,
                       asNumeric = FALSE)
    return(tmp$sign * tmp$modulus)
  }else{
    return(det_LU_Dolittle(Hankel.prime.mat, nbits))
  }

}
