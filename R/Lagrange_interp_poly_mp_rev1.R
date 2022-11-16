Lagrange_interp_poly_mp_rev1 <- function(j, roots.orthog.poly.mp, n, nbits){
  # This module computes the Lagrange interpolation 
  # polynomial described on page 21 of
  #
  # Gautschi, W. (2004) Orthogonal Polynomials:
  # Computation and Approximation.
  #
  # In other words, this module computes the polynomial
  #
  #       Product          x - tau_i
  #      i=1,...n       --------------
  #   i not equal to j   tau_j - tau_i
  # 
  # where tau_1, tau_2, ... , tau_n are the n elements of 
  # roots.orthog.poly.mp
  #
  # Multiple precision is used with nbits.
  #
  # Inputs:
  # j: integer in the set {1, 2, ..., n} (NOT multiple precision)
  # roots.orthog.poly.mp: vector of the roots of the orthogonal
  #                       polynomial (multiple precision)
  # n: number of Gauss quadrature nodes
  # nbits: number of bits in the multiple precision
  #        numbers (NOT multiple precision)
  #
  # Output:
  # Lagrange interpolation polynomial (degree n-1)
  # represented as a multiple precision vector
  # of length 100
  #
  # Written by P.Kabaila in Aug 2022
  
  n.test <- length(roots.orthog.poly.mp)
  if (n.test != n){
    stop("length of orthog.poly.mp incorrect in Lagrange_interp_poly_mp_rev1.R")
  }
  
  if (n > 100){
    stop("n is too large for Lagrange_interp_poly_mp_rev1.R")
  }
  
  prod.mp <- mpfr(c(1,rep(0, 99)), nbits)
  
  for (i in c(1:n)){
    if (i != j){
      denom <- roots.orthog.poly.mp[j] - roots.orthog.poly.mp[i]
      prod.mp <- prod.mp / denom
      prod.shifted.mp <- c(mpfr(0, nbits), prod.mp[c(1:99)])
      prod.mult.mp <- - roots.orthog.poly.mp[i] * prod.mp
      prod.mp <- prod.shifted.mp + prod.mult.mp
    }
  }
  
  prod.mp
  
}