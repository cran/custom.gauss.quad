eval_orthog_poly_mp <- function(x.mp, orthog.poly.mp){
  # This module evaluates the monic orthogonal 
  # polynomial of degree n, where n is a positive
  # integer.
  #
  # Multiple precision is used with nbits.
  #
  # Inputs (both in multiple precision):
  # x.mp: the number at which we want to
  #       evaluate the orthogonal polynomial.
  # orthog.poly.mp: vector of length (n + 1)
  #       consisting of the coefficients of 
  #       the orthogonal polynomial in 
  #       increasing order. 
  #
  # Output:
  # Evaluation of the orthogonal polynomial at x.mp
  # in multiple precision.
  #
  # Written by P.Kabaila in Sept 2021
  
  n <- length(orthog.poly.mp) - 1
  
  powers.mp <- x.mp^c(0:n)
  
  sum(orthog.poly.mp * powers.mp)
  
}
