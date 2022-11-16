sq_Lagrange_interp_poly_mp_rev1 <- 
  function(Lagrange.interp.poly.mp, n, nbits){
  # This module computes the square of the Lagrange interpolation 
  # polynomial described on page 21 of
  #
  # Gautschi, W. (2004) Orthogonal Polynomials:
  # Computation and Approximation.
  #
  # Multiple precision is used with nbits.
  #
  # Inputs:
  # Lagrange.interp.poly.mp: 
  #          Lagrange interpolation polynomial (degree n-1)
  #          represented as a vector of length 100 of 
  #          multiple precision coefficients
  # n: number of Gauss quadrature nodes
  # nbits: number of bits in the multiple precision
  #        numbers (NOT multiple precision)
  #
  # Output:
  # Square of the Lagrange interpolation polynomial
  # (degree 2(n-1)) represented as a multiple precision 
  # vector of length 100
  #
  # Written by P.Kabaila in Sept 2021
  
  len.Lagrange.interp.poly.mp <- length(Lagrange.interp.poly.mp)
  
  if (n > 50){
    stop("n is too large in sq_Lagrange_interp_poly_mp_rev1.R")
  }
  
  poly.sq.mp <- mpfr(rep(0, len.Lagrange.interp.poly.mp), nbits)
  
  for (i in c(1:len.Lagrange.interp.poly.mp)){
    for (j in c(1:len.Lagrange.interp.poly.mp)){
      if ((i + j - 1) <= len.Lagrange.interp.poly.mp)
        poly.sq.mp[i + j - 1] <- poly.sq.mp[i + j - 1] + 
          Lagrange.interp.poly.mp[i] * Lagrange.interp.poly.mp[j]
    }
  }
  
  poly.sq.mp
  
}
