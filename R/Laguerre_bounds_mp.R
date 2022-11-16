Laguerre_bounds_mp <- function(orthog.poly.mp){
  # This module computes the Laguerre bounds
  # on the roots of the monic orthogonal polynomial
  # orthog.poly.mp of degree n, where n is a 
  # positive integer. 
  #
  # It includes an error message if b
  # is the square root of a negative 
  # number.
  #
  # multiple precision is used with nbits
  # 
  # Input:
  # orthog.poly.mp: monic polynomial orthogonal
  #                 polynomial
  #
  # Output:
  # A 2-vector with the lower and upper
  # Laguerre bounds in multiple precision.
  # 
  # Written by P.Kabaila in Aug 2022
  
  # n is the degree of orthog.poly.mp
  n <- length(orthog.poly.mp) - 1
  
  # We use the notation a1, a2 and b
  # of Jensen and Styan (1999)
  a1 <- orthog.poly.mp[n]
  a2 <- orthog.poly.mp[n-1]
  centre <- - a1 / n
  tmp <- (n - 1) * a1^2 - 2 * n * a2
  if (tmp < 0){
    stop("b is imaginary")
  }
  
  b <- sqrt(tmp) / n
  half.width <- b * sqrt(n - 1)
  
  c(centre - half.width,
    centre + half.width)  
  
}