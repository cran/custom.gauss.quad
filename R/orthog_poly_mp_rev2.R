orthog_poly_mp_rev2 <- function(which.f, n, nbits){
  # This module computes the monic orthogonal 
  # polynomial of degree n, where n is a positive
  # integer, using the three-term recursion 
  # on page 10 of
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
  # Inputs (which are not in multiple precision)
  # which.f: a list specifying f 
  # n: degree of polynomial, a positive integer
  #    this is also the number of Gauss quadrature
  #    nodes 
  # m: degrees of freedom of chi-distribution
  # nbits: number of bits in the multiple precision
  #        numbers
  #
  # Output:
  # Monic orthogonal polynomial of degree n
  # in the form of a multiple precision vector
  # of length n+1
  #
  # Written by P.Kabaila in Sept 2022
  
  if (n > 99){
    stop("n is too large in orthog_poly_mp_rev2.R")
  }
  
  poly.oldest <- mpfr(rep(0, 100), nbits)
  poly.old <- mpfr(c(1,rep(0, 99)), nbits)
  
  alpha.0 <- alpha_mp_rev2(which.f, 0, nbits)
  beta.0 <- beta_mp_rev2(which.f, 0, nbits)
  r.shifted.poly.old <- c(mpfr(0, nbits), poly.old[c(1:99)])
  poly.new <- r.shifted.poly.old - 
    alpha.0 * poly.old - beta.0 * poly.oldest
  
  if (n == 1){
    return(poly.new[c(1:2)])
  }
  
  for (k in c(1:(n - 1))){
    
    poly.oldest <- poly.old
    poly.old <- poly.new
    
    alpha.k <- alpha_mp_rev2(which.f, k, nbits)
    beta.k <- beta_mp_rev2(which.f, k, nbits)
    r.shifted.poly.old <- c(mpfr(0, nbits), poly.old[c(1:99)])
    poly.new <- r.shifted.poly.old - 
      alpha.k * poly.old - beta.k * poly.oldest
  }
  
  poly.new[c(1:(n+1))]
  
}