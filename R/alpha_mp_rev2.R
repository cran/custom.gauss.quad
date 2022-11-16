alpha_mp_rev2 <- function(which.f, k, nbits){
  # This module computes alpha_k as defined
  # by (2.1.4) on p.54 of 
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
  # k: index of alpha_k, a nonnegative integer
  # m: degrees of freedom
  # nbits: number of bits in the multiple precision
  #        numbers
  #
  # Output:
  # alpha_k in multiple precision.
  #
  # Written by P.Kabaila in Sept 2022
  
  Delta.kp1.prime <- Hankel_prime_det_mp_rev2(which.f, (k+1), nbits)
  Delta.kp1 <- Hankel_det_mp_rev2(which.f, (k+1), nbits)
  term1 <- Delta.kp1.prime / Delta.kp1
  
  Delta.k.prime <- Hankel_prime_det_mp_rev2(which.f, k, nbits)
  Delta.k <- Hankel_det_mp_rev2(which.f, k, nbits)
  term2 <- Delta.k.prime / Delta.k
  
  term1 - term2
  
}
