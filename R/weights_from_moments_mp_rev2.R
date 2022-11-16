weights_from_moments_mp_rev2 <-
  function(which.f, roots.orthog.poly.mp, n, nbits){
    # This module computes the Gauss quadrature
    # weights using the formula in the proof
    # of Theorem 1.46 on p.23 of
    #
    # Gautschi, W. (2004) Orthogonal Polynomials:
    # Computation and Approximation,
    #
    # where the integral is evaluated using
    # moments evaluated in multiple precision.
    #
    # The nonnegative integrable function f
    # is specified by the list which.f with
    # the following 3 components:
    #  (i)  name,
    #  (ii) support specified by a 2-vector
    #    of the endpoints of the interval,
    #  (iii) parameter vector
    #
    # Inputs
    # which.f: a list specifying f
    # roots.orthog.poly.mp: vector of the roots of the orthogonal
    #                       polynomial (multiple precision)
    # m: degrees of freedom of the chi-distribution
    #    (NOT multiple precision)
    # n: number of Gauss quadrature nodes
    # nbits: number of bits in the multiple precision
    #        numbers (NOT multiple precision)
    #
    # Output
    # n-vector with the Gauss quadrature weights
    # represented as a multiple precision vector
    #
    # Written by P.Kabaila in Sept 2022

    n.test <- length(roots.orthog.poly.mp)
    if (n.test != n){
      char.str1 <- "length of roots.orthog.poly.mp is incorrect"
      char.str2 <- "in weights_from_moments_mp_rev2.R"
      message.char.str <- paste(char.str1, char.str2)
      stop(message.char.str)
    }

    Gauss.quad.wts <- mpfr(rep(0, n), nbits)

    for (j in c(1:n)){
      Lagrange.interp.poly.mp <-
        Lagrange_interp_poly_mp_rev1(j, roots.orthog.poly.mp, n, nbits)
      sq.Lagrange.interp.poly.mp <-
        sq_Lagrange_interp_poly_mp_rev1(Lagrange.interp.poly.mp, n, nbits)

      wt.mp <- mpfr(0, nbits)
      for (k in c(1: (2*n-1))){
        coeff.mp <- sq.Lagrange.interp.poly.mp[k]
        wt.mp <- wt.mp + coeff.mp *  moments(which.f, (k-1), nbits)
      }

      Gauss.quad.wts[j] <- wt.mp

    }

    Gauss.quad.wts

  }
