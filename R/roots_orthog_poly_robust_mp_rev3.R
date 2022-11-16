roots_orthog_poly_robust_mp_rev3 <- 
  function(which.f, orthog.poly.mp, n, nbits){
    # This module evaluates the roots, to within tol, 
    # of the monic orthogonal polynomial of degree n, 
    # where n is a positive integer.
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
    # Inputs:
    # which.f: a list specifying f     
    # orthog.poly.mp: vector of length (n + 1)
    #       consisting of the coefficients in multiple
    #       precision of the orthogonal polynomial in 
    #       increasing order. 
    # n: number of Gauss quadrature nodes
    # nbits: number of bits in the multiple precision
    #        numbers
    #
    # Output:
    # A vector of length n with multiple precision 
    # elements that consists of the roots of 
    # the monic orthogonal polynomial.
    #
    # Written by P.Kabaila in Sept 2022
    
    n.test <- length(orthog.poly.mp) - 1
    if (n.test != n){
      stop(
        "length of orthog.poly.mp incorrect in roots_orthog_poly_robust_mp_rev2")
    }
    
    Laguerre.bounds.mp <- Laguerre_bounds_mp(orthog.poly.mp)
    
    # According to Theorem 1.46 of Gautschi (2004),
    # the nodes are contained in the interior of the
    # support interval. 
    support.vec <- which.f$support
    support.low <- support.vec[1]
    support.up <- support.vec[2]
    
    if (support.low == -Inf){
      low.mp <- Laguerre.bounds.mp[1]
    }else{
      low.mp <- max(mpfr(support.low, nbits), Laguerre.bounds.mp[1])
    }
    
    if (support.up == Inf){
      up.mp <- Laguerre.bounds.mp[2]
    }else{
      up.mp <- min(mpfr(support.up, nbits), Laguerre.bounds.mp[2])
    }
    
    epsilon <- (up.mp - low.mp) / mpfr(100000000, nbits)
    low.mp <- low.mp - epsilon
    up.mp <- up.mp + epsilon
    
    zero.mp <- mpfr(0, nbits)
    
    # Code check
    # cat("asNumeric lower bound=", asNumeric(low.mp),
    #   ",  asNumeric upper bound=", asNumeric(up.mp))
    
    multiplier <- 100
    grid.vec.mp <- seq(low.mp, up.mp, length.out = multiplier * n)
    
    orthog.poly.mp.vec <- mpfr(rep(0, multiplier*n), nbits)
     
    for (i in c(1: (multiplier*n))){
      orthog.poly.mp.vec[i] <- 
        eval_orthog_poly_mp(grid.vec.mp[i], orthog.poly.mp)
    }
    
    # Code check
    # grid.vec <- asNumeric(grid.vec.mp)
    # orthog.poly.vec <- asNumeric(orthog.poly.mp.vec)
    # plot(grid.vec, orthog.poly.vec, type="l",
    #    xlab="x", ylab="orthog poly(x)", ylim=c(-2,2))
    # abline(h=0, col="red")
    # mtext(paste("n=", n))
    
    # Code check
    # browser()
    
    lower.endpt.indices.vec <- NULL
    for (i in c(1: (multiplier*n-1))){
      if (orthog.poly.mp.vec[i] * orthog.poly.mp.vec[i+1] 
          < mpfr(0, nbits)){
        lower.endpt.indices.vec <- c(lower.endpt.indices.vec, i)
      }
    }
    
    if(length(lower.endpt.indices.vec) < n){
      stop("Too few intervals bracketing roots")
    }
    
    tol <- 2^(-nbits)
    roots.vec <- mpfr(rep(0, n), nbits)
    for (i in c(1:n)){
      lower.mp <- grid.vec.mp[lower.endpt.indices.vec[i]]
      upper.mp <- grid.vec.mp[lower.endpt.indices.vec[i] + 1]
      temp <- unirootR(eval_orthog_poly_mp, 
                       c(lower.mp, upper.mp),
                       orthog.poly.mp, tol = tol)
      roots.vec[i] <- temp$root
    }
    
    roots.vec
    
  }