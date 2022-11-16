det_LU_Dolittle <- function(A, nbits){
  # This module computes the determinant of
  # the square matrix A whose elements are
  # Rmpfr numbers.
  # It does this by first computing the LU
  # decomposition of A, using 
  # a translation to Rmpfr numbers of the 
  # module lu.decomposition.R from the R
  # package matrixcalc which implements a 
  # Dolittle LU decomposition.
  #
  # Multiple precision is used with nbits.
  #
  # Inputs
  # A: square matrix with n rows
  # nbits: number of bits in the multiple precision
  #        numbers
  # 
  # Output
  # A list whose first element is L and whose
  # second element is U.
  #
  # Written by P.Kabaila in Aug 2022
  
  tmp.vec <- dim(A)
  if (tmp.vec[1] != tmp.vec[2]){
    stop("Matrix A is not square in det_LU_Dolittle.R")
  }
  n <- tmp.vec[1]
  
  U <- mpfr(rep(0, n^2), nbits)
  dim(U) <- c(n, n)
  
  L <- mpfr(rep(0, n^2), nbits)
  dim(L) <- c(n, n)
  for (i in c(1:n)){
    for (j in c(1:n)){
      if (i == j){
        L[i, j] <- mpfr(1, nbits)
      }
    }
  }
  
  for (i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for (j in 1:n) {
      U[i, j] <- A[i, j]
      if (im1 > 0) {
        for (k in 1:im1) {
          U[i, j] <- U[i, j] - L[i, k] * U[k, j]
        }
      }
    }
    if (ip1 <= n) {
      for (j in ip1:n) {
        L[j, i] <- A[j, i]
        if (im1 > 0) {
          for (k in 1:im1) {
            L[j, i] <- L[j, i] - L[j, k] * U[k, i]
          }
        }
        if (U[i, i] == 0) 
          stop("Matrix A is singular in det_LU_Dolittle.R")
        L[j, i] <- L[j, i] / U[i, i]
      }
    }
  }
  
  prod.cum <- mpfr(1, nbits)
  for (i in 1:n){
    prod.cum <- prod.cum * L[i, i] * U[i,i]
  }
  
  prod.cum
  
} 