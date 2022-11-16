# test_commands.R

# Test the module moments.R

# m <- 160
# which.f <- list(name="scaled.chi.pdf", support=c(0, Inf), parameters=m)
# moments(which.f, 4, 120)
#
# # Test the module alpha_mp_rev2.R
#
# m <- 160
# which.f <- list(name="scaled.chi.pdf", support=c(0, Inf), parameters=m)
# nbits <- 160
# k <- 5
#
# alpha_mp_rev2(which.f, k, nbits)
# The result in the pre-package version is the following:
# 1 'mpfr' number of precision  160   bits
# [1] 1.0141101898090958814506055050387277378354968732815

# Test the module beta_mp_rev2.R

# m <- 160
# which.f <- list(name="scaled.chi.pdf", support=c(0, Inf), parameters=m)
# nbits <- 160
# k <- 5
#
# beta_mp_rev2(which.f, k, nbits)
# The result in the pre-package version is the following:
# 1 'mpfr' number of precision  160   bits
# [1] 0.015504086930667336743019744252084960132933473630352

#-------------------------
# Repeat an example of a check against known results
# Gauss Hermite quadrature

# which.f <- list(name="Hermite",
#                  support=c(-Inf, Inf), parameters=NULL)
#
# n <- 3

# The following commands result in an error if
# Build > Load All is used.
# Consequently, I comment out these commands when
# using Build > Load All and then uncomment it
# when I want to carry out the test computations.

# system.time(output.name <- custom(which.f, n))
#
#  output.name$mat.timings
#  options(digits = 2)
#  output.name$max.abs.diffs.nodes
#  output.name$sum.abs.diffs.weights
#  output.name$L.nodes
#  output.name$L.weights
#
#  len.nbits.vec <- length(output.name$nbits.vec)
#
# # The numbers in blue come from Keisan Casio website
# # (with 50 digits). I copied these using ctrl-C and
# # pasted them using ctrl-V
#
# output.name$list.Gauss.nodes[[len.nbits.vec]][1]
#   -1.2247448713915890490986420373529456959829737403283
#
# output.name$list.Gauss.nodes[[len.nbits.vec]][2]
#   0
#
# output.name$list.Gauss.weights[[len.nbits.vec]][1]
#   0.2954089751509193378830279138901908637995915760204
#
# output.name$list.Gauss.weights[[len.nbits.vec]][2]
#   1.1816359006036773515321116555607634551983663040816
#

#------------------------------

# m <- 160
# which.f <- list(name="scaled.chi.pdf", support=c(0, Inf), parameters=m)
# n <- 1
# nbits.first <- ceiling(60 + 6.5 * n)
# nbits.vec <- seq(nbits.first, nbits.first + 136, length.out = 5)
# len.nbits.vec <- length(nbits.vec)
#
# cat("nbits.vec = (", nbits.vec[1], ",", nbits.vec[2], ",",
#     nbits.vec[3], ",", nbits.vec[4], ",", nbits.vec[5],")", "\n")
#
# if (n==1){
#   list.Gauss.nodes <- NULL
#   list.Gauss.weights <- NULL
#   for (j in c(1:len.nbits.vec)){
#     nbits <- nbits.vec[j]
#     mu0 <- moments(which.f, 0, nbits)
#     mu1 <- moments(which.f, 1, nbits)
#     Gauss.nodes <- mu1 / mu0
#     list.Gauss.nodes <- append(list.Gauss.nodes, Gauss.nodes)
#     Gauss.weights <- mu0
#     list.Gauss.weights <- append(list.Gauss.weights, Gauss.weights)
#   }
# }
#
# list.Gauss.nodes
# list.Gauss.weights
#
# list.Gauss.nodes[[5]]
#
# system.time(gauss.list <- custom(which.f, n))
# options(digits = 17)
# gauss.list
#
#
# m <- 160
# which.f <- list(name="scaled.chi.pdf", support=c(0, Inf), parameters=m)
# n <- 1
# old <- options(digits = 17)
# system.time(gauss.list <- custom(which.f, n))
# gauss.list
# options(old)
# options("digits")
#
# m <- 160
# which.f <- list(name="scaled.chi.pdf", support=c(0, Inf), parameters=m)
# n <- 1
# system.time(gauss.list <- custom(which.f, n))
# old <- options(digits = 17)
# gauss.list
#
#
#
# library(statmod)
#
# alpha.GGL <- 0
# n <- 1
# which.f <- list(name="Generalized.Laguerre", support=c(0, Inf),
#                 parameters=alpha.GGL)
# gauss.list <- custom(which.f, n)
# temp <- gauss.quad(n, kind = "laguerre", alpha = alpha.GGL)
#
# gauss.list$asNumeric.nodes
# temp$nodes
#
# gauss.list$asNumeric.weights
# temp$weights
#
#
#
# alpha.GGL <- 1
# n <- 1
# which.f <- list(name="Generalized.Laguerre", support=c(0, Inf),
#                 parameters=alpha.GGL)
# gauss.list <- custom(which.f, n)
# temp <- gauss.quad(n, kind = "laguerre", alpha = alpha.GGL)
#
# gauss.list$asNumeric.nodes
# temp$nodes
#
# gauss.list$asNumeric.weights
# temp$weights
#
#
# alpha.GGL <- 2
# n <- 1
# which.f <- list(name="Generalized.Laguerre", support=c(0, Inf),
#                 parameters=alpha.GGL)
# gauss.list <- custom(which.f, n)
# temp <- gauss.quad(n, kind = "laguerre", alpha = alpha.GGL)
#
# gauss.list$asNumeric.nodes
# temp$nodes
#
# gauss.list$asNumeric.weights
# temp$weights
#
# old <- options(digits = 17)
# alpha.GGL <- -1/2
# n <- 1
# which.f <- list(name="Generalized.Laguerre", support=c(0, Inf),
#                 parameters=alpha.GGL)
# gauss.list <- custom(which.f, n)
# temp <- gauss.quad(n, kind = "laguerre", alpha = alpha.GGL)
#
# gauss.list$asNumeric.nodes
# temp$nodes
#
# gauss.list$asNumeric.weights
# temp$weights
#
# alpha.GGL <- -1/4
# n <- 1
# which.f <- list(name="Generalized.Laguerre", support=c(0, Inf),
#                 parameters=alpha.GGL)
# gauss.list <- custom(which.f, n)
# temp <- gauss.quad(n, kind = "laguerre", alpha = alpha.GGL)
#
# gauss.list$asNumeric.nodes
# temp$nodes
#
# gauss.list$asNumeric.weights
# temp$weights
#
# alpha.GGL <- -1/8
# n <- 1
# which.f <- list(name="Generalized.Laguerre", support=c(0, Inf),
#                 parameters=alpha.GGL)
# gauss.list <- custom(which.f, n)
# temp <- gauss.quad(n, kind = "laguerre", alpha = alpha.GGL)
#
# gauss.list$asNumeric.nodes
# temp$nodes
#
# gauss.list$asNumeric.weights
# temp$weights
#
#
#
# options(old)
#

# which.f <- list(name="chemistry.example", support=c(0, Inf),
#                                  parameters=NULL)
# n <- 15
#
# old <- options(digits = 17)
# system.time(gauss.list <- custom(which.f, n))
# gauss.list
#
#
#
#
#
#
#
# options(old)
