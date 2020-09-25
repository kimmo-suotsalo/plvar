# gvar_estimate_parameters.R
# Copyright (C) 2020 Johan Pensar
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/.

#==========================================================================================================
# Calculates the MLEs of a sparse Gaussian GVAR in an iterative procedure:
# --- Iteration 1: Calculate the MLEs of the lag parameters given a precision/covariance matrix.
# --- Iteration 2: Calculate the MLE of the sparse precision matrix given the lag parameters.
# The method iterates between 1 and 2 until convergence, or, in practice, until the updated after going 
# through 1 and 2 the likelihood is close enough to the previous likelihood. For me details about the
# method, see Lutkepohl (2005), New Introduction to Multiple Time Series Analysis.
# ---------------------------------------------------------------------------------------------------------
# INPUT:  - data: n-by-d data matrix where n is the sample size and d the number of variables.
#         - B_amat: d-by-(k*d) matrix containing the lag structure (k = lag length, d = number of var).
#         - Q_amat: d-by-d matrix containing the contemporaneous structure.
#         - tol: tolerance threshold specifying convergence of the method (default = 1e-6).
# OUTPUT: a list with the following objects:
#         - B_hat: estimated lag matrix.
#         - Q_hat: estimated precision matrix.
#         - log_lh: log-likelihood at convergence.
#
# NOTE: the function uses the "mixggm" package for fitting the sparse precision matrix.
#
# Author: Johan Pensar (1 July, 2020)
#==========================================================================================================

# Load required package
library("mixggm")

# Main function
gvar_estimate_parameters <- function(data, B_amat, Q_amat, tol = 1e-6){
  n <- nrow(data)
  d <- ncol(data)
  k <- ncol(B_amat)/d
  d <- ncol(data)
  R <- construct_R(B_amat, k)
  Z <- construct_Z(data, k)
  y <- construct_y(data, k)
  Y <- t(data[-(1:k),])
  Qhat <- solve(diag(d))
  
  diff = Inf
  curr_ll = -Inf
  while(abs(diff) > tol){
    alpha_hat <- R%*%solve(t(R)%*%((Z%*%t(Z))%x%Qhat)%*%R)%*%t(R)%*%(Z%x%Qhat)%*%y
    B <- matrix(alpha_hat,d,d*k)
    res <- (B%*%Z-Y)
    tmp <-fitGGM(data = t(res),
                 graph = Q_amat,
                 model = "concentration",
                 regularize = FALSE,
                 verbose = FALSE)
    Qhat <- tmp$omega
    new_ll <- tmp$loglik
    diff <- new_ll-curr_ll
    print(diff)
    curr_ll <- new_ll
  }
  return(list("B_hat" = B, "Q_hat" = Qhat, "loglh" = curr_ll))
}

# Constructs R matrix (encoding the sparsity restrictions)
construct_R <- function(A, k){
  d <- nrow(A)
  alpha <- unlist(A)
  ind <- which(alpha == TRUE)
  m <- length(ind)
  R <- matrix(0, k*(d^2), m) 
  pos <- 1
  for (i in 1:m){
    R[ind[i],pos] <- 1
    pos <- pos+1
  }
  return(R)
}

# Constructs Z matrix (lagged data matrix)
construct_Z <- function(data, k){
  n <- nrow(data)
  tmp <- vector("list",k)
  for (i in 1:k){
    tmp[[i]] <- data[(k+1-i):(n-i),]   
  }
  L <- t(do.call(cbind,tmp))
  return(L)
}

# Constructs y vector (non-lagged data matrix in vector form)
construct_y <- function(data, k){
  y <- as.vector(t(data[-(1:k),]))
  return(y)
}