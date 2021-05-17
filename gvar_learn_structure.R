# gvar_learn_structure.R
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
# Learns a sparse Gaussian GVAR structure in a two-step procedure based on the FM(P)L:
# --- Step 1: Identify lag-parents using the FML + an eBIC-type prior in combination with a greedy 
#             hill-climbing algorithm. The lag length (k) is estimated as the k under which the lag 
#             structure score is maximized.
# --- Step 2: Estimate lag-effect parameters using restricted least squares and calculate the residuals. 
#             Estimate the sparsity pattern of precision matrix using the FMPL + an eBIC-type prior in 
#             combination with a greedy hill-climbing algorithm.
# ---------------------------------------------------------------------------------------------------------
# INPUT:  - data: n-by-d data matrix where n is the sample size and d the number of variables.
#         - K: maximum lag length, the method iterates through 1:K and selects the optimal lag length.
#         - gamma: tuning parameter in the eBIC-type sparsity prior (default = 0.5).
# OUTPUT: a list with the following objects:
#         - B_amat: array containing the estimated lag structure for the different k values.
#         - B_score: scores of the estimated lag structures for the different k values.
#         - est_k: estimated k, the k value with the maximum lag structure score.
#         - Q_amat: array containing the estimated contemporaneous structures for the different k values.
#         - Q_score: scores of the estimated contemporaneous structures for the different k values.
#
# Author: Johan Pensar (1 July, 2020)
#==========================================================================================================

# Main function
gvar_learn_structure <- function(data, K = 5, gamma = 0.5, k = NULL){
  
  d <- ncol(data)
  n <- nrow(data)
  B_amat <- array(0, dim = c(d,d*K,K))
  Q_amat <- array(0, dim = c(d,d,K))
  B_score <- array(0,K)
  Q_score <- array(0,K)
  
  X <- make_X(data, K)
  X <- apply(X,2,function(x) x-mean(x))
  
  S <- t(X)%*%X
  neff <- nrow(X)
  
  if (is.null(k)) {
    k_candidates <- 1:K
  } else {
    k_candidates <- k:k
  }
  
  for (k in k_candidates){
    
    # Learn temporal structure
    B_ols <- matrix(0, d, d*k)
    for (node in 1:d){
      ind <- c(1:(k*d),(K*d)+node)
      tmp <- hc(k*d+1, S[ind,ind], neff, gamma)
      par <- tmp$par
      B_score[k] <- B_score[k]+tmp$score
      B_amat[node,par,k] <- 1
      if (length(par) > 0){
        B_ols[node, par] <- estimate_B_ols(node+K*d, par, X)
      }
    }
    
    # Learn contemporaneous structure
    res <- X[,-(1:(K*d))]-(X[,1:(k*d)]%*%t(B_ols))
    S_res <- t(res)%*%res
    for (node in 1:d){
      tmp <- hc(node, S_res, neff, gamma)
      nbr <- tmp$par
      Q_score[k] <- Q_score[k]+tmp$score
      Q_amat[node,nbr,k] <- 1
      Q_amat[nbr,node,k] <- 1
    }
  }
  return(list("Q_amat" = Q_amat, "B_amat" = B_amat, "Q_score" = Q_score, "B_score" = B_score, "est_k" = which.max(B_score)))
  #return(list("B_ols" = B_ols, "B_amat" = B_amat))
}

# Constructs the lagged data matrix
make_X <- function(data, k){
  n <- nrow(data)
  tmp <- vector("list", k+1)
  for (i in 1:k){
    tmp[[i]] <- data[(k-i+1):(n-i),]
  }
  tmp[[k+1]] <- data[(k+1):n,]
  X <- do.call(cbind, tmp)
  return(X)
}

# Finds the "optimal" parent set using a HC search
hc <- function(node, S, n, gamma){
  curr_par <- NULL
  curr_score <- calc_ml(node, curr_par, S, n, gamma)
  cand_score <- array(NaN, ncol(S))
  cand_score[node] <- -Inf
  add <- TRUE
  while(add == TRUE & length(curr_par) < (n-1)){
    ind <- which(is.nan(cand_score))
    for (i in 1:length(ind)){
      cand_score[ind[i]] <- calc_ml(node, c(curr_par,ind[i]), S, n, gamma)
    }
    pos <- which.max(cand_score)
    if (cand_score[pos] > curr_score){
      curr_par <- c(curr_par,pos)
      curr_score <- cand_score[pos]
      cand_score[pos] <- -Inf
      cand_score[cand_score != -Inf] <- NaN
      add <- TRUE
    } else {
      add <- FALSE
    }
    par_size <- length(curr_par)
    if (par_size > 2 & add == TRUE){
      del <- TRUE
      while(del == TRUE){
        del <- FALSE
        del_cand_score <- array(NaN, par_size)
        for (i in 1:par_size){
          del_cand_score[i] <- calc_ml(node, curr_par[-i], S, n, gamma)
        }
        pos <- which.max(del_cand_score)
        if (del_cand_score[pos] > curr_score){
          cand_score[curr_par[pos]] <- -Inf
          curr_par <- curr_par[-pos]
          curr_score <- del_cand_score[pos]
          if (length(curr_par) > 2){
            del <- TRUE
          } 
        }
      }
    }
  }
  return(list("par" = curr_par, "score" = curr_score))
}

# Calculates the FML of a family
calc_ml <- function(node, par, S, n, gamma){
  fam <- c(node,par)
  dpar <- length(fam)-1
  ml <- -log(pi)*(n-1)*0.5+
        lgamma(0.5*(n+dpar))-
        lgamma(0.5*(dpar+1))-
        log(n)*(dpar+0.5)-
        (log(det(S[fam,fam, drop = FALSE]))-log(det(S[par,par, drop = FALSE])))*(n-1)*0.5
  ml <- ml-dpar*log(ncol(S)-1)*gamma
  return(ml)
}

# Calculates the LS estimates of a node with some lag-parents
estimate_B_ols <- function(node, par, X){
  x <- X[,par]
  y <- X[,node]
  b <- solve(t(x)%*%x)%*%t(x)%*%y
  return(b)
}
