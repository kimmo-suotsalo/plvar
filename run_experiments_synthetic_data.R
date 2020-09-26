# run_experiments_synthetic_data.R
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

# Script for running GVAR simulation in Puhti

d <- 80

library(foreach)
library(doParallel)

#setup parallel backend to use all assigned cores
cores <- detectCores()
cl <- makeCluster(cores[1])
registerDoParallel(cl)


run_iter <- function(iter, d){
  
  tic1 <- Sys.time()
  
  .libPaths(c("/projappl/corander/project_rpackages",.libPaths()))  

  library(SparseTSCGM)
  library(longitudinal)
  source("gvar_learn_structure.R")
  
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
  ks <- c(1,2)
  ns <- c(50,100,200,400,800)
  
  Q_hd <- array(0, dim = c(3, length(ns), length(ks)))
  Q_prec <- array(0, dim = c(3, length(ns), length(ks)))
  Q_rec <- array(0, dim = c(3, length(ns), length(ks)))
  
  B_hd <- array(0, dim = c(3, length(ns), length(ks)))
  B_prec <- array(0, dim = c(3, length(ns), length(ks)))
  B_rec <- array(0, dim = c(3, length(ns), length(ks)))
  
  run_time <- array(0, dim = c(3,length(ns), length(ks)))
  k_estimated <- array(0, dim = c(3,length(ns), length(ks)))
  
  for (k in ks){
      l <- 0
      cont <- TRUE
      while(cont){
        set.seed(1337*iter+l)
        mod <- sim.data(model = paste("ar", k, sep =""),
                        time = max(ns),
                        n.obs = 2, 
                        n.var = d,
                        prob0 = (3*d)/(k*(d^2)), 
                        network = "random")
        data <- mod$data1[c(TRUE,FALSE),]
        # To exclude unstable models
        print(max(abs(data)))
	if (k == 1){
  	  A <- t(mod$gamma)
	} else {
  	  A <- rbind(t(mod$gamma), diag(d*k)[1:(d*(k-1)),])
	}
	if (any(abs(eigen(A)$values)>=1)){
  	  l <- l+1
	} else {
	  cont <- FALSE
	}
      }
      B_amat_true <- t(mod$gamma != 0)
      Q_amat_true <- (mod$theta != 0)-diag(d)
      
      for (j in 1:length(ns)){
        print(paste("Iteration: ",iter," - Lag: ", k, " - Number of vars: ",d," - Sample size: ",ns[j],sep = ""))
        tic2 <- Sys.time()
        cg <- gvar_learn_structure(data[1:ns[j],], K = 5, gamma = 0.5)
        run_time[1,j,k] <- difftime(Sys.time(), tic2, units = "secs")
        k_est <- cg$est_k
        k_estimated[1,j,k] <- k_est        

        Q_hd[1,j,k] <- sum(Q_amat_true != cg$Q_amat[,,k_est])/2
        Q_prec[1,j,k] <- sum(Q_amat_true+cg$Q_amat[,,k_est] == 2)/sum(cg$Q_amat[,,k_est] != 0)
        Q_rec[1,j,k] <- sum(Q_amat_true+cg$Q_amat[,,k_est] == 2)/sum(Q_amat_true != 0)
        if (k_est == k){
          B_hd[1,j,k] <- sum(B_amat_true != cg$B_amat[,1:(d*k_est),k_est])
          B_prec[1,j,k] <- sum(B_amat_true+cg$B_amat[,1:(d*k_est),k_est] == 2)/sum(cg$B_amat[,1:(d*k_est),k_est] != 0)
          B_rec[1,j,k] <- sum(B_amat_true+cg$B_amat[,1:(d*k_est),k_est] == 2)/sum(B_amat_true != 0)
        } else if (k_est > k) {
          B_hd[1,j,k] <- sum(B_amat_true != cg$B_amat[,1:(d*k),k_est])+sum(cg$B_amat[,(d*k+1):(d*k_est),k_est] != 0)
          B_prec[1,j,k] <- sum(B_amat_true+cg$B_amat[,1:(d*k),k_est] == 2)/sum(cg$B_amat[,1:(d*k_est),k_est] != 0)
          B_rec[1,j,k] <- sum(B_amat_true+cg$B_amat[,1:(d*k),k_est] == 2)/sum(B_amat_true != 0)
        } else {
          B_hd[1,j,k] <- sum(B_amat_true[,1:(d*k_est)] != cg$B_amat[,1:(d*k_est),k_est])+sum(B_amat_true[,(d*k_est+1):(d*k)] != 0)
          B_prec[1,j,k] <- sum(B_amat_true[,1:(d*k_est)]+cg$B_amat[,1:(d*k_est),k_est] == 2)/sum(cg$B_amat[,1:(d*k_est),k_est] != 0)
          B_rec[1,j,k] <- sum(B_amat_true[,1:(d*k_est)]+cg$B_amat[,1:(d*k_est),k_est] == 2)/sum(B_amat_true != 0)
        }
        
        tic2 <- Sys.time()
        cg <- tryCatch(quiet(sparse.tscgm(data = as.longitudinal(data[1:ns[j],]), 
                                 lam1 = NULL, 
                                 lam2 = NULL, 
                                 nlambda = NULL, 
                                 model= paste("ar", k, sep =""), 
                                 penalty = "scad",
                                 optimality="bic_mod",
                                 control=list(maxit.out = 10, maxit.in = 100, silent = TRUE))),
                       error = function(e) {print(paste("Error in optimization")); NaN})
        run_time[2,j,k] <- difftime(Sys.time(), tic2, units = "secs")

        if (class(cg) != "sparse.tscgm"){
          Q_hd[2,j,k] <- NaN
          Q_prec[2,j,k] <- NaN
          Q_rec[2,j,k] <- NaN
          B_hd[2,j,k] <- NaN
          B_prec[2,j,k] <- NaN
          B_rec[2,j,k] <- NaN
          run_time[2,j,k] <- NaN
        } else {
          Q_amat <- ((cg$theta-diag(d)) != 0)
          Q_hd[2,j,k] <- sum(Q_amat_true != Q_amat)/2
          Q_prec[2,j,k] <- sum(Q_amat_true+Q_amat == 2)/sum(Q_amat != 0)
          Q_rec[2,j,k] <- sum(Q_amat_true+Q_amat == 2)/sum(Q_amat_true != 0)
          B_amat <- t(cg$gamma != 0)
          B_hd[2,j,k] <- sum(B_amat_true != B_amat)
          B_prec[2,j,k] <- sum(B_amat_true+B_amat == 2)/sum(B_amat != 0)
          B_rec[2,j,k] <- sum(B_amat_true+B_amat == 2)/sum(B_amat_true != 0)
        }
        
        tic2 <- Sys.time()
        cg <- tryCatch(quiet(sparse.tscgm(data = as.longitudinal(data[1:ns[j],]), 
                                          lam1 = NULL, 
                                          lam2 = NULL, 
                                          nlambda = NULL, 
                                          model= paste("ar", k, sep =""), 
                                          penalty = "lasso",
                                          optimality="bic_mod",
                                          control=list(maxit.out = 10, maxit.in = 100, silent = TRUE))),
                       error = function(e) {print(paste("Error in optimization")); NaN})
        run_time[3,j,k] <- difftime(Sys.time(), tic2, units = "secs")
        
        if (class(cg) != "sparse.tscgm"){
          Q_hd[3,j,k] <- NaN
          Q_prec[3,j,k] <- NaN
          Q_rec[3,j,k] <- NaN
          B_hd[3,j,k] <- NaN
          B_prec[3,j,k] <- NaN
          B_rec[3,j,k] <- NaN
          run_time[3,j,k] <- NaN
        } else {
          Q_amat <- ((cg$theta-diag(d)) != 0)
          Q_hd[3,j,k] <- sum(Q_amat_true != Q_amat)/2
          Q_prec[3,j,k] <- sum(Q_amat_true+Q_amat == 2)/sum(Q_amat != 0)
          Q_rec[3,j,k] <- sum(Q_amat_true+Q_amat == 2)/sum(Q_amat_true != 0)
          B_amat <- t(cg$gamma != 0)
          B_hd[3,j,k] <- sum(B_amat_true != B_amat)
          B_prec[3,j,k] <- sum(B_amat_true+B_amat == 2)/sum(B_amat != 0)
          B_rec[3,j,k] <- sum(B_amat_true+B_amat == 2)/sum(B_amat_true != 0)
        }
      }
    
  }
  print(Sys.time()-tic1)
  return(list("B_hd" = B_hd, "B_prec" = B_prec, "B_rec" = B_rec, "Q_hd" = Q_hd, "Q_prec" = Q_prec, "Q_rec" = Q_rec, "run_time" = run_time, "k_estimated" = k_estimated))
}


res <- foreach (iter = 1:20, .errorhandling = "pass") %dopar% {
  run_iter(iter, d)
}


# Save results
saveRDS(res, file = paste("sim_results_d", d, ".rds", sep = ""))

# close cluster
stopCluster(cl)
