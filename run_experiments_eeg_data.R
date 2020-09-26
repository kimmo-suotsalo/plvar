# run_experiments_eeg_data.R
# Copyright (C) 2020 Kimmo Suotsalo
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

library(edf)
library(vars)
library(SparseTSCGM)
library(longitudinal)

source("gvar_learn_structure.R")
source("gvar_estimate_parameters.R")

# Set the pathname of the output file
output_pathname <- "./chbmit_output.csv"

# Set the names of the estimation methods
estimation_methods <- c("PLVAR", "OLS", "PLVAR2", "OLS2", "SCAD2", "LASSO2")

# Set the filenames of the EEG recordings
recording_filenames <- c(
  "chb01/chb01_01.edf", "chb02/chb02_01.edf", "chb03/chb03_01.edf",
  "chb04/chb04_01.edf", "chb05/chb05_01.edf", "chb06/chb06_01.edf",
  "chb07/chb07_01.edf", "chb08/chb08_02.edf", "chb09/chb09_01.edf",
  "chb10/chb10_01.edf", "chb11/chb11_01.edf", "chb12/chb12_06.edf",
  "chb13/chb13_02.edf", "chb14/chb14_01.edf", "chb15/chb15_02.edf",
  "chb16/chb16_01.edf", "chb17/chb17a_03.edf", "chb18/chb18_02.edf",
  "chb19/chb19_02.edf", "chb20/chb20_01.edf", "chb22/chb22_01.edf",
  "chb23/chb23_06.edf"
)

# Set the names of the EEG channels
channel_names <- c(
  "FP1_F7", "F7_T7", "T7_P7", "P7_O1",
  "FP1_F3", "F3_C3", "C3_P3", "P3_O1",
  "FP2_F4", "F4_C4", "C4_P4", "P4_O2",
  "FP2_F8", "F8_T8", "T8_P8", "P8_O2",
  "FZ_CZ",  "CZ_PZ",
  "T7_FT9", "FT9_FT10", "FT10_T8"
)

# Get the number of variables in the VAR model
d <- length(channel_names)

# Set an offset at the beginning of the recording to exclude the first 60 s
# where the signal quality may be low; the sampling frequency is 256 Hz, hence
# we set the offset to 15360 samples
offset <- 15360

# Set the number of samples to be used for training
n_train <- 256

# Set the number of samples to be used for testing
n_test <- 512

# Write a header line in the output file
sink(output_pathname, append = TRUE)
cat(
  "\n",
  formatC(
    paste("recording_filename", ",", sep=""), width = 22, format = "s"
  )
)
for (estimation_method in estimation_methods) {
  cat(
    formatC(
      paste("t_", estimation_method, ",", sep=""), width = 13, format = "s"
    ),
    formatC(
      paste("k_", estimation_method, ",", sep=""), width = 11, format = "s"
    ),
    formatC(
      paste("nd_", estimation_method, ",", sep=""), width = 11, format = "s"
    ),
    formatC(
      paste("nu_", estimation_method, ",", sep=""), width = 11, format = "s"
    ),
    formatC(
      paste("MSE_", estimation_method, sep=""), width = 12, format = "s"
    ),
    sep = ""
  )
  if (estimation_method == estimation_methods[length(estimation_methods)]) {
    cat("\n")
  } else {
    cat(",")
  }
}
sink()

# Process the EEG recordings one by one
for (recording_filename in recording_filenames) {
  
  # Read the EDF file
  cat("Reading file ", recording_filename, " ... ", sep = "")
  recording <- read.edf(
    paste("./physionet.org/files/chbmit/1.0.0/", recording_filename, sep = "")
  )
  cat("Done.", "\n\n")
  
  # Get the signals that correspond to the selected channels
  signals <- recording$signal[channel_names]
  if (length(signals[!is.na(names(signals))]) != length(channel_names)) {
    stop("Couldn't find all the selected EEG channels")
  }
  
  # Initialize a matrix for both the training and the test data
  n = n_train + n_test
  Y <- matrix(nrow = n, ncol = d)
  colnames(Y) <- channel_names

  # Fill in the matrix columns with n samples from the selected signals
  for (j in 1:d) {
    Y[,j] <- signals[[j]]$data[offset + (1:n)]
  }
  
  # Fit a linear model to each channel and subtract it to obtain de-trended,
  # zero-centered data
  x <- 1:n
  mdl <- lm(Y ~ x)
  for (j in 1:d) {
    slope <- mdl$coefficients[2,j]
    intercept <- mdl$coefficients[1,j]
    Y[,j] <- Y[,j] - (slope * x + intercept)
  }
  
  # Split the data into a training set and a test set
  Y_train <- Y[1:n_train,]
  Y_test <- Y[(n_train + 1):n,]
  
  # Write the name of the current recording to the output file
  sink(output_pathname, append = TRUE)
  cat(formatC(recording_filename, width = 22, format = "s"), ",", sep = "")
  sink()
  
  # Run the estimation methods one by one
  for (estimation_method in estimation_methods) {
    
    cat("Running", estimation_method, "... ")
    
    # Initialize a matrix for the predicted data
    Y_pred <- matrix(nrow = n_test, ncol = d)
    colnames(Y_pred) <- channel_names
    
    if (estimation_method == "PLVAR") {
      
      t_start <- Sys.time()
      
      # Learn the graph structure
      gvar <- gvar_learn_structure(Y_train, K = 30, k = NULL)
      
      # Estimate the lag length
      min_score <- min(gvar$B_score)
      max_score <- max(gvar$B_score)
      thresh_score <- min_score + 0.9 * (max_score - min_score)
      k_est <- min(which(gvar$B_score > thresh_score))
      
      # Extract the adjacency matrices that correspond to the estimated lag
      # length
      B_amat <- gvar$B_amat[,1:(k_est*d),k_est]
      Q_amat <- gvar$Q_amat[,,k_est]
      
      # Estimate the model parameters using the graph structure implied by the
      # adjacency matrices 
      gvar_mod <- gvar_estimate_parameters(Y_train, B_amat, Q_amat)
      B_hat <- gvar_mod$B_hat
      C_hat <- 0
      
      # Measure the elapsed time
      t_diff <- difftime(Sys.time(), t_start, units = "secs")
      
      # Calculate the number of directed edges in the estimated model
      nd <- sum(B_amat)
      
      # Calculate the number of undirected edges in the estimated model
      nu <- sum(Q_amat)
      
      B_amat_PLVAR <- B_amat
      Q_amat_PLVAR <- Q_amat
      
    } else if (estimation_method == "OLS") {
      
      t_start <- Sys.time()
      
      # Calculate AICs for different lag lengths
      AICs <- VARselect(Y_train, lag.max = 30)$criteria["AIC(n)",]
      
      # Convert -Inf and NaN values to very small integers
      AICs[AICs == -Inf] <- -1e300
      AICs[is.nan(AICs)] <- -1e300
      
      # Estimate the lag length
      min_score <- min(AICs)
      max_score <- max(AICs)
      thresh_score <- min_score + 0.9 * (max_score - min_score)
      k_est <- min(which(AICs < thresh_score))
      
      # Estimate the model parameters
      mdl <- VAR(y = Y_train, p = k_est)
      B_hat <- Bcoef(mdl)[,1:(k_est*d)]
      C_hat <- Bcoef(mdl)[,k_est*d+1]
      
      # Measure the elapsed time
      t_diff <- difftime(Sys.time(), t_start, units = "secs")
      
      # Calculate the number of directed edges in the estimated model
      nd <- sum(B_hat != 0)
      
      # Calculate the number of undirected edges in the estimated model
      nu <- nd
      
    } else if (estimation_method == "PLVAR2") {
      
      t_start <- Sys.time()
      
      # Learn the graph structure
      gvar <- gvar_learn_structure(Y_train, K = 30, k = 2)
      
      # Set the lag length
      k_est <- 2
      
      # Extract the adjacency matrices that correspond to the lag length
      B_amat <- gvar$B_amat[,1:(k_est*d),k_est]
      Q_amat <- gvar$Q_amat[,,k_est]
      
      # Estimate the model parameters using the graph structure implied by the
      # adjacency matrices
      gvar_mod <- gvar_estimate_parameters(Y_train, B_amat, Q_amat)
      B_hat <- gvar_mod$B_hat
      C_hat <- 0
      
      # Measure the elapsed time
      t_diff <- difftime(Sys.time(), t_start, units = "secs")
      
      # Calculate the number of directed edges in the estimated model
      nd <- sum(B_amat)
      
      # Calculate the number of undirected edges in the estimated model
      nu <- sum(Q_amat)
      
      B_amat_PLVAR2 <- B_amat
      Q_amat_PLVAR2 <- Q_amat
      
    } else if (estimation_method == "OLS2") {
      
      t_start <- Sys.time()
      
      # Set the lag length
      k_est <- 2
      
      # Estimate the model parameters
      mdl <- VAR(y = Y_train, p = k_est)
      B_hat <- Bcoef(mdl)[,1:(k_est*d)]
      C_hat <- Bcoef(mdl)[,k_est*d+1]
      
      # Measure the elapsed time
      t_diff <- difftime(Sys.time(), t_start, units = "secs")
      
      # Calculate the number of directed edges in the estimated model
      nd <- sum(B_hat != 0)
      
      # Calculate the number of undirected edges in the estimated model
      nu <- nd
      
    } else if (estimation_method == "SCAD2") {
      
      t_start <- Sys.time()
      
      # Learn the graph structure and estimate the model parameters
      cg <- sparse.tscgm(
        data = as.longitudinal(Y_train),
        model = "ar2",
        penalty = "scad",
        optimality = "bic_mod",
        control = list(silent = TRUE)
      )
      B_hat <- t(cg$gamma)
      C_hat <- 0
      
      # Set the lag length
      k_est <- 2
      
      # Measure the elapsed time
      t_diff <- difftime(Sys.time(), t_start, units = "secs")
      
      # Calculate the number of directed edges in the estimated model
      nd <- sum(cg$gamma != 0)
      
      # Calculate the number of undirected edges in the estimated model
      nu <- sum(cg$theta - diag(d) != 0)
      
      B_amat_SCAD2 <- 1*(B_hat != 0)
      Q_amat_SCAD2 <- 1*(cg$theta != 0)
      
    } else if (estimation_method == "LASSO2") {
      
      t_start <- Sys.time()
      
      # Learn the graph structure and estimate the model parameters
      cg <- sparse.tscgm(
        data = as.longitudinal(Y_train),
        model = "ar2",
        penalty = "lasso",
        optimality = "bic_mod",
        control = list(silent = TRUE)
      )
      B_hat <- t(cg$gamma)
      C_hat <- 0
      
      # Set the lag length
      k_est <- 2
      
      # Measure the elapsed time
      t_diff <- difftime(Sys.time(), t_start, units = "secs")
      
      # Calculate the number of directed edges in the estimated model
      nd <- sum(cg$gamma != 0)
      
      # Calculate the number of undirected edges in the estimated model
      nu <- sum(cg$theta - diag(d) != 0)
      
      B_amat_LASSO2 <- 1*(B_hat != 0)
      Q_amat_LASSO2 <- 1*(cg$theta != 0)
      
    }
    
    # Make one-step predictions from the test data; for the first prediction
    # steps, use training data when necessary
    for (i in 1:n_test) {
      Y_pred[i,] <- 0
      for (lag_no in 1:k_est) {
        lag_matrix <- B_hat[,(1+(lag_no-1)*d):(lag_no*d)]
        if (i - lag_no < 1) {
          Y_pred[i,] <-
            Y_pred[i,] + lag_matrix %*% Y_train[n_train + i - lag_no,]
        } else {
          Y_pred[i,] <- Y_pred[i,] + lag_matrix %*% Y_test[i - lag_no,]
        }
      }
      Y_pred[i,] <- Y_pred[i,] + C_hat
    }
    
    cat("Done.", "\n")
    
    # Output the elapsed time
    print(t_diff)
    
    # Output the (estimated or set) lag length
    cat("k  ", k_est, "\n")
        
    # Output the number directed edges in the estimated model
    cat("nd ", nd, "\n")
    
    # Output the number undirected edges in the estimated model
    cat("nu ", nu, "\n")
    
    # Calculate and output the mean squared error over all the channels
    MSE <- mean((Y_test - Y_pred)^2)
    cat("MSE", MSE, "\n\n")
    
    # Write the results to the output file
    sink(output_pathname, append = TRUE)
    cat(
      formatC(as.double(t_diff), digits = 4, width = 12, format = "f"),
      formatC(k_est, width = 10, format = "d"),
      formatC(as.double(nd), width = 10, format = "d"),
      formatC(as.double(nu), width = 10, format = "d"),
      formatC(as.double(MSE), digits = 4, width = 12, format = "f"),
      sep = ","
    )
    if (estimation_method == estimation_methods[length(estimation_methods)]) {
      cat("\n")
    } else {
      cat(",")
    }
    sink()

  }
  
  # Write the estimated adjacency matrices to disk
  save(
    B_amat_PLVAR,
    Q_amat_PLVAR,
    B_amat_PLVAR2,
    Q_amat_PLVAR2,
    B_amat_SCAD2,
    Q_amat_SCAD2,
    B_amat_LASSO2,
    Q_amat_LASSO2,
    file = paste(substr(recording_filename, 7, 14), ".rda", sep="")
  )

}