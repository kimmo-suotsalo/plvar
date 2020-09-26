# run_example.R
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

# A short example illustrating how to use the code.

### 1. GENERATE MODEL AND DATA

# We use the "SparseTSCGM" package to generate example data.
library(SparseTSCGM)

# Generate AR(2) model: 10 variables, 200 observations, random network
d <- 10
k <- 2

set.seed(2)
mod <- sim.data(model = paste("ar", k, sep = ""),
                time = 200,
                n.obs = 2, 
                n.var = d,
                prob0 = (3*d)/(k*(d^2)), 
                network = "random")
data <- mod$data1[c(TRUE,FALSE),]
true_B <- t(mod$gamma)
true_Q <- mod$theta
true_B_amat <- (true_B != 0)*1
true_Q_amat <- (true_Q != 0)-diag(d)

# Check that generated model is stable, i.e. no divergence towards -/+ infinity
print(max(abs(data)))

# Proper way of checking stability
if (k == 1){
  A <- t(mod$gamma)
} else {
  A <- rbind(t(mod$gamma), diag(d*k)[1:(d*(k-1)),])
}

if (any(abs(eigen(A)$values)>=1)){
  stop("VAR model is not stable, try changing the seed!")
}


### 2. ESTIMATE STRUCTURE

# Load structure learning function
source("gvar_learn_structure.R")

# Learn structure, set K=5 as max lag length considered during search
gvar <- gvar_learn_structure(data, K = 5)

# Plot the lag length scores
plot(gvar$B_score)

# The estimated lag length is the one maximizing the above score, should be 2
k_est <- gvar$est_k

# Extract structures corresponding to the estimated lag length
B_amat <- gvar$B_amat[,1:(k_est*10),k_est] # B_amat contains K*d columns where everything after k*d are zero, therefore the first k_est*10 columns are selected
Q_amat <- gvar$Q_amat[,,k_est]

# Calculate Hamming distance to true structures
hd_temp <- sum(B_amat != true_B_amat)
hd_contemp <- sum(Q_amat != true_Q_amat)

### 3. ESTIMATE MODEL PARAMETERS

# Load parameter estimation function
source("gvar_estimate_parameters.R")

# Compute MLEs
gvar_mod <- gvar_estimate_parameters(data, B_amat, Q_amat)
B_hat <- gvar_mod$B_hat
Q_hat <- gvar_mod$Q_hat

# Compare true parameters and estimated parameters
print(B_hat)
print(true_B)
print(paste("MSE between true and estimated lag parameters:", mean((true_B-B_hat)^2)))

print(Q_hat)
print(true_Q)
print(paste("MSE between true and estimated precision matrices:", mean((true_Q-Q_hat)^2)))
