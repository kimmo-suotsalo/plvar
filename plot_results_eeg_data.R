# plot_results_eeg_data.R
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

library(RColorBrewer)

# Read the results into a data frame
results <- read.csv("./chbmit_output.csv")

# Set the method names to be displayed in the plots
method_names = c(
  "LS",
  "PLVAR",
  "LS2",
  "PLVAR2",
  "LASSO2",
  "SCAD2"
)

# Set the method colors to be displayed in the plots
method_colors = rep(c(brewer.pal(n = 6, name = 'Dark2'),"grey","grey"))

# Make a plot with two subplots
par(mfrow = c(1,2))

# Draw the first subplot
arr <- cbind(
  results$nd_OLS, results$nd_PLVAR, results$nd_OLS2,
    results$nd_PLVAR2, results$nd_LASSO2, results$nd_SCAD2,
    NaN, NaN,
  results$nu_OLS, results$nu_PLVAR, results$nu_OLS2,
    results$nu_PLVAR2, results$nu_LASSO2, results$nu_SCAD2,
    NaN, NaN,
  results$MSE_OLS, results$MSE_PLVAR, results$MSE_OLS2,
    results$MSE_PLVAR2, results$MSE_LASSO2, results$MSE_SCAD2
)
boxplot(
  arr,
  col = method_colors,
  border = method_colors,
  ylab = "Value",
  xlab = "Metric",
  xaxt = "n",
  yaxt = "n",
  main = "Estimated VAR models",
  log = "y",
  ylim = c(10, 1000000)
)
axis(
  1,
  at = c(3.5, 11.5, 19.5),
  labels = c(expression("n"[t]), expression("n"[c]), "MSE")
)
axis(
  2,
  at = c(10, 100, 1000, 10000, 100000, 1000000),
  labels = c("10", "100", "1000", "10000", "100000", "1000000")
)

# Draw the second subplot
arr <- cbind(
  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
  results$t_OLS, results$t_PLVAR, results$t_OLS2,
    results$t_PLVAR2, results$t_LASSO2, results$t_SCAD2,
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN
)
boxplot(
  arr,
  col = method_colors,
  border = method_colors,
  ylab = "Seconds",
  xaxt = "n",
  yaxt = "n",
  main = "Computation times",
  log = "y",
  ylim = c(0.1, 10000)
)
legend(
  x = "bottomright",
  legend = method_names,
  col = method_colors,
  lty = 1,
  lwd = 10,
  cex = 1,
  bty = "n"
)
axis(
  2,
  at = c(0.1, 1, 10, 100, 1000, 10000),
  labels = c("0.1", "1", "10", "100", "1000", "10000")
)

# Load an example of the estimated adjacency matrices
load("./chb06_01.rda")

# Get the number of variables
d <- nrow(Q_amat_PLVAR)

# Make a plot with three subplots
par(
  mfrow = c(3,1),
  mar = c(2, 1, 3, 1)
)

# Plot the PLVAR2 matrices
image(
  col = gray(1:0),
  z = t(apply(cbind(Q_amat_PLVAR2, B_amat_PLVAR2), 2, rev)),
  main = "PLVAR2",
  xaxt = "n",
  yaxt = "n"
)
box()
abline(v = 0.330)
abline(v = 0.667)
axis(
  1,
  at = c(0.167, 0.5, 0.833),
  labels = c(expression(Omega), expression("A"[1]), expression("A"[2]))
)

# Plot the LASSO2 matrices
image(
  col = gray(1:0),
  z = t(apply(cbind(Q_amat_LASSO2, B_amat_LASSO2), 2, rev)),
  main = "LASSO2",
  xaxt = "n",
  yaxt = "n"
)
box()
abline(v = 0.330)
abline(v = 0.667)
axis(
  1,
  at = c(0.167, 0.5, 0.833),
  labels = c(expression(Omega), expression("A"[1]), expression("A"[2]))
)

# Plot the SCAD2 matrices
image(
  col = gray(1:0),
  z = t(apply(cbind(Q_amat_SCAD2, B_amat_SCAD2), 2, rev)),
  main = "SCAD2",
  xaxt = "n",
  yaxt = "n"
)
box()
abline(v = 0.330)
abline(v = 0.667)
axis(
  1,
  at = c(0.167, 0.5, 0.833),
  labels = c(expression(Omega), expression("A"[1]), expression("A"[2]))
)