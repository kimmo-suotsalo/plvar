# plot_results_synthetic_data.R
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

library(RColorBrewer)

# Read in results

k <- 2
ds <- c(20,40,80)
m <- vector("list",3)

for (i in 1:length(m)){
  d <- ds[i]
  res <- readRDS(paste("results/sim_results_d", d, ".rds", sep = ""))
  tmp <- vector("list",6)
  tmp[[4]] <- t(sapply(res, function(x) rbind(x$Q_rec[,,k],array(NaN,c(2,5)))))
  tmp[[2]] <- t(sapply(res, function(x) rbind(x$Q_prec[,,k],array(NaN,c(2,5)))))
  tmp[[3]] <- t(sapply(res, function(x) rbind(x$B_rec[,,k],array(NaN,c(2,5)))))
  tmp[[1]] <- t(sapply(res, function(x) rbind(x$B_prec[,,k],array(NaN,c(2,5)))))
  tmp[[5]] <- t(sapply(res, function(x) rbind(x$run_time[,,k],array(NaN,c(2,5)))))
  tmp[[6]] <- t(sapply(res, function(x) rbind(x$k_estimated[1,,k])))
  m[[i]] <- tmp
}

title <- c("Temporal","Contemporaneous")

# Precision
pdf(paste("simulation_precision.pdf", sep = ""), width = 8, height = 4.2)
par(mfrow = c(2,3), mai = c(0.3, 0.5, 0.3, 0.2))

for (j in 1:2){
  for (i in 1:3){
    boxplot(m[[i]][[j]][,1:23], 
            #col = rep(c("steelblue","darkred","purple","grey","grey")), 
            col = rep(c(brewer.pal(n = 3, name = 'Dark2'),"grey","grey")), 
            border = rep(c(brewer.pal(n = 3, name = 'Dark2'),"grey","grey")),
            #xlab = "Sample size", 
            xaxt = "n", 
            yaxt = "n",
            frame = FALSE,
            #main = title[j,i],
            #col.main = "gray25",
            ylim = c(0, 1),
            notch = FALSE)
    box(col="gray25",lwd = 1)
    if (j == 1 & i == 1){
      legend(x = "bottomright", legend = c("PLVAR", "SCAD", "LASSO"), col = brewer.pal(n = 3, name = 'Dark2'), lty = 1, lwd = 5, cex = 1, bty = "n") 
    }    
    axis(1, at = c(2,7,12,17,22), labels = c("50","100","200","400","800"), col="gray25", col.ticks="gray25", col.axis="gray25")
    axis(2, at = c(0,0.5,1), labels = c("0.0","0.5","1.0"), col="gray25", col.ticks="gray25", col.axis="gray25")
    if (j == 1){
      mtext(text = paste("d = ",ds[i], sep = ""),  side = 3, line = 1, cex = 0.8)
    }
    if (i == 1){
      mtext(text = title[j], side = 2, line = 2.7, cex = 0.8)
    }
  }
}
dev.off()  

# Recall
pdf(paste("simulation_recall.pdf", sep = ""), width = 8, height = 4.2)
par(mfrow = c(2,3), mai = c(0.3, 0.5, 0.3, 0.2))

for (j in 1:2){
  for (i in 1:3){
    boxplot(m[[i]][[j+2]][,1:23], 
            #col = rep(c("steelblue","darkred","purple","grey","grey")), 
            col = rep(c(brewer.pal(n = 3, name = 'Dark2'),"grey","grey")), 
            border = rep(c(brewer.pal(n = 3, name = 'Dark2'),"grey","grey")),
            #xlab = "Sample size", 
            xaxt = "n", 
            yaxt = "n",
            frame = FALSE,
            #main = title[j,i],
            #col.main = "gray25",
            ylim = c(0, 1),
            notch = FALSE)
    box(col="gray25",lwd = 0.5)
    if (j == 1 & i == 1){
      legend(x = "bottomright", legend = c("PLVAR", "SCAD", "LASSO"), col = brewer.pal(n = 3, name = 'Dark2'), lty = 1, lwd = 5, cex = 1, bty = "n") 
    }
    axis(1, at = c(2,7,12,17,22), labels = c("50","100","200","400","800"), col="gray25", col.ticks="gray25", col.axis="gray25")
    axis(2, at = c(0,0.5,1), labels = c("0.0","0.5","1.0"), col="gray25", col.ticks="gray25", col.axis="gray25")
    if (j == 1){
      mtext(text = paste("d = ",ds[i], sep = ""),  side = 3, line = 1, cex = 0.8)
    }
    if (i == 1){
      mtext(text = title[j], side = 2, line = 2.7, cex = 0.8)
    }
  }
}
dev.off() 

# Running time
pdf(paste("simulation_time.pdf", sep = ""), width = 8, height = 2.1)
par(mfrow = c(1,3), mai = c(0.3, 0.5, 0.3, 0.2))

title <- c("Complete network, d = 20", "Complete network, d = 40", "Complete network, d = 80")
s <- 10^(0:4)
for (i in 1:3){
  boxplot(m[[i]][[5]][,1:23],
          log = "y",
          col.main = "gray25",
          #col = rep(c("  steelblue","darkred","purple","grey","grey")), 
          col = rep(c(brewer.pal(n = 3, name = 'Dark2'),"grey","grey")), 
          border = rep(c(brewer.pal(n = 3, name = 'Dark2'),"grey","grey")),
          #ylab = "net", 
          #xlab = "Sample size", 
          xaxt = "n",
          yaxt = "n",
          frame = FALSE,
          #main = title[i],
          ylim = c(0.1, s[2+i]*2),
          notch = FALSE)
  box(col="gray25", lwd = 1)
  if (i == 1){
    legend(x = "right", legend = c("PLVAR", "SCAD", "LASSO"), col = brewer.pal(n = 3, name = 'Dark2'), lty = 1, lwd = 5, cex = 1, bty = "n") 
  }
  axis(1, at = c(2,7,12,17,22), labels = c("50","100","200","400","800"), col="gray25", col.ticks="gray25", col.axis="gray25")
  axis(2, at = s[1:(i+2)], col="gray25", col.ticks="gray25", col.axis="gray25")
  mtext(text = paste("d = ",ds[i], sep = ""), side = 3, line = 1, cex = 0.8)
  if (i == 1){
    mtext(text = "Complete network", side = 2, line = 2.7, cex = 0.8)
  }
}
dev.off() 


library(pracma)


# Estimated lag lengths
mm <- vector("list",3)
for (j in 1:3){
  tmp <- matrix(0,5,5)
  rownames(tmp) <- 1:5
  colnames(tmp) <- c(50,100,200,400,800)
  for (i in 1:5){
    s <- table(m[[j]][[6]][,i])
    tmp[,i] <- histc(m[[j]][[6]][,i],1:5)$cnt/sum(s)
  }
  mm[[j]] <- tmp
}

pdf(paste("simulation_lag_length.pdf", sep = ""), width = 8, height = 2.1)
par(mfrow = c(1,3), mai = c(0.3, 0.5, 0.3, 0.2))

for (i in 1:3){
  barplot(t(mm[[i]]),
          xlab="Estimated lag length", 
          col=brewer.pal(6, "Blues")[2:6],
          #legend = paste("n =",colnames(mm[[i]])), 
          beside=TRUE)
  mtext(text = paste("d = ",ds[i], sep = ""), side = 3, line = 1, cex = 0.8)
  if (i == 1){
    legend(x = "topright", legend = paste("n =",colnames(mm[[i]])), col = brewer.pal(n = 6, name = 'Blues')[2:6], lty = 1, lwd = 5, cex = 1, bty = "n") 
  }
  
}
dev.off()