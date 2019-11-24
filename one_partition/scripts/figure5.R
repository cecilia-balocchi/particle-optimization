## Code to produce Figure 5 in manuscript ##

methods <- c("lam_1", "lam_10", "lam_100", "fixed_0", "fixed_1", "fixed_n", "km", "sc", "map")
estimation_perf <- matrix(nrow = length(methods), ncol = 6, dimnames = list(methods, c()))
selection_perf <- matrix(nrow = length(methods), ncol = 6, dimnames = list(methods, c()))

elbow_method <- function(x){
  x <- x/x[1]
  possible_k <- 3:9
  diff1 <- rep(NA, times = 7)
  diff2 <- rep(NA, times = 7)
  
  for(k_ix in 1:7){
    diff1[k_ix] <- x[possible_k[k_ix]] - x[possible_k[k_ix]-1]
    diff2[k_ix] <- x[possible_k[k_ix]+1] - 2*x[possible_k[k_ix]] + x[possible_k[k_ix]-1]
  }
  kappa <- abs(diff2)/(1 + diff1^2)^(3/2)
  return(possible_k[which.max(kappa)])
}

for(sim_number in 1:6){
  load(paste0("results/sim", sim_number, "_results.RData"))
  estimation_perf["lam_1",sim_number] <- mean(get(paste0("sim", sim_number, "_lam_1_results"))[,"RMSE"])
  selection_perf["lam_1", sim_number] <- mean(get(paste0("sim", sim_number, "_lam_1_results"))[,"RAND_top"])
  
  estimation_perf["lam_10",sim_number] <- mean(get(paste0("sim", sim_number, "_lam_10_results"))[,"RMSE"])
  selection_perf["lam_10", sim_number] <- mean(get(paste0("sim", sim_number, "_lam_10_results"))[,"RAND_top"])
  
  estimation_perf["lam_100",sim_number] <- mean(get(paste0("sim", sim_number, "_lam_100_results"))[,"RMSE"])
  selection_perf["lam_100", sim_number] <- mean(get(paste0("sim", sim_number, "_lam_100_results"))[,"RAND_top"])
  
  estimation_perf["fixed_0", sim_number] <- mean(get(paste0("sim", sim_number, "_fixed_0_results"))[,"RMSE"])
  selection_perf["fixed_0", sim_number] <- mean(get(paste0("sim", sim_number, "_fixed_0_results"))[,"RAND"])
  
  estimation_perf["fixed_1", sim_number] <- mean(get(paste0("sim", sim_number, "_fixed_1_results"))[,"RMSE"])
  selection_perf["fixed_1", sim_number] <- mean(get(paste0("sim", sim_number, "_fixed_1_results"))[,"RAND"])
  
  estimation_perf["fixed_n", sim_number] <- mean(get(paste0("sim", sim_number, "_fixed_n_results"))[,"RMSE"])
  selection_perf["fixed_n", sim_number] <- mean(get(paste0("sim", sim_number, "_fixed_n_results"))[,"RAND"])
  
  ### for k-means use the "elbow" method to get optimal k
  km_top_k <- apply(get(paste0("sim", sim_number, "_km_ss")), FUN = elbow_method, MARGIN = 1)
  tmp_rmse <- rep(NA, times = 20)
  tmp_rand <- rep(NA, times = 20)
  for(r in 1:20){
    tmp_rmse[r] <- get(paste0("sim", sim_number, "_km_results"))[r,paste0("RMSE_", km_top_k[r])]
    tmp_rand[r] <- get(paste0("sim", sim_number, "_km_results"))[r,paste0("RAND_", km_top_k[r])]
  }
  estimation_perf["km", sim_number] <- mean(tmp_rmse)
  selection_perf["km", sim_number] <- mean(tmp_rand)
  
  
  ### for spectral clustering, just use minimum sum of squares
  sc_top_k <- apply(get(paste0("sim", sim_number, "_sc_ss")), FUN = which.min, MARGIN = 1)
  tmp_rmse <- rep(NA, times = 20)
  tmp_rand <- rep(NA, times = 20)
  for(r in 1:20){
    tmp_rmse[r] <- get(paste0("sim", sim_number, "_sc_results"))[r, paste0("RMSE_", sc_top_k[r])]
    tmp_rand[r] <- get(paste0("sim", sim_number, "_sc_results"))[r, paste0("RAND_", sc_top_k[r])]
  }
  estimation_perf["sc", sim_number] <- mean(tmp_rmse)
  selection_perf["sc", sim_number] <- mean(tmp_rand)
  
  estimation_perf["map", sim_number] <- mean(get(paste0("sim", sim_number, "_map_results"))[,"RMSE"])
  selection_perf["map", sim_number] <- mean(get(paste0("sim", sim_number, "_map_results"))[,"RAND"])
}

png("figures/one_partition_performance.png", width = 8, height = 4, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex = 0.8, cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.8, mfrow = c(1,2))
plot(1, type = "n", ylim = c(0,1.2), xlim = c(0,6),
     main = "Estimation Performance", xlab = "Cluster Separation", ylab = "RMSE")
lines(5:0, estimation_perf["lam_1",], col = 'green')
lines(5:0, estimation_perf["sc",], col = 'red')
lines(5:0, estimation_perf["km",], col = 'blue')
lines(5:0, estimation_perf["fixed_n",], col = 'black')
lines(5:0, estimation_perf["fixed_1",], col = 'black')
points(5:0, estimation_perf["fixed_1",], col = 'black', pch = 5, cex = 0.8)
points(5:0, estimation_perf["fixed_n",], col = 'black', pch = 18, cex = 0.8)
points(5:0, estimation_perf["km",], col = 'blue', pch = 17)
points(5:0, estimation_perf["sc",], col = 'red', pch = 15)
points(5:0, estimation_perf["lam_1",], pch = 16, col = 'green', cex = 0.8)
legend("topleft", legend = c("Our procedure", "KM", "SC", "1-Cluster", "N-Clusters"), 
       pch = c(16, 17, 15, 5, 18), col = c("green", "blue", "red", "black", "black"),cex = 0.7, bty = "n")

plot(1, type = "n", ylim = c(0,1), xlim = c(0, 6), 
     main = "Partition Selection", xlab = "Cluster Separation", ylab = "Rand Index")
lines(5:0, selection_perf["lam_1",], col = 'green')
lines(5:0, selection_perf["sc",], col = 'red')
lines(5:0, selection_perf["km",], col = 'blue')
lines(5:0, selection_perf["fixed_n",], col = 'black')
lines(5:0, selection_perf["fixed_1",], col = 'black')

points(5:0, selection_perf["fixed_1",], col = 'black', pch = 5, cex = 0.8)
points(5:0, selection_perf["fixed_n",], col = 'black', pch = 18, cex = 0.8)
points(5:0, selection_perf["km",], col = 'blue', pch = 17)
points(5:0, selection_perf["sc",], col = 'red', pch = 15)
points(5:0, selection_perf["lam_1",], pch = 16, col = 'green', cex = 0.8)
legend("bottomleft", legend = c("Our procedure", "KM", "SC", "1-Cluster", "N-Clusters"), 
       pch = c(16, 17, 15, 5, 18), col = c("green", "blue", "red", "black", "black"),cex = 0.7, bty = "n")
dev.off()

