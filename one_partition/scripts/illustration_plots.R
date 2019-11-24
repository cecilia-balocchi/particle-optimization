# Make mock-ups of the figures for Section 4 using data from the illustration
source("scripts/partition_functions.R")
source("scripts/plot_partition.R")
load("data/alphas.RData")
load("data/partitions.RData")

N <- 400
T <- 12

methods <- c("lam_1", "lam_10", "lam_100", "fixed_0", "fixed_1", "fixed_n", "km", "sc", "map")
estimation_perf <- matrix(nrow = length(methods), ncol = 6, dimnames = list(methods, c()))
selection_perf <- matrix(nrow = length(methods), ncol = 6, dimnames = list(methods, c()))


for(sim_number in 1:6){
  if(file.exists(paste0("results/illustration_sim", sim_number, ".RData"))){
    load(paste0("results/illustration_sim", sim_number, ".RData"))
    lam_1_results <- get(paste0("lam_1_results_sim", sim_number))
    lam_10_results <- get(paste0("lam_10_results_sim", sim_number))
    lam_100_results <- get(paste0("lam_100_results_sim", sim_number))
    fixed_0_results <- get(paste0("fixed_0_results_sim", sim_number))
    fixed_1_results <- get(paste0("fixed_1_results_sim", sim_number))
    fixed_n_results <- get(paste0("fixed_n_results_sim", sim_number))
    map_results <- get(paste0("map_results_sim", sim_number))
    km_results <- get(paste0("km_results_sim", sim_number))
    sc_results <- get(paste0("sc_results_sim", sim_number))
    

    estimation_perf["lam_1", sim_number] <- mean(lam_1_results[,"RMSE"], na.rm  = TRUE)
    estimation_perf["lam_10", sim_number] <- mean(lam_10_results[,"RMSE"], na.rm = TRUE)
    estimation_perf["lam_100", sim_number] <- mean(lam_100_results[,"RMSE"], na.rm = TRUE)
    estimation_perf["fixed_0", sim_number] <- mean(fixed_0_results[,"RMSE"], na.rm = TRUE)
    estimation_perf["fixed_1", sim_number] <- mean(fixed_1_results[,"RMSE"], na.rm = TRUE)
    estimation_perf["fixed_n", sim_number] <- mean(fixed_n_results[,"RMSE"], na.rm = TRUE)
    estimation_perf["map", sim_number] <- mean(map_results[,"RMSE"], na.rm = TRUE)
    estimation_perf["km", sim_number] <- mean(km_results[,"RMSE_top"], na.rm = TRUE)
    estimation_perf["sc", sim_number] <- mean(sc_results[,"RMSE_top"], na.rm = TRUE)
    
    
    selection_perf["lam_1", sim_number] <- mean(lam_1_results[,"RAND_top"], na.rm  = TRUE)
    selection_perf["lam_10", sim_number] <- mean(lam_10_results[,"RAND_top"], na.rm = TRUE)
    selection_perf["lam_100", sim_number] <- mean(lam_100_results[,"RAND_top"], na.rm = TRUE)
    selection_perf["fixed_0", sim_number] <- mean(fixed_0_results[,"RAND"], na.rm = TRUE)
    selection_perf["fixed_1", sim_number] <- mean(fixed_1_results[,"RAND"], na.rm = TRUE)
    selection_perf["fixed_n", sim_number] <- mean(fixed_n_results[,"RAND"], na.rm = TRUE)
    selection_perf["map", sim_number] <- mean(map_results[,"RAND"], na.rm = TRUE)
    selection_perf["km", sim_number] <- mean(km_results[,"RAND_top"], na.rm = TRUE)
    selection_perf["sc", sim_number] <- mean(sc_results[,"RAND_top"], na.rm = TRUE)
    
  }
  
  
}

# For the purposes of illustration
selection_perf["fixed_0",] <- rand_index(N, gamma_0, gamma_0)
selection_perf["fixed_1",] <- rand_index(N, gamma_1, gamma_0)
selection_perf["fixed_n",] <- rand_index(N, gamma_n, gamma_0)

png("figures/estimation_performance.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex = 0.8, cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.8)
plot(1, type = "n", ylim = c(0, 1.2), xlim = c(0, 6), 
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
legend("topleft", legend = c(expression("Our procedure"~lambda==1), "KM", "SC", "1-Cluster", "N-Clusters"), 
       pch = c(16, 17, 15, 5, 18), col = c("green", "blue", "red", "black", "black"),cex = 0.8, bty = "n")
dev.off()

png("figures/selection_performance.png", width = 4.5, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), cex = 0.8, cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.8)
plot(1, type = "n", ylim = c(0,1), xlim = c(0, 6), 
     main = "Partition Recovery", xlab = "Cluster Separation", ylab = "Rand Index")
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
legend("bottomleft", legend = c(expression("Our procedure"~lambda==1), "KM", "SC", "1-Cluster", "N-Clusters"), 
       pch = c(16, 17, 15, 5, 18), col = c("green", "blue", "red", "black", "black"),cex = 0.8, bty = "n")
dev.off()

###############
# Visualize the particle set for the different deltas
# Focus on Delta = 5 (sim_number = 1), Delta = 3 (sim_number = 3), Delta = 1 (sim_number = 5)


png("figures/top_particles.png", width = 6.5, height = 13/10 * 6.5, units = "in", res = 600)

lay_m <- matrix(nrow = 13, ncol = 8)
lay_m[1,] <- 1
lay_m[2:3, 1:2] <- 2 # true alphas sim1
lay_m[2:3, 3:4] <- 3 # top particle lam_1sim1
lay_m[2:3, 5:6] <- 4 # 2nd best particle lam_1sim1
lay_m[2:3, 7:8] <- 5 # 3rd best particle lam1_sim1
lay_m[4:5, 1:2] <- 6 # true alphas sim1
lay_m[4:5, 3:4] <- 7 # top particle lam100_sim1
lay_m[4:5, 5:6] <- 8 # 2nd best particle lam100_sim1
lay_m[4:5, 7:8] <- 9 # 3rd best particle lam100_sim1
lay_m[6:7, 1:2] <- 10 # true alpha sim3
lay_m[6:7, 3:4] <- 11 # top particle lam_1_sim3
lay_m[6:7, 5:6] <- 12 # 2nd best particle lam_1_sim3
lay_m[6:7, 7:8] <- 13 # 3rd best particle lam_1_sim3
lay_m[8:9, 1:2] <- 14 # true alphas sim3
lay_m[8:9, 3:4] <- 15 # top particle lam_100_sim3
lay_m[8:9, 5:6] <- 16 # 2nd best particle lam_100_sim3
lay_m[8:9, 7:8] <- 17 # 3rd best particle lam_100_sim3
lay_m[10:11, 1:2] <- 18 # true alphas sim5
lay_m[10:11,3:4] <- 19 # top particle lam_1_sim5
lay_m[10:11, 5:6] <- 20 # 2nd best particle lam_1_sim5
lay_m[10:11, 7:8] <- 21 # 3rd best particle lam_1_sim5
lay_m[12:13, 1:2] <- 22 # true alphs sim5
lay_m[12:13, 3:4] <- 23 # top particle lam_100_sim5
lay_m[12:13, 5:6] <- 24 # 2nd best particle lam_100_sim5
lay_m[12:13, 7:8] <- 25 # 3rd best particle lam_100_sim5


# legend
layout(lay_m, widths = rep(1, times = 8), heights = c(0.8, rep(1, times = 12)))


par(mar = c(0.75,0.5,2,0.5), bty = "n")
plot(1, type = "n", xlim = c(0, 1), ylim = c(0,3),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend_seq <- seq(par("usr")[1], par("usr")[2], length = 500)
for(leg_ix in 1:499){
  rect(legend_seq[leg_ix], par("usr")[3], legend_seq[leg_ix+1], par("usr")[4],
       border = NA, col = rgb(colorRamp(col_list, bias = 1)((leg_ix-1)/500)/255))
}
text(x = seq(par("usr")[1]+0.015, par("usr")[2]-0.015, length = 9), y = par("usr")[4]+1.5, cex = 1, 
     labels = round(seq(-max_value, max_value, length = 9), digits = 2), xpd = TRUE)

rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])

# true alphas sim 1
par(mar = c(0.75, 2, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_1[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s"), side = 3, line = 0, cex = 0.8)

# top partition
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value) 
mtext(paste0("log-post = ",round(lam_1_sim1$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)

# 2nd best partition
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
mtext(paste0("log-post = ",round(lam_1_sim1$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)


# 3rd best partition
par(mar = c(0.75, 0.5, 0.75, 2))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
mtext(paste0("log-post = ",round(lam_1_sim1$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)


# true alphas sim 1
par(mar = c(0.75, 2, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_1[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s"), side = 3, line = 0, cex = 0.8)
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[1]], W, values = lam_100_sim1$alpha_hat_particle[,1], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim1$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[2]], W, values = lam_100_sim1$alpha_hat_particle[,2], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim1$log_post[2], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 2.0))
plot_partition_grid(lam_100_sim1$particles[[3]], W, values = lam_100_sim1$alpha_hat_particle[,3], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim1$log_post[3], digits = 3)), side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 2, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_3[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s"), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[1]], W, values = lam_1_sim3$alpha_hat_particle[,1], max_value = max_value)
mtext(paste0("log-post = ", round(lam_1_sim3$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[2]], W, values = lam_1_sim3$alpha_hat_particle[,2], max_value = max_value)
mtext(paste0("log-post = ", round(lam_1_sim3$log_post[2], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 2))
plot_partition_grid(lam_1_sim3$particles[[3]], W, values = lam_1_sim3$alpha_hat_particle[,3], max_value = max_value)
mtext(paste0("log-post = ", round(lam_1_sim3$log_post[3], digits = 3)), side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 2, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_3[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s"), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[1]], W, values = lam_100_sim3$alpha_hat_particle[,1], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim3$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[2]], W, values = lam_100_sim3$alpha_hat_particle[,2], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim3$log_post[2], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 2))
plot_partition_grid(lam_100_sim3$particles[[3]], W, values = lam_100_sim3$alpha_hat_particle[,3], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim3$log_post[3], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 2, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_5[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s"), side = 3, line = 0, cex = 0.8)
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[1]], W, values = lam_1_sim5$alpha_hat_particle[,1], max_value = max_value)
mtext(paste0("log-post = ", round(lam_1_sim5$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[2]], W, values = lam_1_sim5$alpha_hat_particle[,2], max_value = max_value)
mtext(paste0("log-post = ", round(lam_1_sim5$log_post[2], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 2))
plot_partition_grid(lam_1_sim5$particles[[3]], W, values = lam_1_sim5$alpha_hat_particle[,3], max_value = max_value)
mtext(paste0("log-post = ", round(lam_1_sim5$log_post[2], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 2, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_5[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s"), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[1]], W, values = lam_100_sim5$alpha_hat_particle[,1], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim5$log_post[1], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[2]], W, values = lam_100_sim5$alpha_hat_particle[,2], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim5$log_post[2], digits = 3)), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 2))
plot_partition_grid(lam_100_sim5$particles[[3]], W, values = lam_100_sim5$alpha_hat_particle[,3], max_value = max_value)
mtext(paste0("log-post = ", round(lam_100_sim5$log_post[3], digits = 3)), side = 3, line = 0, cex = 0.8)

dev.off()



## An alternative plot


png("figures/top_partitions2.png", width = 6, height = 4.5, units = "in", res = 600)
lay_m2 <- matrix(nrow = 9, ncol = 12)
lay_m2[1,] <- 1
lay_m2[2:3,1] <- 2 # blank
lay_m2[2:3, 2:3] <- 3 # true alphas sim 1
lay_m2[2:3,4:5] <- 4 # blank
lay_m2[2:3,6:7] <- 5 # true alphas sim 2
lay_m2[2:3, 8:9] <- 6 # blank
lay_m2[2:3, 10:11] <- 7 # true alphs sim 3
lay_m2[2:3, 12] <- 8 # blank

lay_m2[4:5, 1:2] <- 9 # particle 1 lam_1_sim1
lay_m2[4:5, 3:4] <- 10 # particle 1 lam_100_sim1
lay_m2[4:5, 5:6] <- 11 # particle 1 lam_1_sim3
lay_m2[4:5, 7:8] <- 12 # particle 1 lam_100_sim3
lay_m2[4:5, 9:10] <- 13 # particle 1 lam_1_sim5
lay_m2[4:5, 11:12] <- 14 # particle 1 lam_100_sim5

lay_m2[6:7, 1:2] <- 15
lay_m2[6:7, 3:4] <- 16
lay_m2[6:7, 5:6] <- 17
lay_m2[6:7, 7:8] <- 18
lay_m2[6:7, 9:10] <- 19
lay_m2[6:7, 11:12] <- 20

lay_m2[8:9, 1:2] <- 21
lay_m2[8:9, 3:4] <- 22
lay_m2[8:9, 5:6] <- 23
lay_m2[8:9, 7:8] <- 24
lay_m2[8:9, 9:10] <- 25
lay_m2[8:9, 11:12] <- 26

layout(lay_m2, widths = rep(1, times = 12), heights = rep(1, times = 9))

par(mar = c(0.75,0.5,2,0.5), bty = "n")
plot(1, type = "n", xlim = c(0, 1), ylim = c(0,3),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend_seq <- seq(par("usr")[1], par("usr")[2], length = 500)
for(leg_ix in 1:499){
  rect(legend_seq[leg_ix], par("usr")[3], legend_seq[leg_ix+1], par("usr")[4],
       border = NA, col = rgb(colorRamp(col_list, bias = 1)((leg_ix-1)/500)/255))
}
text(x = seq(par("usr")[1]+0.015, par("usr")[2]-0.015, length = 9), y = par("usr")[4]+1.5, cex = 1, 
     labels = round(seq(-max_value, max_value, length = 9), digits = 2), xpd = TRUE)

rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])

plot(1, type = "n", xaxt = "n", yaxt = "n")
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_1[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (high separation)"), side = 3, line = 0, cex = 0.8)

plot(1, type = "n", xaxt = "n", yaxt = "n")
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_3[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (medium separation)"), side = 3, line = 0, cex = 0.8)

plot(1, type = "n", xaxt = "n", yaxt = "n")
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_5[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (low separation)"), side = 3, line = 0, cex = 0.8)

plot(1, type = "n", xaxt = "n", yaxt = "n") # last blank 

# Plot Top particles
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 1", side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[1]], W, values = lam_100_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 1", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[1]], W, values = lam_1_sim3$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim3$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 1", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[1]], W, values = lam_100_sim3$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim3$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 1", side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[1]], W, values = lam_1_sim5$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim5$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 1", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[1]], W, values = lam_100_sim5$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim5$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 1", side = 3, line = 0, cex = 0.8)


# Plot 2nd particle
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 2", side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[2]], W, values = lam_100_sim1$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[2]], W, values = lam_1_sim3$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim3$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[2]], W, values = lam_100_sim3$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim3$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[2]], W, values = lam_1_sim5$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim5$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[2]], W, values = lam_100_sim5$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim5$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 2", side = 3, line = 0, cex = 0.8)

# Plot Particle 3
par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==1), side = 1, line = 0, cex = 0.8)

par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[3]], W, values = lam_100_sim1$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==100), side = 1, line = 0, cex = 0.8)


par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[3]], W, values = lam_1_sim3$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim3$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==1), side = 1, line = 0, cex = 0.8)

par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[3]], W, values = lam_100_sim3$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim3$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==100), side = 1, line = 0, cex = 0.8)


par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[3]], W, values = lam_1_sim5$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim5$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==1), side = 1, line = 0, cex = 0.8)

par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[3]], W, values = lam_100_sim5$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim5$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==100), side = 1, line = 0, cex = 0.8)

dev.off()
