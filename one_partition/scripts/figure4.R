source("scripts/partition_functions.R")
source("scripts/plot_partition.R")
load("data/alphas.RData")
load("data/partitions.RData")

N <- 400
T <- 12

for(sim_number in 1:6){
  if(file.exists(paste0("results/illustration_sim", sim_number, ".RData"))){
    load(paste0("results/illustration_sim", sim_number, ".RData"))
  }
}

png("figures/top_partitions.png", width = 6, height = 4.5, units = "in", res = 600)
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

par(mar = c(1,0.5,2,0.5), bty = "n")
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
par(mar = c(0.5, 0.5, 1.25, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_1[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (high separation)"), side = 3, line = 0, cex = 0.8)

plot(1, type = "n", xaxt = "n", yaxt = "n")
par(mar = c(0.5, 0.5, 1.25, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_3[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (medium separation)"), side = 3, line = 0, cex = 0.8)

plot(1, type = "n", xaxt = "n", yaxt = "n")
par(mar = c(0.5, 0.5, 1.25, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_5[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (low separation)"), side = 3, line = 0, cex = 0.8)

plot(1, type = "n", xaxt = "n", yaxt = "n") # last blank 

# Plot Top particles
par(mar = c(0.75, 1.25, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 1", side = 3, line = 0, cex = 0.8)
mtext("Particle 1", side = 2, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[1]], W, values = lam_100_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 1", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[1]], W, values = lam_1_sim3$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim3$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 1", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[1]], W, values = lam_100_sim3$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim3$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 1", side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[1]], W, values = lam_1_sim5$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim5$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 1", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[1]], W, values = lam_100_sim5$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim5$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 1", side = 3, line = 0, cex = 0.8)


# Plot 2nd particle
par(mar = c(0.75, 1.25, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 2", side = 2, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[2]], W, values = lam_100_sim1$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[2]], W, values = lam_1_sim3$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim3$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[2]], W, values = lam_100_sim3$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim3$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[2]], W, values = lam_1_sim5$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim5$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 2", side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[2]], W, values = lam_100_sim5$alpha_hat_particle[,2], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim5$log_post[2], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 2", side = 3, line = 0, cex = 0.8)

# Plot Particle 3
par(mar = c(1, 1.25, 0.75, 0.5))
plot_partition_grid(lam_1_sim1$particles[[1]], W, values = lam_1_sim1$alpha_hat_particle[,1], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Particle 3", side = 2, line = 0, cex = 0.8)
mtext(expression(lambda==1), side = 1, line = 0, cex = 0.8)

par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim1$particles[[3]], W, values = lam_100_sim1$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==100), side = 1, line = 0, cex = 0.8)


par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim3$particles[[3]], W, values = lam_1_sim3$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim3$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==1), side = 1, line = 0, cex = 0.8)

par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim3$particles[[3]], W, values = lam_100_sim3$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim3$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==100), side = 1, line = 0, cex = 0.8)


par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_1_sim5$particles[[3]], W, values = lam_1_sim5$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim5$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==1), side = 1, line = 0, cex = 0.8)

par(mar = c(1, 0.5, 0.75, 0.5))
plot_partition_grid(lam_100_sim5$particles[[3]], W, values = lam_100_sim5$alpha_hat_particle[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim5$log_post[3], digits = 2)), side = 3, line = 0, cex = 0.8)
#mtext("Particle 3", side = 3, line = 0, cex = 0.8)
mtext(expression(lambda==100), side = 1, line = 0, cex = 0.8)

dev.off()
