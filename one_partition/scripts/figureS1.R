### Additional plots ###
source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

load("data/partitions.RData")
load("data/alphas.RData")

load("results/illustration_sim1.RData")
load("results/illustration_sim3.RData")
load("results/illustration_sim5.RData")

png("figures/kmeans_sc_illustration.png", width = 6, height = 2.5, units = "in", res = 600)
lay_m2 <- matrix(nrow = 5, ncol = 12)
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

layout(lay_m2, widths = rep(1, times = 12), heights = rep(1, times = 5))

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
mtext(expression("True"~alpha~"'s (high separation)"), side = 3, line = 0, cex = 0.7)

plot(1, type = "n", xaxt = "n", yaxt = "n")
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_3[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (medium separation)"), side = 3, line = 0, cex = 0.7)

plot(1, type = "n", xaxt = "n", yaxt = "n")
par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_5[,1], max_value = max_value)
mtext(expression("True"~alpha~"'s (low separation)"), side = 3, line = 0, cex = 0.7)

plot(1, type = "n", xaxt = "n", yaxt = "n") # last blank 

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(km_sim1$particles[[3]], W, values = km_sim1$alpha[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("K-Means", side = 3, line = 0, cex = 0.7)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(sc_sim1$particles[[4]], W, values = sc_sim1$alpha[,4], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Spectral Clustering", side = 3, line = 0, cex = 0.7)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(km_sim3$particles[[3]], W, values = km_sim3$alpha[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("K-Means", side = 3, line = 0, cex = 0.7)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(sc_sim3$particles[[4]], W, values = sc_sim3$alpha[,4], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Spectral Clustering", side = 3, line = 0, cex = 0.7)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(km_sim5$particles[[3]], W, values = km_sim5$alpha[,3], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_1_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("K-Means", side = 3, line = 0, cex = 0.7)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(sc_sim5$particles[[4]], W, values = sc_sim5$alpha[,4], max_value = max_value)
#mtext(paste0("log-post = ", round(lam_100_sim1$log_post[1], digits = 2)), side = 3, line = 0, cex = 0.8)
mtext("Spectral Clustering", side = 3, line = 0, cex = 0.7)

dev.off()