# Generate data
library(igraph)
library(MASS)
source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

n_side <- 20
N <- 400
g <- make_lattice(length = n_side, dim = 2)

W <- as_adj(g, type = "both")
W <- matrix(W, nrow = N)

quad_1 <- rep(11:20, times = 10) + 20*rep(10:19, each = 10)
quad_2 <- rep(1:10, times = 10) + 20*rep(10:19, each = 10)
quad_3 <- rep(1:10, times = 10) + 20*rep(0:9, each = 10)
quad_4 <- rep(11:20, times = 10) + 20*rep(0:9, each = 10)

# For the large grid, we need to define the actual clusters
cluster_1 <- c(165, 166, 184:187, 204:207, 225, 226)
cluster_2 <- c(quad_2, quad_3)
cluster_2 <- cluster_2[!cluster_2 %in% cluster_1]
cluster_3 <- quad_4
cluster_4 <- quad_1

gamma_0 <- list(cluster_1, cluster_2, cluster_3, cluster_4)
gamma_1 <- list(c(quad_1, quad_2, quad_3, quad_4))
gamma_2 <- list(c(quad_1, quad_2), c(quad_3, quad_4))
gamma_3 <- list(quad_1, quad_2, c(quad_3, quad_4))
gamma_4 <- list(quad_1, quad_2, quad_3, quad_4)
gamma_n <- as.list(1:N)


seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

sigma2 <- 1
rho <- 0.95
#a1 <- 1/10
a1 <- 1/5
tmp <- partition_config(gamma_0, W)
K <- tmp$K
config <- tmp$config

n_sim <- 100
delta <- c(5, 3, 2, 1, 0.75, 0)
for(sim_number in 1:6){
  
  alphabar <- c(0, -1, 1, 0) * delta[sim_number]
  alpha <- matrix(NA, ncol= n_sim, nrow= N)
  for(k in 1:K){
    n_k <- config[k]
    W_k <- W[gamma_0[[k]], gamma_0[[k]]]
    D <- diag(rowSums(W_k))
    W_star_k <- D - W_k
    Omega_alpha <- rho * W_star_k + (1 - rho) * diag(n_k)
    Sigma_alpha <- solve(Omega_alpha)
    set.seed(seed_seq[sim_number] + k)
    
    alpha[gamma_0[[k]],] <- t(matrix(mvrnorm(n = n_sim, mu = rep(alphabar[k], times = n_k), Sigma = a1 * sigma2 * Sigma_alpha), nrow = n_sim))
  }
  assign(paste0("alpha_", sim_number), alpha)
}

max_value <- ceiling(max(abs(c(alpha_1, alpha_2, alpha_3, alpha_4, alpha_5))))

save(gamma_0, gamma_1, gamma_2, gamma_3, gamma_4, gamma_n, W, file = "data/partitions.RData")
save(sigma2, rho, a1, alpha_1, alpha_2, alpha_3, alpha_4, alpha_5,
     max_value, file ="data/alphas.RData")
#######################################
# Plot the alpha's for each pair (alpha1+alpha5, etc.)

lay_m <- matrix(nrow = 5, ncol = 6)
lay_m[1,] <- 1
lay_m[2:3, 1:2] <- 2
lay_m[2:3, 3:4] <- 3
lay_m[2:3, 5:6] <- 4
lay_m[4:5, 1:2] <- 5
lay_m[4:5, 3:4] <- 6
lay_m[4:5, 5:6] <- 7


png("figures/alphas_full.png", width = 6, height = 5, units = "in", res = 600)
layout(lay_m, widths = rep(1, times = 6), heights = c(0.5, rep(1, times = 4)))
# legend
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

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_1[,1], max_value = max_value)
mtext(expression(Delta == 5), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_2[,40], max_value = max_value)
mtext(expression(Delta == 3), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_3[,21], max_value = max_value)
mtext(expression(Delta == 2), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_4[,76], max_value = max_value)
mtext(expression(Delta == 1), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_5[,1], max_value = max_value)
mtext(expression(Delta == 0.75), side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_6[,31], max_value = max_value)
mtext(expression(Delta == 0), side = 3, line = 0, cex = 0.8)

dev.off()


lay_m2 <- matrix(nrow = 3, ncol = 6)
lay_m2[1,] <- 1
lay_m2[2:3, 1:2] <- 2
lay_m2[2:3, 3:4] <- 3
lay_m2[2:3, 5:6] <- 4



png("figures/alphas_reduced.png", width = 6, height = 6/2.5, units = "in", res = 600)
layout(lay_m2, widths = rep(1, times = 6), heights = c(0.5, rep(1, times = 2)))
# legend
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

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_1[,1], max_value = max_value)
mtext("High separation", side = 3, line = 0, cex = 0.8)

par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_3[,21], max_value = max_value)
mtext(expression(Delta == 2), side = 3, line = 0, cex = 0.8)


par(mar = c(0.75, 0.5, 0.75, 0.5))
plot_partition_grid(gamma_0, W, values = alpha_6[,31], max_value = max_value)
mtext(expression(Delta == 0), side = 3, line = 0, cex = 0.8)


dev.off()





png("figures/new_alphas_hist.png", width = 4, height = 4, units = "in", res = 300)
layout(matrix(1:4, 2,2, byrow = TRUE), widths = rep(1, times = 2), heights = rep(1, times = 2))
# legend
par(mar = c(3,3,2,1), bty = "n")
hist(alpha_1[,1], xlim = c(-1.1, 1.1)*max_value, breaks = 100, main = expression(Delta == 5), ylim = c(0, 25))
#abline(v = c(-5,0, 5), lty = 2, col = 'red')

par(mar = c(3,3,2,1), bty = "n")
hist(alpha_5[,1], xlim = c(-1.1, 1.1)*max_value, breaks = 100, main = expression(Delta == 3), ylim = c(0, 25))
#abline(v = c(-2,0, 2), lty = 2, col = 'red')

par(mar = c(3,3,2,1), bty = "n")
hist(alpha_2[,1], xlim = c(-1.1, 1.1)*max_value, breaks = 100, main = expression(Delta == 2), ylim = c(0, 25))
#abline(v = c(-1,0, 1), lty = 2, col = 'red')

par(mar = c(3,3,2,1), bty = "n")
hist(alpha_3[,1], xlim = c(-1.1, 1.1)*max_value, breaks = 100, main = expression(Delta == 1), ylim = c(0, 25))
#abline(v = c(-0.75,0, 0.75), lty = 2, col = 'red')

dev.off()



