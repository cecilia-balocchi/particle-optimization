source("scripts/plot_functions.R")
load("data/alphabetas.RData")
load("data/partitions.RData")

## plot examples of alphas and betas
alpha_plot <- cbind(alpha_1[,1], alpha_2[,1], alpha_3[,1])
beta_plot <- cbind(beta_1[,1], beta_2[,1], beta_3[,1])

# par("mar") = c(bottom, left, top, right)
# par("usr") = c(x1, x2, y1, y2) user coordinates of the plotting region
# layout.show(8)

png("figures/true_partitions_alphabeta.png", width = 9, height = 6, units = "in", res = 600)

layout.matrix <- matrix(c(1,4,2,5,3,6,7,8), nrow = 2, ncol = 4)
layout(mat = layout.matrix,
       heights = c(2,2), # Heights of the two rows
       widths = c(2,2,2,0.5)) # Widths of the two columns


par(mar = c(2,3,2,1))
plot_partition_grid(gamma_alpha, W, values = alpha_plot[,1], max_value = max(alpha_plot), min_value = min(alpha_plot), title = "High separation - alpha", colorscheme = "PRGn", title_cex = 1.5)
plot_partition_grid(gamma_alpha, W, values = alpha_plot[,2], max_value = max(alpha_plot), min_value = min(alpha_plot), title = "Moderate separation - alpha", colorscheme = "PRGn", title_cex = 1.5)
plot_partition_grid(gamma_alpha, W, values = alpha_plot[,3], max_value = max(alpha_plot), min_value = min(alpha_plot), title = "Low separation - alpha", colorscheme = "PRGn", title_cex = 1.5)
plot_partition_grid(gamma_beta, W, values = beta_plot[,1], max_value = max(beta_plot), min_value = min(beta_plot), title = "High separation - beta", title_cex = 1.5)
plot_partition_grid(gamma_beta, W, values = beta_plot[,2], max_value = max(beta_plot), min_value = min(beta_plot), title = "Moderate separation - beta", title_cex = 1.5)
plot_partition_grid(gamma_beta, W, values = beta_plot[,3], max_value = max(beta_plot), min_value = min(beta_plot), title = "Low separation - beta", title_cex = 1.5)

par(mar = c(3,1,3,1), bty = "n")

col_list <- rev(brewer.pal(n = 5, name = "PRGn"))
plot(1, type = "n", xlim = c(0, 0.2), ylim = c(0,1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend_seq <- seq(par("usr")[3], par("usr")[4], length = 500)
for(leg_ix in 1:499){
  rect(par("usr")[1]/3,legend_seq[leg_ix], par("usr")[2]/3, legend_seq[leg_ix+1],
       border = NA, col = rgb(colorRamp(col_list, bias = 1)((leg_ix-1)/500)/255))
}
text(x =par("usr")[1]+0.15, y = seq(par("usr")[3]+0.015, par("usr")[4]-0.015, length = 9), cex = 1.2, 
     labels = round(seq(min(alpha_plot), max(alpha_plot), length = 9), digits = 2), xpd = TRUE)
rect(par("usr")[1]/3, par("usr")[3], par("usr")[2]/3, par("usr")[4])


col_list <- rev(brewer.pal(n = 5, name = "RdBu"))
plot(1, type = "n", xlim = c(0, 0.2), ylim = c(0,1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend_seq <- seq(par("usr")[3], par("usr")[4], length = 500)
for(leg_ix in 1:499){
  rect(par("usr")[1]/3,legend_seq[leg_ix], par("usr")[2]/3, legend_seq[leg_ix+1],
       border = NA, col = rgb(colorRamp(col_list, bias = 1)((leg_ix-1)/500)/255))
}
text(x =par("usr")[1]+0.15, y = seq(par("usr")[3]+0.015, par("usr")[4]-0.015, length = 9), cex = 1.2, 
     labels = round(seq(min(beta_plot), max(beta_plot), length = 9), digits = 2), xpd = TRUE)
rect(par("usr")[1]/3, par("usr")[3], par("usr")[2]/3, par("usr")[4])

dev.off()



