tmp = load("data/alphabetas.RData")
tmp = load("data/partitions.RData")
source("scripts/plot_functions.R")

sim_number <- 3 #1
batch <- 1
batch_size <- 4

if(sim_number == 1){
  sep_str <- "hs"
} else if(sim_number == 2){
  sep_str <- "ms"
}else if(sim_number == 3){
  sep_str <- "ls"
}

ALPHA <- get(paste0("alpha_", sim_number))
BETA <- get(paste0("beta_", sim_number))
r = 2
alpha <- ALPHA[, r + (batch-1)*batch_size]
beta <- BETA[, r + (batch-1)*batch_size]

# max separation parameters
hs_alpha <- alpha_1[, r + (batch-1)*batch_size]
hs_beta <- beta_1[, r + (batch-1)*batch_size]

tmp = load(paste0("results/counts_newX_sim", sim_number, "_", batch, ".RData"))

lam_output <- map
j_max <- which.max(lam_output$w)
gammaA_top <- lam_output$particle_set_A[[j_max]]
gammaB_top <- lam_output$particle_set_B[[j_max]]

### PLOT TRUE FIRST
plot_partition_grid(gamma_alpha, W, alpha,
                    max_value = max(c(alpha, lam_output$alpha_particle[,j_max], hs_alpha)),
                    min_value = min(c(alpha, lam_output$alpha_particle[,j_max], hs_alpha)))
plot_partition_grid(gamma_beta, W, beta,
                    max_value = max(c(beta, lam_output$beta_particle[,j_max], hs_beta)),
                    min_value = min(c(beta, lam_output$beta_particle[,j_max], hs_beta)))

### PLOT MAP
plot_partition_grid(gammaA_top, W, lam_output$alpha_particle[,j_max], 
                    max_value = max(c(alpha, lam_output$alpha_particle[,j_max], hs_alpha)),
                    min_value = min(c(alpha, lam_output$alpha_particle[,j_max], hs_alpha)))
plot_partition_grid(gammaB_top, W, lam_output$beta_particle[,j_max], 
                    max_value = max(c(beta, lam_output$beta_particle[,j_max], hs_beta)),
                    min_value = min(c(beta, lam_output$beta_particle[,j_max], hs_beta)))

## note:lam_1 has several overlapping: 
## (both for sim_number = 1, 2)
# lam_output <- lam_1$output
## this is specifically for sim_number = 1
# table(lam_output$w)
# 0.0536397696131343 0.0984073604233352  0.101232908271177 
# 1                  1                  1 
# 0.103362781142942  0.114953027988821 
# 5                  2
## so let's not plot those.
## lam_10 has no overlapping particles
## we will include these results in paper (as they are more interesting)

lam_output <- lam_10$output

png(paste0("figures/laplace_",sep_str,"_true_partitions_alphabeta.png"), width = 9, height = 4.5, units = "in", res = 600)
par(mfrow = c(1,2), mar = c(1.5,1,1.5,0.5))
plot_partition_grid(gamma_alpha, W, alpha,
                    max_value = max(c(alpha, lam_output$alpha_particle, hs_alpha)),
                    min_value = min(c(alpha, lam_output$alpha_particle, hs_alpha)),
                    title = paste("alpha"),title_cex = 1.5, 
                    colorscheme = "PRGn")
plot_partition_grid(gamma_beta, W, beta,
                    max_value = max(c(beta, lam_output$beta_particle, hs_beta)),
                    min_value = min(c(beta, lam_output$beta_particle, hs_beta)),
                    title = paste("beta"),title_cex = 1.5,)
dev.off()

layout.mat <- matrix(c(1,2,0,3,4,5,6,0,7,8,9,10,0,11,12,13,14,0,15,16,17,18,0,19,20), 
                     nrow = 5, ncol = 5, byrow = T)
png(paste0("figures/laplace_",sep_str,"_estimated_partitions_lam10.png"), width = 8, height = 10, units = "in", res = 600)
layout(mat = layout.mat,
       heights = c(1,1,1,1,1),    # Heights of the two rows
       widths = c(1,1,0.1,1,1)) 
for(j in order(lam_output$w, decreasing = T)){
  gammaA_top <- lam_output$particle_set_A[[j]]
  gammaB_top <- lam_output$particle_set_B[[j]]
  
  par(mar = c(1.3,0.2,1,0.2))
  plot_partition_grid(gammaA_top, W, lam_output$alpha_particle[,j],
                      max_value = max(c(alpha, lam_output$alpha_particle, hs_alpha)),
                      min_value = min(c(alpha, lam_output$alpha_particle, hs_alpha)), 
                      title = paste("alpha"), title_cex = 1.3, 
                      colorscheme = "PRGn")
  par(mar = c(1.3,0.2,1,0.2))
  plot_partition_grid(gammaB_top, W, lam_output$beta_particle[,j],
                      max_value = max(c(beta, lam_output$beta_particle, hs_beta)),
                      min_value = min(c(beta, lam_output$beta_particle, hs_beta)),
                      title = paste("beta - w =",round(lam_output$w[j],3)),
                      title_cex = 1.3)
}
dev.off()
