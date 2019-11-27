source("scripts/plot_functions.R")

# Make MLE plots
Y <- read.csv("data/Y_tracts.csv", header = F)
X <- read.csv("data/X_tracts.csv", header = F)
Y <- as.matrix(Y)
X <- as.matrix(X)
X <- X - rowMeans(X)
betas.mle <- numeric(n_tr)
for(i in 1:n_tr)
  betas.mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
alphas.mle <- rowMeans(Y) - log(2)
# note we need to subtract log(2) because the inverse  hyperbolic since transformation (ihst) function use
# did not have the shift factor "-log(2)". Instead of re-running the algorithm, we can shift the alpha found.

phillypoly <- tracts
polyfortified <- fortify(phillypoly)
p1 <- plot_borders_ggplot(rep(0, n_tr), phillypoly, w.sym, var = alphas.mle, polyfortified = polyfortified, legend = TRUE, 
                    legend_name = "Mean\nlevel", palette_num = 3, map = googlemap)
p2 <- plot_borders_ggplot(rep(0, n_tr), phillypoly, w.sym, var = betas.mle, polyfortified = polyfortified, legend = TRUE, 
                    legend_name = "Time\ntrend", palette_num = 5, map = googlemap)
ggsave(filename = "alpha_mle.png", plot = p1, device = "png", path = "figures/",
         width = 2, height = 2, units = "cm", scale = 10)
ggsave(filename = "beta_mle.png", plot = p2, device = "png", path = "figures/",
         width = 2, height = 2, units = "cm", scale = 10)
## Strings should give the paths of the 4 output files from the Particle Optimization algorithm
strings <- c("results/EPPrior5_0_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata",
             "results/UnifPrior_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata",
             "results/UnifEPPrior5_0_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata",
             "results/EPUnifPrior5_0_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata")

## Create the plots with the best three particles (Figure 6)
for(prior in 1:4){
  input_string <- strings[prior]
  if(prior == 1){
    plot_label <- "Ewens-Pitman prior"
    output_filename <-  "tract_ep"
  } else if(prior == 2){
    output_filename <-  "tract_unif"
    plot_label <- "Uniform prior"
  } else if(prior == 3){
    output_filename <-  "tract_epunif"
    plot_label <- "EP + Uniform prior"
  } else if(prior == 4){
    output_filename <-  "tract_unifep"
    plot_label <- "Uniform + EP prior"
  }
  plot_list <- plot_particles_diff(input_string, alpha_shift= -log(2))
  ggsave(filename = paste0(output_filename, ".png"),
         plot = arrangeGrob(grobs = plot_list, nrow = 2, widths=c(rep(5/27,5),2/27),
                            left=textGrob(plot_label, rot = 90, gp = gpar(fontsize = 22))),
         device = "png", path = "figures/", width = 5.75, height = 2.5, units = "cm", scale = 6)
}

## Create the plots that combine the differences between the top particle and the remaining  (Figure 7)
for(prior in 1:2){
  input_string <- strings[prior]
  if(prior == 1){
    plot_label <- "Ewens-Pitman prior"
    output_filename <-  "tract_ep"
  } else if(prior == 2){
    output_filename <-  "tract_unif"
    plot_label <- "Uniform prior"
  }
  ret <- combine_diff_one(strings = input_string, k_combined = 10, AorB = TRUE, map = googlemap)
  ggsave(filename = paste0(output_filename,"_newhyp2K6_L10_combinedA.png"), plot = ret[[1]],
         device = "png", path = "figures/", width = 2, height = 2, units = "cm", scale = 6)
  ret <- combine_diff_one(strings = input_string, k_combined = 10, AorB = FALSE, map = googlemap)
  ggsave(filename = paste0(output_filename,"_newhyp2K6_L10_combinedB.png"), plot = ret[[1]],
         device = "png", path = "figures/", width = 2, height = 2, units = "cm", scale = 6)
}

## Create the plot of the best particle for each model (Figure 8)
alphas <- array(NA, dim = c(n_tr, 4))
betas <- array(NA, dim = c(n_tr, 4))
z1s <- array(NA, dim = c(n_tr, 4))
z2s <- array(NA, dim = c(n_tr, 4))
tracts_fort <- fortify(tracts)
ps <- list()
for(prior in 1:4){
  input_string <- strings[prior]
  if(prior == 1){
    plot_label <- "Ewens-Pitman prior"
    output_filename <-  "tract_ep"
  } else if(prior == 2){
    output_filename <-  "tract_unif"
    plot_label <- "Uniform prior"
  } else if(prior == 3){
    output_filename <-  "tract_epunif"
    plot_label <- "EP + Uniform prior"
  } else if(prior == 4){
    output_filename <-  "tract_unifep"
    plot_label <- "Uniform + EP prior"
  }
  load(input_string)
  i = prior
  j <- which(tmp$w==max(tmp$w))
  partitionA <- tmp$particle_set_A[[j]]
  partitionB <- tmp$particle_set_B[[j]]
  alphas[,i] <- tmp$alpha_particle[,j] - log(2)
  betas[,i] <- tmp$beta_particle[,j]
  z1 <- numeric(n_tr)
  for(k in 1:length(partitionA)){
    z1[partitionA[[k]]] <- k
  }
  z2 <- numeric(n_tr)
  for(k in 1:length(partitionB)){
    z2[partitionB[[k]]] <- k
  }
  z1s[,i] <- z1
  z2s[,i] <- z2
}

limitsA <- c(min(alphas), max(alphas))
limitsB <- c(min(betas), max(betas))
for(prior in 1:4){ 
  i=prior
  if(prior == 1){
    output_filename <-  "tract_ep_newhyp2K6_bestpart.png"
  } else if(prior == 2){
    output_filename <-  "tract_unif_newhyp2K6_bestpart.png" 
  } else if(prior == 3){
    output_filename <-  "tract_epunif_newhyp2K6_bestpart.png"
  } else if(prior == 4){
    output_filename <-  "tract_unifep_newhyp2K6_bestpart.png"
  }
  ps[[1]] <- plot_borders_ggplot(clusters_id = z1s[,i], phillypoly = tracts, 
                                 w.sym = w.sym, var = alphas[,i], limits = limitsA, polyfortified = tracts_fort, 
                                 title = paste("Mean level"), palette_num = 3,
                                 legend = TRUE)
  ps[[2]] <- plot_borders_ggplot(clusters_id = z2s[,i], phillypoly = tracts, 
                                 w.sym = w.sym, var = betas[,i], limits = limitsB, polyfortified = tracts_fort,
                                 title = paste("Time trend"), palette_num = 5,
                                 legend = TRUE)
  
  ggsave(filename = output_filename, plot = arrangeGrob(grobs = ps, nrow = 1), device = "png", path = "figures/",
         width = 5, height = 2.5, units = "cm", scale = 6)
}
