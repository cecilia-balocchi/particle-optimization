source("scripts/plot_functions.R")
load("data/shapefile/phillytracts")
load("data/phillygooglemap.rdata")

# get the adjacency matrix from the shapefile
n_tr <- dim(tracts)[1]
list.poly <- poly2nb(tracts)
w <- matrix(data = 0, nrow = n_tr, ncol = n_tr)
for(i in 1:n_tr){
  w[i, list.poly[[i]]] <- 1
}
w.sym <- w + t(w) - w*t(w)
tracts_fort <- fortify(tracts)



input_string <- "results/density_L20priorsep_ep_eta1_lambda100task1.rdata"
load(input_string)
tmp <- lam_100$output  
L <- 20
w_sorted <- sort(tmp$w, decreasing = TRUE, index.return = TRUE)

alphas <- array(NA, dim = c(n_tr, L))
betas <- array(NA, dim = c(n_tr, L))
z1s <- array(NA, dim = c(n_tr, L))
z2s <- array(NA, dim = c(n_tr, L))
ps <- list()

for(jj in 1:L){
  j <- w_sorted$ix[jj]

  partitionA <- tmp$particle_set_A[[j]]
  partitionB <- tmp$particle_set_B[[j]]
  alphas[,jj] <- tmp$alpha_particle[,j] ## - log(2) ## !!!
  betas[,jj] <- tmp$beta_particle[,j]
  z1 <- numeric(n_tr)
  for(k in 1:length(partitionA)){
    z1[partitionA[[k]]] <- k
  }
  z2 <- numeric(n_tr)
  for(k in 1:length(partitionB)){
    z2[partitionB[[k]]] <- k
  }
  z1s[,jj] <- z1
  z2s[,jj] <- z2
}

limitsA <- c(min(alphas), max(alphas))
limitsB <- c(min(betas), max(betas))



# create single plot for alpha and beta (Figure 5)
jj = 1
legend_str <- ""
p1 <- plot_borders_ggplot(clusters_id = z1s[,jj], phillypoly = tracts, 
                               w.sym = w.sym, var = alphas[,jj], limits = limitsA, polyfortified = tracts_fort, 
                               title = paste("Mean level",legend_str), palette_num = 3,
                               legend = TRUE)
p2 <- plot_borders_ggplot(clusters_id = z2s[,jj], phillypoly = tracts, 
                               w.sym = w.sym, var = betas[,jj], limits = limitsB, polyfortified = tracts_fort,
                               title = paste("Time trend",legend_str), palette_num = 5,
                               legend = TRUE)
output_filename <-  paste0("density_tract_ep_ep")
ggsave(filename = paste0(output_filename, "_alpha.png"), plot = p1, device = "png", path = "figures/",
       width = 2.5, height = 2.5, units = "cm", scale = 6)

ggsave(filename = paste0(output_filename, "_beta.png"), plot = p2, device = "png", path = "figures/",
       width = 2.5, height = 2.5, units = "cm", scale = 6)



## this could be used to plot all the particles
for(jj in 1:L){ 
  output_filename <-  paste0("tract_density_ep_ep_L20_part",jj)
  legend_str <- paste("Particle",jj,"with w =",round(w_sorted$x[jj], 3))
  ps[[1]] <- plot_borders_ggplot(clusters_id = z1s[,jj], phillypoly = tracts, 
                                 w.sym = w.sym, var = alphas[,jj], limits = limitsA, polyfortified = tracts_fort, 
                                 title = paste("Mean level",legend_str), palette_num = 3,
                                 legend = TRUE)
  ps[[2]] <- plot_borders_ggplot(clusters_id = z2s[,jj], phillypoly = tracts, 
                                 w.sym = w.sym, var = betas[,jj], limits = limitsB, polyfortified = tracts_fort,
                                 title = paste("Time trend",legend_str), palette_num = 5,
                                 legend = TRUE)
  
  ggsave(filename = paste0(output_filename, ".png"), plot = arrangeGrob(grobs = ps, nrow = 1), device = "png", path = "figures/",
         width = 5, height = 2.5, units = "cm", scale = 6)
}