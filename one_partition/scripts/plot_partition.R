# this will only work for grids 
plot_partition_grid <- function(gamma, A_block, values = NULL, max_value = NULL, title = NA, title_cex = 0.8){
  col_list <- rev(brewer.pal(n = 5, name = "RdBu"))
  
  tmp <- partition_config(gamma, A_block)
  K <- tmp$K
  config <- tmp$config
  n_side <- sqrt(nrow(A_block))
  A_cluster <- tmp$A_cluster
  if(is.null(values)){
    values <- rep(0, times = length(A_block))
  }
  if(is.null(max_value)){
    max_value <- ceiling(max(abs(values)))
  }
  scaled_values <- rescale(values, to = c(0,1), from = max_value*c(-1,1))
  
  # Figure out where the borders are prior to plotting
  borders <- data.frame("unit1" = NA, "unit2" = NA)
  for(k in 1:K){
    cl_neighbors <- which(A_cluster[,k] == 1) # labels of clusters that are adjacenct to cluster k
    if(length(cl_neighbors) > 0){
      for(kk in cl_neighbors){
        
        if(config[k] > 1 & config[kk] > 1){
          tmp_border_index <- which(A_block[gamma[[kk]], gamma[[k]]] == 1, arr.ind = TRUE)
          border_index <- cbind(gamma[[kk]][tmp_border_index[,"row"]], gamma[[k]][tmp_border_index[,"col"]])
          
          # why sort? it provides information about whether it is a horizontal or vertical border
        } else if(config[k] > 1 & config[kk] == 1){ # cluster kk is a singleton
          tmp_border_index <- which(A_block[gamma[[kk]], gamma[[k]]] == 1) # these are column indices
          border_index <- cbind("unit1" = rep(gamma[[kk]], times = length(tmp_border_index)),
                                "unit2" = gamma[[k]][tmp_border_index])
        } else if(config[k] == 1 & config[kk] > 1){ # cluster k is a singleton
          tmp_border_index <- which(A_block[gamma[[kk]], gamma[[k]]] == 1)
          border_index <- cbind("unit1" = gamma[[kk]][tmp_border_index],
                                "unit2" = rep(gamma[[k]], times = length(tmp_border_index)))
        } else if(config[k] == 1 & config[kk] == 1){
          border_index <- cbind("unit1" = gamma[[kk]], "unit2" = gamma[[k]])
        }
        border_index <- t(apply(border_index, FUN = sort, MARGIN = 1)) # this just sorts things
        colnames(border_index) <- c("unit1", "unit2")
        borders <- rbind(borders, border_index)
      }
    }
  }
  draw_borders <- FALSE
  if(nrow(borders) > 1){
    borders <- borders[-1,] # drop the place-holder row
    draw_borders <- TRUE
  }
  plot(1, type = "n", xlim = c(0, n_side), ylim = c(0, n_side), xaxt = "n", yaxt = "n", 
       main = title, cex.main = title_cex, xlab = "", ylab = "", bty = "n")
  for(k in 1:K){
    i_vec <- ceiling(gamma[[k]]/n_side)
    j_vec <- gamma[[k]] - n_side*(i_vec- 1)
    for(ix in 1:config[k]){
      rect(j_vec[ix] - 1, i_vec[ix] - 1, j_vec[ix], i_vec[ix], border = "lightgray", lty = 1,lwd = 0.5,
           col = rgb(colorRamp(col_list, bias = 1)(scaled_values[j_vec[ix] + n_side*(i_vec[ix]-1)])/255))
    }
  }
  rect(0, 0, n_side, n_side, lwd = 2)
  if(draw_borders){
    
    # 1st observation: vertical borders will be between units whose indices differ by 1
    # 2nd observation: horizontal borders will be between units whose indices differ by n_side
    # 3rd observation: thanks to the sorting, borders[,2] - borders[,1] > 0
    n_borders <- nrow(borders)
    for(border_ix in 1:n_borders){
      i_1 <- ceiling(borders[border_ix,1]/n_side)
      j_1 <- borders[border_ix,1] - n_side * (i_1 - 1)
      # remember: (j_1, i_1) is the top-right corner of unit 1
      
      i_2 <- ceiling(borders[border_ix,2]/n_side)
      j_2 <- borders[border_ix,2] - n_side * (i_2 - 1) 
      
      diff <- borders[border_ix, 2] - borders[border_ix,1]
      # (j_2, i_2) is the top-right corder of unit 2
      if(diff == n_side){ # horizontal border
        # check to see that j_1 == j_2
        if(j_1 != j_2){
          print("something funny has happened!")
        }
        lines(x = c(j_1-1, j_1), y = c(i_1, i_1), lwd = 2)
      } else if(diff == 1){ # vertical border
        # check to see that i_1 == i_2
        if(i_1 != i_2){
          print("something funny has happened")
        }
        lines(x = c(j_1, j_1), y = c(i_1, i_1 - 1), lwd = 2)
      } else{ # something has gone horribly wrong
        print("something funny happened!")
      }
      
    }
  }
}


plot_particles_grid <- function(results, A_block, max_value, L_star = NULL, title = FALSE){
  if(is.null(L_star)) L_star <- length(results$particles)
  n_side <- ceiling(sqrt(L_star))
  col_list <- rev(brewer.pal(n = 5, name = "RdBu"))
  
  #png("~/Dropbox/Particle EM Spatial Clustering/PaperCombined/figures/test_layout.png", width = 6, height = 6, units = "in", res = 300)
  
  #png("~/Documents/ParticleOptPaper/figures/test_layout.png", width = 6, height = 6, units = "in", res = 300)
  
  layout(matrix(c(rep(1, times = n_side), 2:(n_side*n_side+1)), n_side+1, n_side, byrow = TRUE),
         heights = c(1/2.5, rep(1, times = n_side)))
  par(mar = c(1,1,1,1))
  plot(1, type = "n", xlim = c(0,1), ylim = c(0, 3), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = rgb(1,0,0,1/3))
  legend_seq <- seq(par("usr")[1], par("usr")[2], length = 500)
  for(leg_ix in 1:499){
    rect(legend_seq[leg_ix], par("usr")[3], legend_seq[leg_ix+1], par("usr")[4],
         border = NA, col = rgb(colorRamp(col_list, bias = 1)((leg_ix-1)/500)/255))
  }
  text(x = seq(par("usr")[1], par("usr")[2], length = 5), y = par("usr")[4]+1.0, cex = 0.9, 
       labels = round(seq(-max_value, max_value, length = 5), digits = 2), xpd = TRUE)
  for(l in 1:L_star){
    
    par(mar = c(1, 0.5, 1, 0.5), mgp = c(1.8, 0.5, 0), cex.main = 0.8)
    #if(results$counts[l] != 1) plot_title <- paste0("count = ", results$counts[l], ",  w = ", round(results$pstar[l], digits = 2), ",  log-post = ", round(results$log_post[l], digits = 2))
    #else plot_title <- paste0("w = ", round(results$pstar[l], digits = 3), ",  log-post = ", round(results$log_post[l], digits = 3))
    
    if(title == TRUE) plot_title <- paste("log-post =", round(results$log_post[l], digits = 2))
    else plot_title <- ""
    par(mar = c(1,0.5,1,0.5), mgp = c(1.8, 0.5, 0))
    plot_partition_grid(results$particles[[l]], A_block = A_block, values = results$alpha[,l], 
                        max_value = max_value, title = plot_title, title_cex = 0.75)
    
  }
  
  #dev.off()
  
  
  
}


plot_particle_trajectory <- function(iter, particle_set_traj, A_block, log_post){
  L_star <- length(particle_set[[iter]])
  n_side <- ceiling(sqrt(L_star))
  par(mar = c(1,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(n_side, n_side), cex.main = 0.7)
  for(l in 1:L_star) plot_partition_grid(particle_set_traj[[iter]][[l]], A_block, title = paste0("log-post = ", round(log_post[l,iter], digits = 3)))
  
}


