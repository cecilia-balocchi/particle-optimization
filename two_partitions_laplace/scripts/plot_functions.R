require(sp)
require(spdep)
require(rgeos)
require(RColorBrewer)
require(scales)
require(gridExtra)
require(ggplot2)
require(dplyr)
require(grid)

plot_partition_grid <- function(gamma, A_block, values = NULL, max_value = NULL, min_value = NULL, 
                                title = NA, title_cex = 0.8, colorscheme = "RdBu")
{
  # other options for colorscheme is PRGn
  col_list <- rev(brewer.pal(n = 5, name = colorscheme))
  
  
  tmp <- get_config(gamma, A_block)
  K <- tmp$num_clusters
  config <- tmp$cluster_size
  n_side <- sqrt(nrow(A_block))
  A_cluster <- tmp$A_cluster
  if(is.null(values)){
    values <- rep(0, times = length(A_block))
  }
  if(is.null(max_value)){
    max_value <- ceiling(max(abs(values)))
  }
  if(is.null(min_value)){
    min_value <- -max_value
  }
  scaled_values <- rescale(values, to = c(0,1), from = c(min_value,max_value))
  
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


plot_clusters_grid <- function(gamma, A_block, title = NA, title_cex = 0.8)
{
  
  tmp <- get_config(gamma, A_block)
  K <- tmp$num_clusters
  config <- tmp$cluster_size
  n_side <- sqrt(nrow(A_block))
  A_cluster <- tmp$A_cluster
  
  if(K < 11){
    col_list <- rev(brewer.pal(n = K, name = "Spectral"))
  } else stop("Too many clusters")
  
  
  values <- rep(NA, times = nrow(A_block))
  for(k in 1:K) values[gamma[[k]]] <- k
  
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
           col = col_list[k])
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

# Configuration of the parttion
get_config <- function(gamma, A_block){
  tmp <- sapply(gamma, FUN = length)
  singletons <- which(tmp == 1)
  non_singletons <- which(tmp != 1)
  K <- length(gamma) # number of clusters
  n <- nrow(A_block) # number of block groups
  n_singletons <- length(singletons)
  n_non_singletons <- length(non_singletons)
  Z <- matrix(0, nrow = n, ncol = n) # pairwise co-allocation matrix
  
  for(k in 1:K){
    Z[gamma[[k]], gamma[[k]]] <- 1
  }
  
  if(K > 1){ # identify which clusters are spatially adjacent
    A_cluster <- matrix(0, nrow = K, ncol = K)
    for(k in 1:(K-1)){
      for(kk in (k+1):K){
        A_cluster[k,kk] <- 1*any(A_block[gamma[[k]], gamma[[kk]]] == 1)
        A_cluster[kk,k] <- A_cluster[k,kk]
      }
    }
  } else{
    A_cluster <- NULL
  }
  return(list(num_clusters = K, cluster_size = tmp, 
              singletons = singletons, non_singletons = non_singletons,
              A_cluster = A_cluster, Z = Z))
}