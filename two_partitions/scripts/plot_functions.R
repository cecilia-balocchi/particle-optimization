require(sp)
require(spdep)
require(rgeos)
require(RColorBrewer)
require(scales)
require(gridExtra)
require(ggplot2)
require(ggmap)
require(dplyr)
require(grid)
# import the shapefile
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


## this function is used to make a plot of the map with the borders of the clusters, and the colors for the parameters of each neighborhood
plot_borders_ggplot <- function(clusters_id, phillypoly, w.sym, var = data$bgs_betas, limits = NULL, polyfortified = NULL, title = NULL, legend = FALSE, legend_name = NULL, palette_num = 5, map = NULL){
  # clusters_id is a vector, each element is the cluster assignment of each region
  # phillypoly is the shapefile
  # w.sym is the adjacency matrix
  # var is the variable we want to plot
  # limits is a vector of length 2 containing the limits of the color scale 
  # polyfortified is the object returned by fortify, so that if the function is repeated we don't need to call fortify all the times
  # title is a string for the plot title
  # legend is a boolean to plot (or not) the legend
  # legend_name is a string for the legend
  # palette_num is color palette number for fill
  # map can be the googlemap (I saved it once because now the google map API is not free anymore, but can be downloaded with ggmap::get_googlemap)
  n_tr <- dim(phillypoly)[1]
  group_ind <- (1:n_tr)
  
  if(is.null(polyfortified)){
    polyfortified <- fortify(phillypoly)
  }
  if(is.null(limits)){
    limits <- c(min(var, na.rm=T), max(var, na.rm=T))
  }
  MI <- limits[1]
  MA <- limits[2]
  mi <- min(var, na.rm=T)
  ma <- max(var, na.rm=T)
  z <- (0-mi)/(ma-mi)
  z_adj<-(ma-mi)/(MA-MI) * z + (mi-MI)/(MA-MI)
  
  bl <- border_list(clusters_id = clusters_id, index = group_ind, w = w.sym)
  if(length(bl)>0){
    SL1 <- border_SL_list(bl, phillypoly)
    SL_complete <- SpatialLines(SL1)
    mySLDF <- SpatialLinesDataFrame(SL_complete, data = data.frame(ID = rep(1,length(SL_complete)), row.names = 1:length(SL_complete)))
  } 
  
  my.data <- data.frame(id = as.numeric(rownames(phillypoly@data)), value =var)
  my.data$id <- as.character(my.data$id)
  plotData <- left_join(polyfortified,my.data, by = "id")
  if(!is.null(map)){
    p1 <- ggmap(map,extent = "device",ylab = "Latitude",
                xlab = "Longitude",color = "bw") + theme_bw() + 
      theme(panel.border = element_blank()) + 
      geom_polygon(data=plotData, aes(x=long,y=lat,group=group, fill = value),
                   color=alpha("black",0.15), alpha =1) 
  } else {
    p1 <- ggplot()+theme_void() + 
      geom_polygon(data=plotData, aes(x=long,y=lat,group=group, fill = value),
                   color=alpha("red",0.15), alpha =1) 
  }
  if(legend){
    p1 <- p1 +
      scale_fill_distiller(type = "div", palette = palette_num, limits = c(MI,MA), values = c(0,z_adj,1), name = legend_name) + 
      theme(plot.title = element_text(hjust = 0.5, vjust = -1, size=18),
                         axis.ticks = element_blank(),axis.text = element_blank(),
                         panel.grid = element_blank(), axis.title = element_blank(),
                         legend.text = element_text(size=14), legend.title = element_text(size=14))
  } else {
    p1 <- p1 +
      scale_fill_distiller(type = "div", palette = palette_num, limits = c(MI,MA), values = c(0,z_adj,1), guide= FALSE) + 
      theme(plot.title = element_text(hjust = 0.5, vjust = -1, size=18),
                         axis.ticks = element_blank(),axis.text = element_blank(),
                         panel.grid = element_blank(), axis.title = element_blank())
  }
  
  if(!is.null(title)){
    p1 <- p1 + labs(title = title)
  }
  if(length(bl)>0){
    p1 <- p1 + geom_line(data = mySLDF, aes(x=long,y=lat,group=group), color=alpha("black",1), size = 0.75)
  }
  return(p1)
}

## this function is used to make the difference plots
plot_borders_ggplot_diff <- function(clusters_id, phillypoly, w.sym, diff_vec = NULL, polyfortified = NULL, weight = NULL, title = NULL, map = NULL){
  # clusters_id is a vector, each element is the cluster assignment of each region
  # phillypoly is the shapefile
  # w.sym is the adjacency matrix
  # diff_vec is the output from find_diff_vec
  # polyfortified is the object returned by fortify, so that if the function is repeated we don't need to call fortify all the times
  # weight can be the particle weight that can be printed in the title
  # title is a string for the plot title
  # map can be the googlemap (I saved it once because now the google map API is not free anymore, but can be downloaded with ggmap::get_googlemap)
  n_tr <- dim(phillypoly)[1]
  group_ind <- (1:n_tr) 
  
  if(!is.null(weight)){
    string <- paste0("w = ",format(weight, digits = 3,scientific = T))
  } else {
    string <- NULL
  }
  if(is.null(polyfortified)){
    polyfortified <- fortify(phillypoly)
  }
  
  bl <- border_list(clusters_id = clusters_id, index = group_ind, w = w.sym)
  if(length(bl)>0){
    SL1 <- border_SL_list(bl, phillypoly)
    SL_complete <- SpatialLines(SL1)
    mySLDF <- SpatialLinesDataFrame(SL_complete, data = data.frame(ID = rep(1,length(SL_complete)), row.names = 1:length(SL_complete)))
  } 
  
  plotData <- polyfortified
  if(!is.null(map)){
    p1 <- ggmap(map,extent = "device",ylab = "Latitude",
                xlab = "Longitude",color = "bw") + theme_bw() +
      geom_polygon(data=plotData, aes(x=long,y=lat,group=group),
                   color=alpha("black",0.15), alpha =0.5, fill = "white")
  } else {
    p1 <- ggplot()+ theme_void() +
      geom_polygon(data=plotData, aes(x=long,y=lat,group=group),
                   color=alpha("black",0.15), alpha =0)
  }
  p1 <- p1 + 
    theme(plot.title = element_text(hjust = 0.5, vjust = -1, size=18),
          axis.ticks = element_blank(),axis.text = element_blank(),
          panel.grid = element_blank(), axis.title = element_blank())
  if(!is.null(title)){
    p1 <- p1 + labs(title = title)
  }
  if(length(bl)>0){
    p1 <- p1 + geom_line(data = mySLDF, aes(x=long,y=lat,group=group), color="black", size = 0.75, alpha = 1)
  }
  if(!is.null(diff_vec)){
    my.data2 <- data.frame(id = as.numeric(rownames(phillypoly@data)), value = diff_vec)
    my.data2$id <- as.character(my.data2$id)
    plotData2 <- left_join(polyfortified,my.data2, by = "id")
    p1 <- p1 + geom_polygon(data=plotData2, aes(x=long,y=lat,group=group, fill = as.factor(value)),
                            color = alpha("black",0), alpha =0.5) + scale_fill_manual(values = c("black"), guide=FALSE)
  }
  return(p1)
}

## Function that from the path of the output plots the three best particles and the difference plots
plot_particles_diff <- function(strings, alpha_shift = 0){
  # strings can be a vector of path strings if you want to combine more than one set of particles and
  # combine the best ones. We will use it with only one string
  # this code plots the best three (3) particles, with their differences in between. 
  # a larger number of particles would be difficult to visualize in a paper.
  sorted_x <- c()
  hhs <- c()
  for(hh in 1:length(strings)){
    load(strings[hh])
    tmp_sorted <- sort(unique(tmp$log_posterior), index.return = T, decreasing = T)
    sorted_x <- c(sorted_x, tmp_sorted$x[1:3])
    hhs <- c(hhs, rep(hh, 3))
  }
  hhs <- hhs[match(unique(sorted_x), sorted_x)]
  sorted_x <- unique(sorted_x)
  doublesorted_ix <- sort(sorted_x, index.return = T, decreasing = T)$ix[1:3]
  n <- dim(tmp$alpha_particle)[1]
  z1s <- array(NA, dim = c(n, 3))
  z2s <- array(NA, dim = c(n, 3))
  alphas <- array(NA, dim = c(n, 3))
  betas <- array(NA, dim = c(n, 3))
  i <- 1
  ps <- list()
  newrow <- 6
  tracts_fort <- fortify(tracts)
  for(ii in doublesorted_ix){
    hh <- hhs[ii]
    load(strings[hh])
    j <- which(tmp$log_posterior == sorted_x[ii])[1]
    cat(ii, " ", hh, " ", j, "\n")
    partitionA <- tmp$particle_set_A[[j]]
    partitionB <- tmp$particle_set_B[[j]]
    alphas[,i] <- tmp$alpha_particle[,j] + alpha_shift
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
    i <- i + 1
  }
  limitsA <- c(min(alphas), max(alphas))
  limitsB <- c(min(betas), max(betas))
  for(i in 1:3){
    z1 <- z1s[,i]
    z2 <- z2s[,i]
    legend <- FALSE
    legend_name <- NULL
    if(i == 3){
      legend <- TRUE
      legend_name <- ""
    }
    title <- paste0("w = ", format(tmp$log_posterior[j], digits = 6))
    plot_tmp <- plot_borders_ggplot(clusters_id = z1, phillypoly = tracts, w.sym = w.sym, 
                                    var = alphas[,i], limits = limitsA, polyfortified = tracts_fort, 
                                    title = paste("Mean level\nParticle",i), palette_num = 3,
                                    legend = legend, legend_name = legend_name)
    if(i == 3){
      leg <- g_legend(plot_tmp)
      ps[[(2*(i-1))+1]] <- plot_tmp+theme(legend.position = 'none')
      ps[[6]] <- leg
      legend_name <- ""
    } else {
      ps[[(2*(i-1))+1]] <- plot_tmp
    }
    plot_tmp <- plot_borders_ggplot(clusters_id = z2, phillypoly = tracts, w.sym = w.sym, 
                                    var = betas[,i], limits = limitsB, polyfortified = tracts_fort,
                                    title = paste("Time trend\nParticle",i), palette_num = 5,
                                    legend = legend, legend_name = legend_name)
    if(i == 3){
      leg <- g_legend(plot_tmp)
      ps[[newrow + (2*(i-1))+1]] <- plot_tmp+theme(legend.position = 'none')
      ps[[newrow + 6]] <- leg
    } else {
      ps[[newrow + (2*(i-1))+1]] <- plot_tmp
    }
  }
  for(i in 1:2){
    z1 <- z1s[,1]
    titlestr <- paste("\nDifference from",1,"to",i+1)
    z2 <- z1s[,i+1]
    PW1 <- pairwise(z1)
    PW2 <- pairwise(z2)
    diff <- PW1 - PW2
    tmp2 <- find_diff_vec(diff)
    tmp2NA <- tmp2
    tmp2NA[tmp2NA == 0] <- NA
    z1 <- rep(1, n)
    ps[[2*i]] <- plot_borders_ggplot_diff(clusters_id = z1, phillypoly = tracts, w.sym = w.sym,
                                          diff_vec = tmp2NA, polyfortified = tracts_fort, 
                                          title = titlestr)
    z1 <- z2s[,1] 
    titlestr <- paste("\nDifference from",1,"to",i+1)
    z2 <- z2s[,i+1]
    PW1 <- pairwise(z1)
    PW2 <- pairwise(z2)
    diff <- PW1 - PW2
    tmp2 <- find_diff_vec(diff)
    tmp2NA <- tmp2
    tmp2NA[tmp2NA == 0] <- NA
    z1 <- rep(1, n)
    ps[[newrow + 2*i]] <- plot_borders_ggplot_diff(clusters_id = z1, phillypoly = tracts, w.sym = w.sym,
                                                   diff_vec = tmp2NA, polyfortified = tracts_fort,
                                                   title = titlestr)
  }
  return(ps)
}

## Combine all the different functions into one plot
combine_diff_one <- function(strings, k_combined = 10, AorB = TRUE, map = NULL){
  # k_combined: how many particles
  # AorB: are we using particleA or particleB (alpha or beta) 
  # map: googlemap from ggmap::get_googlemap
  load(strings)
  # this is written for particleA
  logpost <- numeric(k_combined)
  
  sorted_x <- c()
  tmp_sorted <- sort(unique(tmp$log_posterior), index.return = T, decreasing = T)
  sorted_x <- c(sorted_x, tmp_sorted$x)
  
  sorted_x <- unique(sorted_x)
  doublesorted_ix <- sort(sorted_x, index.return = T, decreasing = T)$ix[1:k_combined]
  n <- dim(tmp$alpha_particle)[1]
  
  diff_overall <- numeric(n)
  ii <- doublesorted_ix[1]
  j <- which(tmp$log_posterior == sorted_x[ii])[1]
  if(AorB){
    partitionA <- tmp$particle_set_A[[j]]
    titlestr <- "Mean level"
  } else {
    partitionA <- tmp$particle_set_B[[j]]
    titlestr <- "Time trend"
  }
  z1 <- numeric(n)
  for(k in 1:length(partitionA)){
    z1[partitionA[[k]]] <- k
  }
  PW1 <- pairwise(z1)
  logpost[1] <- tmp$log_posterior[j]
  for(jj in 2:k_combined){
    ii <- doublesorted_ix[jj]
    j <- which(tmp$log_posterior == sorted_x[ii])[1]
    if(AorB){
      partitionB <- tmp$particle_set_A[[j]]
    } else {
      partitionB <- tmp$particle_set_B[[j]]
    }
    z2 <- numeric(n)
    for(k in 1:length(partitionB)){
      z2[partitionB[[k]]] <- k
    }
    PW2 <- pairwise(z2)
    logpost[jj] <- tmp$log_posterior[j]
    diff <- PW1 - PW2
    tmp2 <- find_diff_vec(diff)
    diff_overall[which(tmp2 == 1)] <- 1
  }
  diff_overall[diff_overall == 0] <- NA
  p <- plot_borders_ggplot_diff(z1, tracts, w.sym, diff_vec = diff_overall, map = map, title = titlestr)
  return(list(plot = p, logposts = logpost))
}


## Given a cluster_id vector (each element of the vector contains a unique numeric identifier of the cluster)
## it returns the list of cluster borders as a SpatialLines object, to plot the partition
border_list <- function(clusters_id, index, w){
  a <- list()
  tmp <- 1
  for(i in 1:length(index)){
    js <- which(w[i,]==1)
    for(j in js){
      if(clusters_id[i] != clusters_id[j]){
        if(j > i){
          a[[tmp]] <- c(index[i],index[j])
          tmp <- tmp + 1
        }
      }
    }
  }
  a
}
border_SL_list <- function(border_list, phillypoly){
  l <- list()
  iter <- 1
  for(n in 1:length(border_list)){
    # i,j have 0-indexing
    i <- border_list[[n]][1]
    j <- border_list[[n]][2]
    # cat(n,":",i,j,"\n")
    A <- SpatialPolygons(list(Polygons(list(phillypoly@polygons[[i]]@Polygons[[1]]),ID=1)))
    B <- SpatialPolygons(list(Polygons(list(phillypoly@polygons[[j]]@Polygons[[1]]),ID=1)))
    tmp <- gIntersection(A,B)
    if(class(tmp)[1] == "SpatialLines"){
      tmp1 <- tmp@lines[[1]]
      tmp1@ID <- as.character(iter)
      l[[iter]] <- tmp1
      iter <- iter + 1
    }
  }
  l
}

### Functions to find the difference between two partitions
## From the cluster assignment vector, we find the pairwise matrix
pairwise <- function(z1){
  n <- length(z1)
  PW <- array(0, dim = c(n,n))
  for(i in unique(z1)){
    index <- which(z1 == i)
    PW[index, index] <- 1
  }
  return(PW)
}
## From the difference of the pairwise matrices of two partitions, find neighborhoods with different cluster assignment
find_diff_vec <- function(diff){
  # diff is a matrix equal to the difference of the pairwise matrices of two partitions
  
  ## Algorithm
  # let z1 and z2 be the cluster assignment vectors of partition1 and partition2
  # diff = pairwise(z1) - pairwise(z2) (it's the difference of pairwise matrices)
  # run through each row i of diff
  # count how many elements are different from 0 in row i; this correspond to how many elements j have z1(j) == z1(i) but z2(j) != z2(i) 
  # or z1(j) != z1(i) but z2(j) == z2(i). Such i certainly changed is cluster assignment from z1 to z2 (but also j)
  # identify the element i that has the largest of such number, because we want to represent the minimal difference between z1 and z2, call it x
  # let's find who else has the same diff row as x, and call it cluster: these are other elements that had the same behavior as x
  # and can be considered in the same cluster (it's not really its cluster, but if x was removed together with "cluster" from an originally larger cluster
  # then we can think of "cluster" as x's cluster)
  # "cluster" (which contains x) is highlighted as part of the difference between the two partitions
  # we remove the -1 or +1 from diff that were associated with cluster
  # and repeat
  n <- dim(diff)[1]
  vec_diff <- numeric(n)
  # flag tells us if we reached the end of the loop
  flag <- TRUE
  while(flag){
    x <- which.max(apply(X = diff, MARGIN = 1, FUN = function(x){ sum(x < 0 | x > 0) }))
    ref <- diff[x,]
    cluster <- which(apply(X = diff, MARGIN = 1, FUN = function(x){ all.equal(x,ref) == TRUE }))
    if(length(cluster) == n){
      flag <- FALSE
    } else {
      vec_diff[cluster] <- 1
      diff[which(ref != 0), cluster] <- 0
      diff[cluster, which(ref != 0)] <- 0
    }
  }
  return(vec_diff)
}


## function to return only the legend (to separate the legend from the plot)
## found here: https://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plot_partition_grid <- function(gamma, A_block, values = NULL, max_value = NULL, min_value = NULL, 
                                title = NA, title_cex = 0.8, colorscheme = "RdBu"){
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


plot_clusters_grid <- function(gamma, A_block, title = NA, title_cex = 0.8){
  
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