rand_index <- function(N, gamma1, gamma2){
  tmp_binder <- binder_loss(N, gamma1, gamma2)
  tmp_rand <- 1 - binder_loss(N, gamma1, gamma2) * 2/(N * (N-1))
  return(tmp_rand)
}

binder_loss <- function(N, gamma1, gamma2){
  
  k1 <- length(gamma1)
  k2 <- length(gamma2)
  counts <- matrix(nrow = k1, ncol = k2)
  # need matrix of counts n_ij counting the number of indices that belong to cluster k in gamma1 and cluster kk in gamma2
  for(k in 1:k1){
    for(kk in 1:k2){
      counts[k,kk] <- length(intersect(gamma1[[k]], gamma2[[kk]]))
    }
  }
  row_sums <- rowSums(counts)
  col_sums <- colSums(counts)
  dist <- 0.5 * (sum( (row_sums)^2) + sum( (col_sums)^2) - 2 * sum(counts^2))
  
  return(dist)
  
}


get_tot_ss <- function(alpha, beta, gamma_A, gamma_B){
  tot_ss <- 0
  K <- length(gamma_A)
  for(k in 1:K) tot_ss <- tot_ss + sum( (alpha[gamma_A[[k]]] - mean(alpha[gamma_A[[k]]]))^2)
  K <- length(gamma_B)
  for(k in 1:K) tot_ss <- tot_ss + sum( (beta[gamma_B[[k]]] - mean(beta[gamma_B[[k]]]))^2 )
  return(tot_ss)
}

rand_adj_vi <- function(N, gamma1, gamma2){
  
  k1 <- length(gamma1)
  k2 <- length(gamma2)
  counts <- matrix(nrow = k1, ncol = k2)
  # need matrix of counts n_ij counting the number of indices that belong to cluster k in gamma1 and cluster kk in gamma2
  for(k in 1:k1){
    for(kk in 1:k2){
      counts[k,kk] <- length(intersect(gamma1[[k]], gamma2[[kk]]))
    }
  }
  row_sums <- rowSums(counts)
  col_sums <- colSums(counts)
  
  dist <- 0.5 * (sum( (row_sums)^2) + sum( (col_sums)^2) - 2 * sum(counts^2))
  tmp_rand <- 1 - dist * 2/(N * (N-1))

  choose_row <- choose(row_sums, 2)
  choose_col <- choose(col_sums, 2)

  correct <- sum(choose_row) * sum(choose_col)/choose(N, 2)
  rand_adj <- (sum(choose(counts, 2)) - correct)/(0.5 * sum(choose_row) + 0.5 * sum(choose_col) - correct)

  base_log <- 2
  vi <- sum( (row_sums/N) * log(row_sums/N, base = base_log) ) + 
        sum( (col_sums/N) * log(col_sums/N, base = base_log) ) - 2* sum( (counts/N) * log(counts/N, base = base_log), na.rm = T )

  return(list(rand = tmp_rand, rand_adj = rand_adj, vi = vi))
  
}

partition_list2z <- function(gamma){
  k_tot <- length(gamma)
  z <- numeric(length(unlist(gamma)))
  count <- 1
  for(k in 1:k_tot){
    z[gamma[[k]]] <- count
    count <- count + 1
  }
  return(z)
}

partition_z2list <- function(z){
  unik_z <- unique(z)
  gamma <- list()
  for(k in 1:length(unik_z)){
    gamma[[k]] <- which(z == unik_z[k])
  }
  return(gamma)
}

Modify_partition <- function(z, W){
  n <- length(z)
  z_unik <- unique(z)
  k_tot <- length(z_unik)
  new_k <- k_tot + 1
  new_z <- z
  for(k in 1:k_tot){
    cl_ind <- which(z == z_unik[k])
    if(length(cl_ind)>1){
      W_cl <- W[cl_ind,cl_ind]
      tmp = Connected_components(W_cl)
      if(tmp$count > 1){
        # modify new_z for the sub clusters identified by the components
        for(comp in 1:tmp$count){
          new_z[cl_ind[which(tmp$components == comp)]] <- new_k
          new_k <- new_k + 1
        }
      }
    }
  }
  new_z_sorted <- new_z
  new_z_unik <- unique(new_z)
  count <- 1
  for(k in new_z_unik){
    new_z_sorted[new_z == k] <- count
    count <- count + 1
  }

  return(new_z_sorted)
}

Connected_components <- function(W_cluster){
  n <- nrow(W_cluster)
  visited <- rep(0, n)
  components <- rep(NA, n)
  count <- 0
  for(v in 1:n){
    if(visited[v] == 0){
      count <- count + 1
      tmp <- DFSUtil(W_cluster, n, v, visited, components, count)
      visited <- tmp$visited
      components <- tmp$components
    }
  }
  return(list(components = components, count = count))
}

DFSUtil <- function(W_cluster, n, v, visited, components, count){
  visited[v] <- 1
  components[v] <- count
  for(i in 1:n){
    if(W_cluster[v,i] == 1){
      if(visited[i] == 0){
        tmp <- DFSUtil(W_cluster, n, i, visited, components, count)
        visited <- tmp$visited
        components <- tmp$components
      }
    }
  } 
  return(list(visited = visited, components = components))
}

spectral_clustering <- function(x, W, k){
  n <- length(x)
  if(n == 1){
    warning("Spectral clustering: n is 1")
  }
  dist_m <- as.matrix(dist(x))
  similarity_m <- exp(-dist_m/(2 * sd(x)))
  W_beta <- diag(n) + similarity_m * W # element-wise multiplication
  Dinv_sqrt <- diag(1/sqrt(rowSums(W_beta)))
  L <- diag(n) - Dinv_sqrt %*% W_beta %*% Dinv_sqrt
  tmp <- eigen(L)
  U <- tmp$vectors[,(n-k+1):n, drop = FALSE]
  cat(paste0("Spectral clustering: k ",k," dim(U)",nrow(U)," ",ncol(U),"\n"))
  U <- diag(1/sqrt(rowSums(U^2))) %*% U
  km <- kmeans(x = U, centers = k, nstart = 25)
  return(km)
}

elbow_method <- function(x, Kmax){
  x <- x/x[1]
  possible_k <- 2:(Kmax-1)
  diff1 <- rep(NA, times = length(possible_k))
  diff2 <- rep(NA, times = length(possible_k))
  
  for(k_ix in 1:length(possible_k)){
    diff1[k_ix] <- x[possible_k[k_ix]] - x[possible_k[k_ix]-1]
    diff2[k_ix] <- x[possible_k[k_ix]+1] - 2*x[possible_k[k_ix]] + x[possible_k[k_ix]-1]
  }
  kappa <- abs(diff2)/(1 + diff1^2)^(3/2)
  return(possible_k[which.max(kappa)])
}

get_ssc_partition <- function(mst_edges, H, params_hat, tol = 1e-9){
  
  # zero entries in H * params_hat mean that beta_i - beta_j 
  # are equal (and so belong in the same cluster)
  # we will create an adjacency matrix
  
  N <- length(params_hat)
  tmp_adj <- matrix(0, nrow = N, ncol = N)
  diffs <- H %*% params_hat
  for(k in 1:length(diffs)){
    if(abs(diffs[k]) < tol){
      tmp_adj[mst_edges[k,1], mst_edges[k,2]] <- 1
      tmp_adj[mst_edges[k,2], mst_edges[k,1]] <- 1
    }
  }
  
  tmp_G <- graph_from_adjacency_matrix(tmp_adj, mode = "undirected")
  tmp_connected <- components(tmp_G)
  
  gamma_hat <- list()
  for(k in 1:tmp_connected$no){
    gamma_hat[[k]] <- sort(which(tmp_connected$membership == k))
  }
  return(gamma_hat)
}