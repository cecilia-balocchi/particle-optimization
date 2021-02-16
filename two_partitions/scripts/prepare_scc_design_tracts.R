# Get the MST for the philly adjacency matrix

library(igraph)

W <- as.matrix(read.csv("data/A_block_density.csv", header = FALSE))
N <- nrow(W)

T <- 12
tmp_time <- scale(1:T, center = TRUE, scale = TRUE)

############
# Get MST & transformed design matrix
#   for full simulations, we need only compute this once and read it in
############

G <- graph_from_adjacency_matrix(W, mode = "undirected")

M <- mst(G)
mst_edges <- get.edges(M, es = E(M)) 

H <- matrix(0, N-1, N)
for(e in 1:nrow(mst_edges)){
  H[e, mst_edges[e,1]] <- 1
  H[e, mst_edges[e,2]] <- -1
}

H_tilde <- rbind(H, rep(1/N, times = N))
H_tilde_inv <- solve(H_tilde) 

X_tilde <- matrix(0, nrow = N*T, ncol = 2*N)
for(i in 1:N){
  tmp_X0 <- matrix(0, nrow = T, ncol = N)
  tmp_X1 <- matrix(0, nrow = T, ncol = N)
  
  tmp_X0[,i] <- 1
  tmp_X1[,i] <- tmp_time
  
  X_tilde[( (i-1)*T + 1):(i*T),1:N] <- tmp_X0 %*% H_tilde_inv
  X_tilde[( (i-1)*T + 1):(i*T), (N+1):(2*N)] <- tmp_X1 %*% H_tilde_inv
  
}
penalty_factor <- rep(1, times = 2*N)
penalty_factor[c(N, 2*N)] <- 0

#########################################
# Prepare the training and testing data
#########################################
Y_train_orig <- as.matrix(read.csv(file = "data/Y_density.csv", header = FALSE))
Y_train <- rep(NA, times = N*T)
for(i in 1:N){
  Y_train[(1 + (i-1)*T):(i*T)] <- as.vector(Y_train_orig[i,1:T])
}
  

rm(Y_train_orig)

# Read in data for 2018
load("data/crime_tract2019.rdata")
Y_test <- crime_tract2019$tr.density[crime_tract2019$year == 2018]

save(Y_train, Y_test, N, T, X_tilde, H, H_tilde, H_tilde_inv, mst_edges, penalty_factor, file = "data/ssc_philly.RData")


