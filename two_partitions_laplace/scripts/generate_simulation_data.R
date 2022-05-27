###############################
# Create two partitions
#  gamma_alpha will have 10 clusters, many of which are singletons situated in larger clusters
#  gamma_beta will be the 5-cluster flag partition

# Idea:when generating data from a CAR--within--cluster model, we will pick cluster centers
# in such a way that singletons in gamma_alpha will always have the same mean
# also clusters 1 and 6 will have the same mean
# cluster 4 and 8 will have same mean
# cluster 10 will have mean 0.
###############################
library(igraph)
library(MASS)

n_side <- 10
N <- n_side^2
g <- make_lattice(length = n_side, dim = 2)

W <- as_adj(g, type = "both")
W <- matrix(W, nrow = N)

# Define the four quadrants
quad_1 <- rep(1:5, times = 5) + 10*rep(0:4, each = 5)
quad_2 <- rep(6:10, times = 5) + 10*rep(0:4, each = 5)
quad_3 <- rep(1:5, times = 5) + 10*rep(5:9, each = 5)
quad_4 <- rep(6:10, times = 5) + 10*rep(5:9, each = 5)

alpha_cluster1 <- c(1,2,3,11,13,21,22,23)
alpha_cluster2 <- c(12)
alpha_cluster3 <- c(34,35,44,25,26)
alpha_cluster4 <- c(100,90,89,78,79,80)
alpha_cluster5 <- quad_4[!(quad_4 %in% c(alpha_cluster4))]
alpha_cluster6 <- c(73)
alpha_cluster7 <- (1:100)[!(1:100) %in% c(alpha_cluster1, alpha_cluster2, alpha_cluster3,
                                          alpha_cluster4, alpha_cluster5, alpha_cluster6)]


gamma_alpha <- list(alpha_cluster1, alpha_cluster2, alpha_cluster3, alpha_cluster4, alpha_cluster5,
                    alpha_cluster6, alpha_cluster7)

# plot_clusters_grid(gamma_alpha, W)

beta_cluster1 <- c(22,23,32,33,43,44,53,54,62,63,72,73)
beta_cluster2 <- quad_2
beta_cluster3 <- quad_4
beta_cluster4 <- (1:100)[!(1:100) %in% c(beta_cluster1,beta_cluster2,beta_cluster3)]

gamma_beta <- list(beta_cluster1, beta_cluster2,beta_cluster3,beta_cluster4)

# plot_clusters_grid(gamma_beta, W)

gamma_1 <- list(1:N)
gamma_n <- as.list(1:N)

save(W, gamma_alpha, gamma_beta, 
     gamma_1, gamma_n, N, 
     file = "data/partitions.RData")

#### Generate data ####
tot_sim <- 3            ## how many configurations ( = length(Delta_alpha))

seed_seq <- c(129, 724, 603, 212, 1123, 2391) # c(129, 724, 603, 212, 1123, 2391, 815, 1947)
n_sim <- 100
Delta_alpha <- c(2, 1, 0)
Delta_beta <- c(1, 0.5, 0)
rho <- 0.95
a1 <- 1/8 #1/10
b1 <- 1/8 # 1/3
sigma <- 0.25
sigma2 <- sigma^2
## Note: in practice, within cluster sd will be 0.25/sqrt(8) =~ 0.09

for(sim_number in 1:tot_sim){

  ## these values were tested with simulation2
  alphabar <- 3.5 + c(-1, 1, 0.5, 0, -0.5, 
                      1, 0) * Delta_alpha[sim_number]
  betabar <- c(-1, 0,-1, 1) * Delta_beta[sim_number]
  
  ## each column of alpha is for a different simulation (a total of n_sim) within a configuration (1:tot_sim)
  alpha <- matrix(NA, nrow = N, ncol = n_sim)
  beta <- matrix(NA, nrow = N, ncol = n_sim)
  
  for(k in 1:length(gamma_alpha)){
    cl_index <- gamma_alpha[[k]]
    n_k <- length(cl_index)
    if(n_k==1){
      Sigma_alpha <- 1/(1-rho)
    } else {
      W_k <- W[cl_index, cl_index]
      D <- diag(rowSums(W_k))
      W_star_k <- D - W_k
      Omega_alpha <- rho * W_star_k + (1 - rho) * diag(n_k)
      Sigma_alpha <- solve(Omega_alpha)
    }
    set.seed(seed_seq[sim_number] + k)
    
    alpha[cl_index,] <- t(matrix(mvrnorm(n = n_sim, 
                                         mu = rep(alphabar[k], times = n_k), 
                                         Sigma = a1 * sigma2 * Sigma_alpha), 
                                 nrow = n_sim))
  }
  for(k in 1:length(gamma_beta)){
    cl_index <- gamma_beta[[k]]
    n_k <- length(cl_index)
    if(n_k==1){
      Sigma_alpha <- 1/(1-rho)
    } else {
      W_k <- W[cl_index, cl_index]
      D <- diag(rowSums(W_k))
      W_star_k <- D - W_k
      Omega_alpha <- rho * W_star_k + (1 - rho) * diag(n_k)
      Sigma_alpha <- solve(Omega_alpha)
    }
    set.seed(seed_seq[sim_number] + k)
    
    beta[cl_index,] <- t(matrix(mvrnorm(n = n_sim, 
                                        mu = rep(betabar[k], times = n_k), 
                                        Sigma = b1 * sigma2 * Sigma_alpha), 
                                nrow = n_sim))
  }
  assign(paste0("alpha_", sim_number), alpha)
  assign(paste0("beta_", sim_number), beta)
}

t <- 12
tmp_time <- scale(1:t, center = TRUE, scale = TRUE)


X_i <- cbind(rep(1, times = N*t),rep(tmp_time, times = N))
X <- matrix(rep(tmp_time, each = N), ncol = length(tmp_time), nrow = N)


varlist <- c("sigma2", "rho", "a1", "b1",
             "X_i", "X", 
             paste0("alpha_",1:tot_sim),paste0("beta_",1:tot_sim))
save(list = varlist,
     file ="data/alphabetas.RData")
