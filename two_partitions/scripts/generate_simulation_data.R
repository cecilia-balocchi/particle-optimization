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

n_side <- 20
N <- 400
g <- make_lattice(length = n_side, dim = 2)

W <- as_adj(g, type = "both")
W <- matrix(W, nrow = N)


# Define the four quadrants
quad_1 <- rep(11:20, times = 10) + 20*rep(10:19, each = 10)
quad_2 <- rep(1:10, times = 10) + 20*rep(10:19, each = 10)
quad_3 <- rep(1:10, times = 10) + 20*rep(0:9, each = 10)
quad_4 <- rep(11:20, times = 10) + 20*rep(0:9, each = 10)

# Define the 10 clusters for alpha 
alpha_cluster1 <- rep(1:5, times = 5) + 20 * rep(0:4, each = 5) # lower left 5 x 5 grid
alpha_cluster1 <- alpha_cluster1[alpha_cluster1 != 43]
alpha_cluster2 <- c(43)
alpha_cluster3 <- 361:400
alpha_cluster4 <- seq(320, 380, length = 4)
alpha_cluster3 <- alpha_cluster3[!alpha_cluster3 %in% alpha_cluster4]
alpha_cluster5 <- c(211)
alpha_cluster6 <- quad_1[!quad_1 %in% c(alpha_cluster3, alpha_cluster4, alpha_cluster5)]
alpha_cluster7 <- c(305)
alpha_cluster8 <- rep(6:9, times = 4) + 20 * rep(5:8, each = 4)
alpha_cluster9 <- c(54)
alpha_cluster10 <- (1:400)[!(1:400) %in% c(alpha_cluster1, alpha_cluster2, alpha_cluster3,
                                           alpha_cluster4, alpha_cluster5, alpha_cluster6,
                                           alpha_cluster7, alpha_cluster8, alpha_cluster9)]

gamma_alpha <- list(alpha_cluster1, alpha_cluster2, alpha_cluster3, alpha_cluster4, alpha_cluster5,
                   alpha_cluster6, alpha_cluster7, alpha_cluster8, alpha_cluster9, alpha_cluster10)

beta_cluster1 <- c(165, 166, 184:187, 204:207, 225, 226)
beta_cluster2 <- c(quad_2, quad_3)
beta_cluster2 <- beta_cluster2[!beta_cluster2 %in% beta_cluster1]
beta_cluster3 <- quad_4
beta_cluster4 <- quad_1

gamma_beta <- list(beta_cluster1, beta_cluster2, beta_cluster3, beta_cluster4)

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

for(sim_number in 1:tot_sim){

  ## these values were tested with simulation2
  alphabar <- 3.5 + c(-1, 1, 0.5, 0, -1, 
                      -0.5, -0.5, 0.5, 1, 0) * Delta_alpha[sim_number]
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
