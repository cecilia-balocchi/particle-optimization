rm(list = ls())
library(igraph) # to generate the graph
library(MASS)

###############################
### Generate the data first ###
###############################

# Generate a 10x10 grid graph
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

# Define the clusters
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

beta_cluster1 <- c(22,23,32,33,43,44,53,54,62,63,72,73)
beta_cluster2 <- quad_2
beta_cluster3 <- quad_4
beta_cluster4 <- (1:100)[!(1:100) %in% c(beta_cluster1,beta_cluster2,beta_cluster3)]
gamma_beta <- list(beta_cluster1, beta_cluster2,beta_cluster3,beta_cluster4)

zA <- numeric(N)
for(i in 1:length(gamma_alpha)){
  zA[ gamma_alpha[[i]] ] <- i
}
zB <- numeric(N)
for(i in 1:length(gamma_beta)){
  zB[ gamma_beta[[i]] ] <- i
}

# If you want to plot, you can use the function "plot_clusters_grid" in the Github repo in
# particle-optimization/two_partitions/scripts/plot_functions.R
# plot_clusters_grid(gamma_alpha, W)
# plot_clusters_grid(gamma_beta, W)

# Generate data (alpha and beta)
# (this generates more datasets than necessary)
tot_sim <- 3  # how many configurations ( = length(Delta_alpha)), i.e. high/medium/low cluster separation

seed_seq <- c(129, 724, 603, 212, 1123, 2391) 
n_sim <- 100
Delta_alpha <- c(2, 1, 0)
Delta_beta <- c(1, 0.5, 0)
rho <- 0.95
a1 <- 1/8 
b1 <- 1/8 
sigma <- 0.25
sigma2 <- sigma^2  # Note: in practice, within cluster sd will be 0.25/sqrt(8) =~ 0.09

for(sim_number in 1:tot_sim){
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

# Generate data (X and Y)
sim_number <- 1  # configuration of cluster separation (high/medium/low)
ALPHA <- get(paste0("alpha_", sim_number))
BETA <- get(paste0("beta_", sim_number))
r = 2  # choose one of the simulated datasets
alpha <- ALPHA[, r]
beta <- BETA[, r]

t <- 12
set.seed(seed_seq[sim_number] * 100 + r)
X <- matrix(rnorm(N*t, sd = 2), ncol = t, nrow = N)
Y <- matrix(nrow = N, ncol = t)
for(i in 1:N) Y[i,] <- rnorm(t, mean = alpha[i] + beta[i] * X[i,], 
                             sd = sqrt(sigma2))

##############################
### Run the PARTOPT method ###
##############################
library(Rcpp)
library(RcppArmadillo)
library(mcclust) # useful for metrics between partitions
library(PARTOPT)

L <- 5
lam_10 <- PARTOPT(Y, X, W, L = L)

lam_output <- lam_10$final
alpha_wtd <- rep(0, times = N)
for(ll in 1:L)  alpha_wtd <- alpha_wtd + lam_output$alpha_particle[, ll ] * lam_output$w[ ll ] 
beta_wtd <- rep(0, times = N)
for(ll in 1:L)  beta_wtd <- beta_wtd + lam_output$beta_particle[, ll ] * lam_output$w[ ll ] 

# comparing with the truth
plot(alpha, alpha_wtd); abline(0,1)
plot(beta, beta_wtd); abline(0,1)

# best particle
j_max <- which.max(lam_output$w)
gammaA_top <- lam_output$particle_set_A[[j_max]]
gammaB_top <- lam_output$particle_set_B[[j_max]]

zA_top <- numeric(N)
for(i in 1:length(gammaA_top)){
  zA_top[ gammaA_top[[i]] ] <- i
}
zB_top <- numeric(N)
for(i in 1:length(gammaB_top)){
  zB_top[ gammaB_top[[i]] ] <- i
}

mcclust::arandi(cl1 = zA, cl2 = zA_top)
mcclust::arandi(cl1 = zB, cl2 = zB_top)

# we could re-run the function using the already discovered hyperparameters
h1 <- lam_10$hyperpar[1:4]
h2 <- lam_10$hyperpar[5:6]

lam_10 <- PARTOPT(Y, X, W, L = 10, hyperpar = h1, hyperpar_sigma = h2)

# we can provide a pair of partitions as a first guess
lam_10 <- PARTOPT(Y, X, W, L = 3, gamma_init_A = gammaA_top, gamma_init_B = gammaB_top)

# we can change the value of lambda (larger values will ensure diversity among the particles)
lam_10 <- PARTOPT(Y, X, W, L = 3, lambda = 30)
