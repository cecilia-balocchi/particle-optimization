# MAP simulations for simulation settings 4 -- 6
library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/get_hyperparameters.R")
source("scripts/get_tot_ss.R")

sourceCpp("src/map_cluster.cpp")


load("data/alphas.RData")
load("data/partitions.RData")


N <- 400
T <- 12

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
batch <- as.numeric(args[2])

if(sim_number != 6){
  alpha <- get(paste0("alpha_", sim_number))
} else{
  alpha <- matrix(0, nrow = nrow(alpha_1), ncol = ncol(alpha_1))
}

batch_size <- 2
fixed_metrics <- c("RMSE", "RAND", "K", "LP")
map_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics) + 1, dimnames = list(c(), c(fixed_metrics, "TIME")))

for(r in 1:batch_size){
  print(paste("  Starting r = ", r, "at", Sys.time()))
  set.seed(seed_seq[sim_number] * 100 + 10 * (batch - 1) + r)
  Y <- matrix(nrow = N, ncol = T)
  for(i in 1:N) Y[i,] <- rnorm(T, mean = alpha[i, r + (batch-1)*batch_size], sd = 1)
  hyper_params <- get_hyperparameters(Y, floor(log(N)), rho = 0.9)
  if(!is.null(hyper_params)){
    a1 <- hyper_params[["a1"]]
    a2 <- hyper_params[["a2"]]
    nu_sigma <- hyper_params[["nu_sigma"]]
    lambda_sigma <- hyper_params[["lambda_sigma"]]
    map <- map_partition(Y = Y, A_block = W,
                         a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
    
    map_results[r, "RMSE"] <- sqrt( mean((alpha - map$alpha)^2 ))
    map_results[r, "RAND"] <- rand_index(N, map$particles[[1]], gamma_0)
    map_results[r, "K"] <- length(map$particles[[1]])
    map_results[r, "LP"] <- map$log_post[1]
    map_results[r, "TIME"] <- map$time
  }
}

assign(paste0("map_results_sim", sim_number, "_", batch), map_results)
save(list = paste0("map_results_sim", sim_number, "_", batch),
     file = paste0("results/sim", sim_number, "_map_", batch, ".RData"))
