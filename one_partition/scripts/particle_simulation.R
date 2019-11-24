# Simulation for our method
library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/get_hyperparameters.R")
source("scripts/get_tot_ss.R")

sourceCpp("src/ensm_cluster_mean.cpp")

load("data/alphas.RData")
load("data/partitions.RData")


N <- 400
T <- 12

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
lambda <- as.numeric(args[2])
batch <- as.numeric(args[3])

if(sim_number != 6){
   alpha <- get(paste0("alpha_", sim_number))
} else {
  alpha <- matrix(0, nrow = nrow(alpha_1), ncol = ncol(alpha_1))
}

batch_size <- 2


particle_metrics <- c("RMSE", "L", "RMSE_top", "RAND_top", "K_top", "LP_top", paste0("RMSE_", 1:10), paste0("wtd_RMSE_", 1:10), paste0("LP_", 1:10), 
                      paste0("K_", 1:10), paste0("avg_K_", 1:10), paste0("RAND_", 1:10), paste0("avg_RAND_", 1:10),
                      paste0("LP_", 1:10), "TIME")
results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))

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
    
    fit <- ensm_cluster_mean(Y = Y, A_block = W, init_id = 2, 
                               a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = lambda)
    
    results[r, "RMSE"] <- sqrt(mean( (alpha - fit$alpha)^2))
    results[r, "L"] <- length(fit$particles)
    tmp_K <- sapply(fit$particles, FUN = length)
    tmp_rand <- sapply(fit$particles, FUN = rand_index, N = N, gamma1 = gamma_0)
    results[r, "RMSE_top"] <- sqrt(mean( (alpha - fit$alpha_hat_particle[,1])^2))
    results[r, "K_top"] <- tmp_K[1]
    results[r, "LP_top"] <- fit$log_post[1]
    results[r, "RAND_top"] <- tmp_rand[1]
    for(l in 1:length(fit$particles)){
      results[r, paste0("RMSE_",l)] <- sqrt( mean( (alpha - fit$alpha_hat_particle[,l])^2))
      tmp_probs <- fit$pstar[1:l]/sum(fit$pstar[1:l])
      tmp_alpha <- rep(0, times = N)
      for(ll in 1:l)  tmp_alpha <- tmp_alpha + fit$alpha_hat_particle[,ll] * tmp_probs[ll] 
      results[r, paste0("wtd_RMSE_",l)] <- sqrt( mean( (alpha - tmp_alpha)^2))
      results[r, paste0("K_", l)] <- tmp_K[l]
      results[r, paste0("avg_K_", l)] <- mean(tmp_K[1:l])
      results[r, paste0("RAND_", l)] <- tmp_rand[l]
      results[r, paste0("avg_RAND_",l)] <- mean(tmp_rand[1:l])
      results[r, paste0("LP_", l)] <- fit$log_post[l]
    }
    results[r, "TIME"] <- fit$time
    print(paste("    Finished", r, "at", Sys.time()))
    
  }
}
assign(paste0("lam_", lambda, "_results_sim", sim_number, "_", batch), results)
save(list = paste0("lam_", lambda, "_results_sim", sim_number, "_", batch),
     file = paste0("results/sim", sim_number, "_lam_", lambda, "_", batch, ".RData"))


