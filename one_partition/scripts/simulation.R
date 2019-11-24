library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/get_hyperparameters.R")
source("scripts/get_tot_ss.R")
sourceCpp("src/ensm_cluster_mean.cpp")
sourceCpp("src/map_cluster.cpp")
sourceCpp("src/partition_summary.cpp")
sourceCpp("src/kmeans.cpp")
sourceCpp("src/spectral_clustering.cpp")
load("data/alphas.RData")
load("data/partitions.RData")


N <- 400
T <- 12

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
batch <- as.numeric(args[2])

print(sim_number)
alpha <- get(paste0("alpha_", sim_number))

batch_size <- 4


particle_metrics <- c("RMSE", "L", "RMSE_top", "RAND_top", "K_top", "LP_top", paste0("RMSE_", 1:10), paste0("wtd_RMSE_", 1:10), paste0("LP_", 1:10), 
                      paste0("K_", 1:10), paste0("avg_K_", 1:10), paste0("RAND_", 1:10), paste0("avg_RAND_", 1:10),"TIME")
lam_1_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
lam_10_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
lam_100_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))


# orig_K_top: the partition with the best log-posterior was found when we split into orig_K_top clusters
# RMSE_top, RAND_top, K_top: RMSE, rand index, and number of connected clusters in the partition with the best log posterior
# RMSE_L, RAND_L, K_L: RMSE, rand index, and number of connected clusters for the 10 we find
km_metrics <- c(paste0("RMSE_", 1:10), paste0("K_", 1:10), paste0("RAND_", 1:10), paste0("LP_", 1:10), "TIME")
km_results <- matrix(nrow = batch_size, ncol = length(km_metrics), dimnames = list(c(), km_metrics))
km_ss <- matrix(nrow = batch_size, ncol = 10)

sc_metrics <- c(paste0("RMSE_", 1:10), paste0("K_", 1:10), paste0("RAND_", 1:10), paste0("LP_", 1:10), "TIME")
sc_results <- matrix(nrow = batch_size, ncol = length(sc_metrics), dimnames = list(c(), sc_metrics))
sc_ss <- matrix(nrow = batch_size, ncol  = 10)

fixed_metrics <- c("RMSE", "RAND", "K", "LP")
fixed_0_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_1_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_n_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
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
    
    lam_1 <- ensm_cluster_mean(Y = Y, A_block = W, init_id = 2, 
                               a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 1)
    
    lam_1_results[r, "RMSE"] <- sqrt(mean( (alpha - lam_1$alpha)^2))
    lam_1_results[r, "L"] <- length(lam_1$particles)
    tmp_K <- sapply(lam_1$particles, FUN = length)
    tmp_rand <- sapply(lam_1$particles, FUN = rand_index, N = N, gamma1 = gamma_0)
    lam_1_results[r, "RMSE_top"] <- sqrt(mean( (alpha - lam_1$alpha_hat_particle[,1])^2))
    lam_1_results[r, "K_top"] <- tmp_K[1]
    lam_1_results[r, "LP_top"] <- lam_1$log_post[1]
    lam_1_results[r, "RAND_top"] <- tmp_rand[1]
    for(l in 1:length(lam_1$particles)){
      lam_1_results[r, paste0("RMSE_",l)] <- sqrt( mean( (alpha - lam_1$alpha_hat_particle[,l])^2))
      tmp_probs <- lam_1$pstar[1:l]/sum(lam_1$pstar[1:l])
      tmp_alpha <- rep(0, times = N)
      for(ll in 1:l)  tmp_alpha <- tmp_alpha + lam_1$alpha_hat_particle[,ll] * tmp_probs[ll] 
      lam_1_results[r, paste0("wtd_RMSE_",l)] <- sqrt( mean( (alpha - tmp_alpha)^2))
      lam_1_results[r, paste0("K_", l)] <- tmp_K[l]
      lam_1_results[r, paste0("avg_K_", l)] <- mean(tmp_K[1:l])
      lam_1_results[r, paste0("RAND_", l)] <- tmp_rand[l]
      lam_1_results[r, paste0("avg_RAND_",l)] <- mean(tmp_rand[1:l])
      lam_1_results[r, paste0("LP_", l)] <- lam_1$log_post[l]
    }
    lam_1_results[r, "TIME"] <- lam_1$time
    print(paste("    Found lambda = 1 at", Sys.time()))
    
    lam_10 <- ensm_cluster_mean(Y = Y, A_block = W, init_id = 2, 
                                a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 10)
    
    lam_10_results[r, "RMSE"] <- sqrt(mean( (alpha - lam_10$alpha)^2))
    lam_10_results[r, "L"] <- length(lam_10$particles)
    tmp_K <- sapply(lam_10$particles, FUN = length)
    tmp_rand <- sapply(lam_10$particles, FUN = rand_index, N = N, gamma1 = gamma_0)
    lam_10_results[r, "RMSE_top"] <- sqrt(mean( (alpha - lam_10$alpha_hat_particle[,1])^2 ))
    lam_10_results[r, "K_top"] <- tmp_K[1]
    lam_10_results[r, "LP_top"] <- lam_10$log_post[1]
    lam_10_results[r, "RAND_top"] <- tmp_rand[1]
    for(l in 1:length(lam_10$particles)){
      lam_10_results[r, paste0("RMSE_",l)] <- sqrt( mean( (alpha - lam_10$alpha_hat_particle[,l])^2))
      tmp_probs <- lam_10$pstar[1:l]/sum(lam_10$pstar[1:l])
      tmp_alpha <- rep(0, times = N)
      for(ll in 1:l)  tmp_alpha <- tmp_alpha + lam_10$alpha_hat_particle[,ll] * tmp_probs[ll] 
      lam_10_results[r, paste0("wtd_RMSE_",l)] <- sqrt( mean( (alpha - tmp_alpha)^2 ))
      lam_10_results[r, paste0("K_", l)] <- tmp_K[l]
      lam_10_results[r, paste0("avg_K_", l)] <- mean(tmp_K[1:l])
      lam_10_results[r, paste0("RAND_", l)] <- tmp_rand[l]
      lam_10_results[r, paste0("avg_RAND_",l)] <- mean(tmp_rand[1:l])
      lam_10_results[r, paste0("LP_", l)] <- lam_10$log_post[l]
    }
    lam_10_results[r, "TIME"] <- lam_10$time
    print(paste("    Finished lambda = 10 at", Sys.time()))
    
    
    lam_100 <- ensm_cluster_mean(Y = Y, A_block = W, init_id = 2,
                                 a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)
    lam_100_results[r, "RMSE"] <- sqrt(mean( (alpha - lam_100$alpha)^2 ))
    lam_100_results[r, "L"] <- length(lam_100$particles)
    tmp_K <- sapply(lam_100$particles, FUN = length)
    tmp_rand <- sapply(lam_100$particles, FUN = rand_index, N = N, gamma1 = gamma_0)
    lam_100_results[r, "RMSE_top"] <- sqrt(mean( (alpha - lam_100$alpha_hat_particle[,1])^2 ))
    lam_100_results[r, "K_top"] <- tmp_K[1]
    lam_100_results[r, "LP_top"] <- lam_100$log_post[1]
    lam_100_results[r, "RAND_top"] <- tmp_rand[1]
    for(l in 1:length(lam_100$particles)){
      lam_100_results[r, paste0("RMSE_",l)] <- sqrt( mean( (alpha - lam_100$alpha_hat_particle[,l])^2 ))
      tmp_probs <- lam_100$pstar[1:l]/sum(lam_100$pstar[1:l])
      tmp_alpha <- rep(0, times = N)
      for(ll in 1:l)  tmp_alpha <- tmp_alpha + lam_100$alpha_hat_particle[,ll] * tmp_probs[ll] 
      lam_100_results[r, paste0("wtd_RMSE_",l)] <- sqrt( mean( (alpha - tmp_alpha)^2 ))
      lam_100_results[r, paste0("K_", l)] <- tmp_K[l]
      lam_100_results[r, paste0("avg_K_", l)] <- mean(tmp_K[1:l])
      lam_100_results[r, paste0("RAND_", l)] <- tmp_rand[l]
      lam_100_results[r, paste0("avg_RAND_",l)] <- mean(tmp_rand[1:l])
      lam_100_results[r, paste0("LP_", l)] <- lam_100$log_post[l]
    }
    lam_100_results[r, "TIME"] <- lam_100$time    
    print(paste("    Finished lambda = 100 at", Sys.time()))
    
    
    fixed_0 <- partition_summary(Y, W, gamma_init = gamma_0, 
                                 a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
    
    fixed_0_results[r,"RMSE"] <- sqrt( mean( (alpha - fixed_0$alpha)^2 ))
    fixed_0_results[r, "RAND"] <- rand_index(N, fixed_0$particles[[1]], gamma_0)
    fixed_0_results[r, "K"] <- length(fixed_0$particles)
    fixed_0_results[r, "LP"] <- fixed_0$log_post
    print(paste("    Finished fixed_0 at", Sys.time()))
    
    fixed_1 <- partition_summary(Y, W, gamma_init = gamma_1,
                                 a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
    fixed_1_results[r,"RMSE"] <- sqrt( mean( (alpha - fixed_1$alpha)^2 ))
    fixed_1_results[r, "RAND"] <- rand_index(N, fixed_1$particles[[1]], gamma_0)
    fixed_1_results[r, "K"] <- length(fixed_1$particles)
    fixed_1_results[r, "LP"] <- fixed_1$log_post
    print(paste("    Finished fixed_1 at", Sys.time()))
    
    fixed_n <- partition_summary(Y, W, gamma_init = gamma_n, 
                                 a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
    fixed_n_results[r,"RMSE"] <- sqrt( mean( (alpha - fixed_n$alpha)^2 ))
    fixed_n_results[r, "RAND"] <- rand_index(N, fixed_n$particles[[1]], gamma_0)
    fixed_n_results[r, "K"] <- length(fixed_n$particles)
    fixed_n_results[r, "LP"] <- fixed_n$log_post
    print(paste("    Finished fixed_n at", Sys.time()))
    
    
    km <- kmeans_particle(Y, W, max_split = 10, 
                          a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
    
    #top_ix <- which.max(km$log_post)
    #km_results[r, "RMSE_top"] <- sqrt(mean( (alpha - km$alpha[,top_ix])^2 ))
    #km_results[r, "RAND_top"] <- rand_index(N, km$particles[[top_ix]], gamma_0)
    #km_results[r, "K_top"] <- length(km$particles[[top_ix]])
    #km_results[r, "LP_top"] <- km$log_post[top_ix[]]
    #km_results[r, "orig_K_top"] <- top_ix
    for(l in 1:10){
      km_results[r, paste0("RMSE_",l)] <- sqrt(mean( (alpha - km$alpha[,l])^2 ))
      km_results[r, paste0("K_",l)] <- length(km$particles[[l]])
      km_results[r, paste0("RAND_",l)] <- rand_index(N, km$particles[[l]], gamma_0)
      km_results[r, paste0("LP_",l)] <- km$log_post[l]
      km_ss[r, l] <- get_tot_ss(rowMeans(Y),km$particles[[l]])
    }
    km_results[r, "TIME"] <- km$time
    print(paste("    Finished km at", Sys.time()))
    
    sc <- spectral_particle(Y, W, max_splits = 10, 
                            a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
    
    #top_ix <- which.max(sc$log_post)
    #sc_results[r, "RMSE_top"] <- sqrt(mean( (alpha - sc$alpha[,top_ix])^2))
    #sc_results[r, "RAND_top"] <- rand_index(N, sc$particles[[top_ix]], gamma_0)
    #sc_results[r, "K_top"] <- length(sc$particles[[top_ix]])
    #sc_results[r, "LP_top"] <- sc$log_post[top_ix[]]
    #sc_results[r, "orig_K_top"] <- top_ix
    for(l in 1:10){
      sc_results[r, paste0("RMSE_",l)] <- sqrt(mean( (alpha - sc$alpha[,l])^2))
      sc_results[r, paste0("K_",l)] <- length(sc$particles[[l]])
      sc_results[r, paste0("RAND_",l)] <- rand_index(N, sc$particles[[l]], gamma_0)
      sc_results[r, paste0("LP_",l)] <- sc$log_post[l]
      sc_ss[r, l] <- get_tot_ss(rowMeans(Y), sc$particles[[l]])
    }
    sc_results[r, "TIME"] <- sc$time
    print(paste("    Finished sc at", Sys.time()))
    
    map <- map_partition(Y = Y, A_block = W,
                         a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
    
    map_results[r, "RMSE"] <- sqrt( mean((alpha - map$alpha)^2 ))
    map_results[r, "RAND"] <- rand_index(N, map$particles[[1]], gamma_0)
    map_results[r, "K"] <- length(map$particles[[1]])
    map_results[r, "LP"] <- map$log_post[1]
    map_results[r, "TIME"] <- map$time
    print(paste("    Finished MAP at", Sys.time()))

  } # closes if checking that the hyper-parameters were not null
}

assign(paste0("lam_1_results_sim", sim_number, "_", batch), lam_1_results)
assign(paste0("lam_10_results_sim", sim_number, "_", batch), lam_10_results)
assign(paste0("lam_100_results_sim", sim_number, "_", batch), lam_100_results)
assign(paste0("fixed_0_results_sim", sim_number, "_", batch), fixed_0_results)
assign(paste0("fixed_1_results_sim", sim_number, "_", batch), fixed_1_results)
assign(paste0("fixed_n_results_sim", sim_number, "_", batch), fixed_n_results)
assign(paste0("km_ss_sim", sim_number, "_", batch), km_ss)
assign(paste0("km_results_sim", sim_number, "_",batch), km_results)
assign(paste0("sc_ss_sim", sim_number, "_", batch), sc_ss)
assign(paste0("sc_results_sim", sim_number, "_", batch), sc_results)

assign(paste0("map_results_sim", sim_number, "_", batch), map_results)

save_list <- paste0(c("lam_1_", "lam_10_", "lam_100_", 
                      "fixed_0_", "fixed_1_", "fixed_n_",
                      "km_", "sc_", "map_"), "results_sim", sim_number, "_", batch)
save_list <- c(save_list, paste0("km_ss_sim", sim_number, "_",batch), paste0("sc_ss_sim", sim_number, "_", batch))
save(list = save_list, file = paste0("results/sim", sim_number, "_", batch, ".RData"))
