# Fixed partition methods
library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/get_hyperparameters.R")
source("scripts/get_tot_ss.R")

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

if(sim_number != 6){
  alpha <- get(paste0("alpha_", sim_number))
} else{
  alpha <- matrix(0, nrow = nrow(alpha_1), ncol = ncol(alpha_1))
}

batch_size <- 2
fixed_metrics <- c("RMSE", "RAND", "K", "LP")

km_metrics <- c(paste0("RMSE_", 1:10), paste0("K_", 1:10), paste0("RAND_", 1:10), paste0("LP_", 1:10), "TIME")
km_results <- matrix(nrow = batch_size, ncol = length(km_metrics), dimnames = list(c(), km_metrics))
km_ss <- matrix(nrow = batch_size, ncol = 10)

sc_metrics <- c(paste0("RMSE_", 1:10), paste0("K_", 1:10), paste0("RAND_", 1:10), paste0("LP_", 1:10), "TIME")
sc_results <- matrix(nrow = batch_size, ncol = length(sc_metrics), dimnames = list(c(), sc_metrics))
sc_ss <- matrix(nrow = batch_size, ncol = 10)

fixed_metrics <- c("RMSE", "RAND", "K", "LP")
fixed_0_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_1_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_n_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))

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

    for(l in 1:10){
      sc_results[r, paste0("RMSE_",l)] <- sqrt(mean( (alpha - sc$alpha[,l])^2))
      sc_results[r, paste0("K_",l)] <- length(sc$particles[[l]])
      sc_results[r, paste0("RAND_",l)] <- rand_index(N, sc$particles[[l]], gamma_0)
      sc_results[r, paste0("LP_",l)] <- sc$log_post[l]
      sc_ss[r, l] <- get_tot_ss(rowMeans(Y), sc$particles[[l]])
    }
    sc_results[r, "TIME"] <- sc$time
    print(paste("    Finished sc at", Sys.time()))

    
  } # closes if checking that the hyper-parameters were not null
}

assign(paste0("fixed_0_results_sim", sim_number, "_", batch), fixed_0_results)
assign(paste0("fixed_1_results_sim", sim_number, "_", batch), fixed_1_results)
assign(paste0("fixed_n_results_sim", sim_number, "_", batch), fixed_n_results)
assign(paste0("km_ss_sim", sim_number, "_", batch), km_ss)
assign(paste0("km_results_sim", sim_number, "_", batch), km_results)
assign(paste0("sc_ss_sim", sim_number, "_", batch), sc_ss)
assign(paste0("sc_results_sim", sim_number, "_", batch), sc_results)

save_list <- paste0(c("fixed_0_", "fixed_1_", "fixed_n_","km_", "sc_"), "results_sim", sim_number, "_", batch)
save_list <- c(save_list, paste0("km_ss_sim", sim_number, "_", batch), paste0("sc_ss_sim", sim_number, "_", batch))
save(list = save_list, file = paste0("results/sim", sim_number, "_fixed_results_", batch, ".RData"))
