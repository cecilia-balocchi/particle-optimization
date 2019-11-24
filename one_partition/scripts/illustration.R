library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/get_hyperparameters.R")
sourceCpp("src/ensm_cluster_mean.cpp")
sourceCpp("src/map_cluster.cpp")
sourceCpp("src/partition_summary.cpp")
sourceCpp("src/kmeans.cpp")

load("data/alphas.RData")
load("data/partitions.RData")


N <- 400
T <- 12

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

batch_size <- 5
batch <- 1
r <- 1

km_methods <- paste0("km", 1:10)
methods <- c("mle", "map", "lam_1", "lam_10", "lam_100", "fixed_0", "fixed_1", "fixed_n")
metrics <- c("MSER", "DLP", "K", "L", "Ptrue", "TIME")

for(sim_number in 1:6){
  alpha <- get(paste0("alpha_", sim_number))
  if(!file.exists(paste0("results/illustration_sim", sim_number, ".RData"))){
    print(paste("  Starting simulation", sim_number, "at", Sys.time()))
    results <- matrix(nrow = length(methods) + 10, ncol = length(metrics), dimnames = list(c(methods, km_methods), metrics))
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
      print(paste("    Found lambda = 1 at", Sys.time()))
      lam_10 <- ensm_cluster_mean(Y = Y, A_block = W, init_id = 2, 
                                  a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 10)
      print(paste("    Found lambda = 10 at", Sys.time()))
      lam_100 <- ensm_cluster_mean(Y = Y, A_block = W, init_id = 2,
                                   a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)
      print(paste("    Found lambda = 100 at", Sys.time()))
      fixed_0 <- partition_summary(Y, W, gamma_init = gamma_0, 
                                   a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
      print(paste("    Found fixed_0 at", Sys.time()))
      fixed_1 <- partition_summary(Y, W, gamma_init = gamma_1,
                                   a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
      print(paste("    Found fixed_1 at", Sys.time()))
      fixed_n <- partition_summary(Y, W, gamma_init = gamma_n, 
                                   a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
      print(paste("    Found fixed_n at", Sys.time()))
      km <- kmeans_particle(Y, W, gamma_1, max_split = 10, 
                            a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
      print(paste("    Found km at", Sys.time()))
      map <- map_partition(Y = Y, A_block = W, gamma_init = gamma_1, 
                           a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)
      print(paste("    Found MAP at", Sys.time()))
      ### Get the performance results
      results["mle", "MSER"] <- mean( (alpha - rowMeans(Y))^2)/mean( (alpha - fixed_0$alpha)^2)
      for(method in methods[!methods == "mle"]){
        fit <- get(method)
        results[method, "MSER"] <- mean( (alpha - fit$alpha)^2)/mean( (alpha - fixed_0$alpha)^2)
        results[method, "DLP"] <- mean( fit$log_post - fixed_0$log_post)
        results[method, "K"] <- mean(sapply(fit$particles, FUN = length))
        results[method, "L"] <- length(fit$particles)
        results[method, "Ptrue"] <- 1*any(sapply(fit$particles, FUN = partition_equal, gamma2 = gamma_0, A_block = W))
        if(method %in% c("lam_1", "lam_10", "lam_100", "map", "km")) results[method, "TIME"] <- fit$time
      }
      
      for(k in 1:10){
        results[paste0("km",k), "MSER"] <- mean( (alpha - km$alpha[,k])^2) / mean( (alpha - fixed_0$alpha)^2)
        results[paste0("km",k), "DLP"] <- km$log_post[k] - fixed_0$log_post
        results[paste0("km",k), "K"] <- length(km$particles[[k]])
        results[paste0("km",k), "Ptrue"] <- 1*partition_equal(km$particles[[k]], gamma2 = gamma_0, A_block = W)
      }
      
      assign(paste0("lam1_sim", sim_number), lam_1)
      assign(paste0("lam10_sim", sim_number), lam_10)
      assign(paste0("lam100_sim", sim_number), lam_100)
      assign(paste0("fixed_0_sim", sim_number), fixed_0)
      assign(paste0("fixed_1_sim", sim_number), fixed_1)
      assign(paste0("fixed_n_sim", sim_number), fixed_n)
      assign(paste0("km_sim", sim_number), km)
      assign(paste0("map_sim", sim_number), map)
      assign(paste0("results_sim", sim_number), results)
      save(list = paste0(c("lam1", "lam10", "lam100", "fixed_0", "fixed_1", "fixed_n",
                           "km", "map", "results"), "_sim", sim_number),
           file = paste0("results/illustration_sim", sim_number, ".RData"))
    }
  }

}
