library(Rcpp)
library(RcppArmadillo)
library(mcclust)
library(cluster)

source("src/partition_functions.R")
source("src/ensm_hyperpar.R")
source("src/ensm_givenhyperpar.R")
source("src/summary_givenhyperpar.R")
source("src/load_anderson2017.R")

## need to work on the data generation process
load("data/alphabetas.RData")
load("data/partitions.RData")

N <- 400
T <- 12

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
batch <- as.numeric(args[2])
## if batch_size = 4, and n_sim = 100, then we need batch to range from 1 to 25

print(sim_number)
ALPHA <- get(paste0("alpha_", sim_number))
BETA <- get(paste0("beta_", sim_number))
KMsplits <- 5
KMsplits_tot <- KMsplits * KMsplits
SCsplits <- 10
SCsplits_tot <- SCsplits * SCsplits

batch_size <- 4

particle_metrics <- c("RAND_A_top", "RAND_B_top","RANDadj_A_top", "RANDadj_B_top", "VI_A_top", "VI_B_top", "LP_top", 
                      "avg_RAND_A", "avg_RAND_B","avg_RANDadj_A", "avg_RANDadj_B","avg_VI_A", "avg_VI_B",
                      "wavg_RAND_A", "wavg_RAND_B","wavg_RANDadj_A", "wavg_RANDadj_B","wavg_VI_A", "wavg_VI_B",
                      "RMSE","RMSE_top", "pred_err", "pred_err_top",
                      "K_A_top","K_B_top","avg_K_A","avg_K_B","wavg_K_A","wavg_K_B",
                      "p_trueA","p_trueB","p_trueAB",
                      "L", "TIME")
lam_1_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
lam_10_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
lam_100_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))

# orig_K_top: the partition with the best log-posterior was found when we split into orig_K_top clusters
# RMSE_top, RAND_top, K_top: RMSE, rand index, and number of connected clusters in the partition with the best log posterior
# RMSE_L, RAND_L, K_L: RMSE, rand index, and number of connected clusters for the 10 we find

km_sc_metrics <- c("RMSE","RMSE_mle", "K_A", "K_B",
                "RAND_A", "RAND_B","RANDadj_A", "RANDadj_B", "VI_A", "VI_B", "pred_err","LP","TIME")
km_metrics <- km_sc_metrics
sc_metrics <- km_sc_metrics

km_results <- matrix(nrow = batch_size, ncol = length(km_metrics), dimnames = list(c(), km_metrics))
sc_results <- matrix(nrow = batch_size, ncol = length(sc_metrics), dimnames = list(c(), sc_metrics))

fixed_metrics <- c("RMSE", "RAND_A","RAND_B","RANDadj_A", "RANDadj_B", "VI_A", "VI_B", "K_A","K_B","pred_err", "LP")
fixed_0_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_1_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_n_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
map_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics) + 1, dimnames = list(c(), c(fixed_metrics, "TIME")))

anderson_metrics <- c("RMSE", 
                      "post_RAND_A","post_RAND_B", "post_RANDadj_A","post_RANDadj_B", "post_VI_A","post_VI_B", "post_K_A","post_K_B",
                      "post_mod_RAND_A","post_mod_RAND_B", "post_mod_RANDadj_A","post_mod_RANDadj_B", "post_mod_VI_A","post_mod_VI_B", "post_mod_K_A","post_mod_K_B",
                      "est_RAND_A","est_RAND_B", "est_RANDadj_A","est_RANDadj_B", "est_VI_A","est_VI_B", "est_K_A","est_K_B",
                      "est_mod_RAND_A","est_mod_RAND_B", "est_mod_RANDadj_A","est_mod_RANDadj_B", "est_mod_VI_A","est_mod_VI_B", "est_mod_K_A","est_mod_K_B",
                      "pred_err", "LP", "TIME")
anderson_ks <- c("33","35","53","55")
anderson_results <- array(NA, dim = c(batch_size,length(anderson_metrics),length(anderson_ks)), dimnames = list(c(), anderson_metrics, anderson_ks))

ihs_transform <- function(x){
  return( log(x + sqrt(x^2 + 1)) - log(2))
}
# x = seq(0,10, length = 100)
# y = ihs_transform(x)
# all.equal(sinh(y + log(2)),x) ## TRUE

t <- 12
t_new <- 13
tmp_time <- 1:t - mean(1:t)
t_new_scaled <- (t_new - mean(1:t))/sd(tmp_time)
# t_new_scaled <- (t_new - mean(1:t))/sqrt(sum(tmp_time^2))

zA <- numeric(N)
for(i in 1:length(gamma_alpha)){
  zA[ gamma_alpha[[i]] ] <- i
}
zB <- numeric(N)
for(i in 1:length(gamma_beta)){
  zB[ gamma_beta[[i]] ] <- i
}

## todo:
miter <- 100
anderson_ndraws <- 1000
L <- 10

# miter <- 1
# anderson_ndraws <- 100
# L <- 2

for(r in 1:batch_size){
  print(paste("  Starting r = ", r, "at", Sys.time()))
  set.seed(seed_seq[sim_number] * 100 + 10 * (batch - 1) + r)
  alpha <- ALPHA[, r + (batch-1)*batch_size]
  beta <- BETA[, r + (batch-1)*batch_size]
  Y <- matrix(nrow = N, ncol = T)
  for(i in 1:N) Y[i,] <- rnorm(T, mean = alpha[i] + beta[i] * X[i,], 
                                  sd = sqrt(sigma2))
  C <- sinh(Y + log(2)) 
  C <- round(C) # this is the data for running Anderson code
  # Y_outsample <- matrix(nrow = N, ncol = 1)
  Y_outsample <- rnorm(N, mean = alpha + beta * t_new_scaled, sd = sqrt(sigma2))
  C_outsample <- sinh(Y_outsample + log(2)) 
  C_outsample <- round(C_outsample)

  if(any(C < 0)){
    warning("Some counts were negative: specifically in locations: ",which(C<0), " with values ", C[which(C<0)])
    C[which(C<0)] <- 0
  }
  if(any(C_outsample < 0)){
    warning("Some (outsample) counts were negative: specifically in locations: ",which(C_outsample<0), " with values ", C_outsample[which(C_outsample<0)])
    C_outsample[which(C_outsample<0)] <- 0
  }

  betas_mle <- numeric(N)
  for(i in 1:N)
    betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
  alphas_mle <- rowMeans(Y)
  
  ## let's run MAP first, and decide the hyper-parameters.
  z <- Sys.time()
  map <- ensm_hyperpar(Y, X, W, L = 1, max_iter = miter, 
                         eta_py = 1.0, sigma_py = 0.0, 
                         lambda = 1, rho = 0.9)
  time1 <- Sys.time() - z 
  map_results[r, "TIME"] <- time1

  ## let's set the hyper-parameters for the next ones:
  hyperpar <- map$hyperpar # a1, a2, b1, b2, alpha_sigma, beta_sigma

  ### If we wanted to compare with the "old hyperparameters" CHECK
  old_hyperpar <- c(map$first$hyper_parameters, hyperpar[5:6])

  map <- map$adjusted

  gammaA_map <- map$particle_set_A[[1]]
  gammaB_map <- map$particle_set_B[[1]]

  map_results[r, "RMSE"] <- sqrt(mean( c((alpha - map$alpha_particle[,1])^2,
                                         (beta - map$beta_particle[,1])^2) ))
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_map)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_map)

  map_results[r, "RAND_A"] <- tmprandA$rand
  map_results[r, "RAND_B"] <- tmprandB$rand
  map_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  map_results[r, "RANDadj_B"] <- tmprandB$rand_adj 
  map_results[r, "VI_A"] <- tmprandA$vi
  map_results[r, "VI_B"] <- tmprandB$vi
  map_results[r, "K_A"] <- length(gammaA_map)
  map_results[r, "K_B"] <- length(gammaB_map)
  map_results[r, "LP"] <- map$log_posterior[1]

  Ypred_top <- map$alpha_particle[,1] + map$beta_particle[,1] * t_new_scaled
  map_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred_top)^2 )) 

  rm(map)
  gc()

  print(paste("    Finished MAP at", Sys.time()))
  #########

  z <- Sys.time()
  lam_1 <- ensm_givenhyperpar(Y, X, W, L = L, max_iter = miter, 
                         hyperpar = hyperpar, gamma_init_A = gammaA_map, gamma_init_B = gammaB_map, 
                         eta_py = 1.0, sigma_py = 0.0, 
                         lambda = 1, rho = 0.9, A_or_B_first = 0.5)
  time1 <- Sys.time() - z 
  lam_1_results[r, "TIME"] <- time1
  
  lam_output <- lam_1$output

 
  
  alpha_wtd <- rep(0, times = N)
  for(ll in 1:L)  alpha_wtd <- alpha_wtd + lam_output$alpha_particle[, ll ] * lam_output$w[ ll ] 
  beta_wtd <- rep(0, times = N)
  for(ll in 1:L)  beta_wtd <- beta_wtd + lam_output$beta_particle[, ll ] * lam_output$w[ ll ] 
  
  Lstar <- length(lam_output$unik_index)
  lam_1_results[r, "RMSE"] <- sqrt(mean( c((alpha - alpha_wtd)^2, (beta - beta_wtd)^2) )) ## todo
  lam_1_results[r, "L"] <- Lstar
  tmp_K_A <- sapply(lam_output$particle_set_A, FUN = length)
  tmp_K_B <- sapply(lam_output$particle_set_B, FUN = length)
  
  
  j_max <- which.max(lam_output$w)
  gammaA_top <- lam_output$particle_set_A[[j_max]]
  gammaB_top <- lam_output$particle_set_B[[j_max]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_top)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_top)
  
  lam_1_results[r, "RMSE_top"] <- sqrt(mean( c((alpha - lam_output$alpha_particle[,j_max])^2,(beta - lam_output$beta_particle[,j_max])^2) ))
  lam_1_results[r, "K_A_top"] <- length(gammaA_top)
  lam_1_results[r, "K_B_top"] <- length(gammaB_top)
  lam_1_results[r, "LP_top"] <- lam_output$log_posterior[j_max]
  lam_1_results[r, "RAND_A_top"] <- tmprandA$rand
  lam_1_results[r, "RAND_B_top"] <- tmprandB$rand
  lam_1_results[r, "RANDadj_A_top"] <- tmprandA$rand_adj
  lam_1_results[r, "RANDadj_B_top"] <- tmprandB$rand_adj
  lam_1_results[r, "VI_A_top"] <- tmprandA$vi
  lam_1_results[r, "VI_B_top"] <- tmprandB$vi

  tmprandAs <- sapply(lam_output$particle_set_A, rand_adj_vi, N = 400, gamma1 = gamma_alpha)
  tmprandBs <- sapply(lam_output$particle_set_B, rand_adj_vi, N = 400, gamma1 = gamma_beta)
  tmprandAs <- matrix(as.numeric(tmprandAs), ncol = L)
  tmprandBs <- matrix(as.numeric(tmprandBs), ncol = L)
  tmprandAs_avg <- rowMeans(tmprandAs)
  tmprandBs_avg <- rowMeans(tmprandBs)
  tmprandAs_wavg <- as.numeric(tmprandAs %*% matrix(lam_output$w, ncol = 1))
  tmprandBs_wavg <- as.numeric(tmprandBs %*% matrix(lam_output$w, ncol = 1))

  lam_1_results[r, "avg_RAND_A"] <- tmprandAs_avg[1]
  lam_1_results[r, "avg_RAND_B"] <- tmprandBs_avg[1]
  lam_1_results[r, "avg_RANDadj_A"] <- tmprandAs_avg[2]
  lam_1_results[r, "avg_RANDadj_B"] <- tmprandBs_avg[2]
  lam_1_results[r, "avg_VI_A"] <- tmprandAs_avg[3]
  lam_1_results[r, "avg_VI_B"] <- tmprandBs_avg[3]

  lam_1_results[r, "wavg_RAND_A"] <- tmprandAs_wavg[1]
  lam_1_results[r, "wavg_RAND_B"] <- tmprandBs_wavg[1]
  lam_1_results[r, "wavg_RANDadj_A"] <- tmprandAs_wavg[2]
  lam_1_results[r, "wavg_RANDadj_B"] <- tmprandBs_wavg[2]
  lam_1_results[r, "wavg_VI_A"] <- tmprandAs_wavg[3]
  lam_1_results[r, "wavg_VI_B"] <- tmprandBs_wavg[3]

  lam_1_results[r, "p_trueA"] <- mean(tmprandAs[1,] == 1)
  lam_1_results[r, "p_trueB"] <- mean(tmprandBs[1,] == 1)
  lam_1_results[r, "p_trueAB"] <- mean((tmprandAs[1,] == 1) & (tmprandBs[1,] == 1))

  tmpKAs <- sapply(lam_output$particle_set_A, FUN = length)
  tmpKBs <- sapply(lam_output$particle_set_B, FUN = length)

  lam_1_results[r, "avg_K_A"] <- mean(tmpKAs)
  lam_1_results[r, "avg_K_B"] <- mean(tmpKBs)
  lam_1_results[r, "wavg_K_A"] <- sum(tmpKAs * lam_output$w)
  lam_1_results[r, "wavg_K_B"] <- sum(tmpKBs * lam_output$w)
   
  Ypred_wtd <- alpha_wtd + beta_wtd * t_new_scaled
  Ypred_top <- lam_output$alpha_particle[,j_max] + lam_output$beta_particle[,j_max] * t_new_scaled
  lam_1_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred_wtd)^2 )) 
  lam_1_results[r, "pred_err_top"] <- sqrt(mean( (Y_outsample - Ypred_top)^2 )) 
  
  rm(lam_1)
  gc()

  print(paste("    Found lambda = 1 at", Sys.time()))

  ########

  z <- Sys.time()
  lam_10 <- ensm_givenhyperpar(Y, X, W, L = L, max_iter = miter, 
                         hyperpar = hyperpar, gamma_init_A = gammaA_map, gamma_init_B = gammaB_map, 
                         eta_py = 1.0, sigma_py = 0.0, 
                         lambda = 10, rho = 0.9, A_or_B_first = 0.5)
  time1 <- Sys.time() - z 
  lam_10_results[r, "TIME"] <- time1
  
  lam_output <- lam_10$output

 
  
  alpha_wtd <- rep(0, times = N)
  for(ll in 1:L)  alpha_wtd <- alpha_wtd + lam_output$alpha_particle[, ll ] * lam_output$w[ ll ] 
  beta_wtd <- rep(0, times = N)
  for(ll in 1:L)  beta_wtd <- beta_wtd + lam_output$beta_particle[, ll ] * lam_output$w[ ll ] 
  
  Lstar <- length(lam_output$unik_index)
  lam_10_results[r, "RMSE"] <- sqrt(mean( c((alpha - alpha_wtd)^2, (beta - beta_wtd)^2) )) ## todo
  lam_10_results[r, "L"] <- Lstar
  tmp_K_A <- sapply(lam_output$particle_set_A, FUN = length)
  tmp_K_B <- sapply(lam_output$particle_set_B, FUN = length)
  
  
  j_max <- which.max(lam_output$w)
  gammaA_top <- lam_output$particle_set_A[[j_max]]
  gammaB_top <- lam_output$particle_set_B[[j_max]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_top)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_top)
  
  lam_10_results[r, "RMSE_top"] <- sqrt(mean( c((alpha - lam_output$alpha_particle[,j_max])^2,(beta - lam_output$beta_particle[,j_max])^2) ))
  lam_10_results[r, "K_A_top"] <- length(gammaA_top)
  lam_10_results[r, "K_B_top"] <- length(gammaB_top)
  lam_10_results[r, "LP_top"] <- lam_output$log_posterior[j_max]
  lam_10_results[r, "RAND_A_top"] <- tmprandA$rand
  lam_10_results[r, "RAND_B_top"] <- tmprandB$rand
  lam_10_results[r, "RANDadj_A_top"] <- tmprandA$rand_adj
  lam_10_results[r, "RANDadj_B_top"] <- tmprandB$rand_adj
  lam_10_results[r, "VI_A_top"] <- tmprandA$vi
  lam_10_results[r, "VI_B_top"] <- tmprandB$vi

  tmprandAs <- sapply(lam_output$particle_set_A, rand_adj_vi, N = 400, gamma1 = gamma_alpha)
  tmprandBs <- sapply(lam_output$particle_set_B, rand_adj_vi, N = 400, gamma1 = gamma_beta)
  tmprandAs <- matrix(as.numeric(tmprandAs), ncol = L)
  tmprandBs <- matrix(as.numeric(tmprandBs), ncol = L)
  tmprandAs_avg <- rowMeans(tmprandAs)
  tmprandBs_avg <- rowMeans(tmprandBs)
  tmprandAs_wavg <- as.numeric(tmprandAs %*% matrix(lam_output$w, ncol = 1))
  tmprandBs_wavg <- as.numeric(tmprandBs %*% matrix(lam_output$w, ncol = 1))

  lam_10_results[r, "avg_RAND_A"] <- tmprandAs_avg[1]
  lam_10_results[r, "avg_RAND_B"] <- tmprandBs_avg[1]
  lam_10_results[r, "avg_RANDadj_A"] <- tmprandAs_avg[2]
  lam_10_results[r, "avg_RANDadj_B"] <- tmprandBs_avg[2]
  lam_10_results[r, "avg_VI_A"] <- tmprandAs_avg[3]
  lam_10_results[r, "avg_VI_B"] <- tmprandBs_avg[3]

  lam_10_results[r, "wavg_RAND_A"] <- tmprandAs_wavg[1]
  lam_10_results[r, "wavg_RAND_B"] <- tmprandBs_wavg[1]
  lam_10_results[r, "wavg_RANDadj_A"] <- tmprandAs_wavg[2]
  lam_10_results[r, "wavg_RANDadj_B"] <- tmprandBs_wavg[2]
  lam_10_results[r, "wavg_VI_A"] <- tmprandAs_wavg[3]
  lam_10_results[r, "wavg_VI_B"] <- tmprandBs_wavg[3]

  lam_10_results[r, "p_trueA"] <- mean(tmprandAs[1,] == 1)
  lam_10_results[r, "p_trueB"] <- mean(tmprandBs[1,] == 1)
  lam_10_results[r, "p_trueAB"] <- mean((tmprandAs[1,] == 1) & (tmprandBs[1,] == 1))

  tmpKAs <- sapply(lam_output$particle_set_A, FUN = length)
  tmpKBs <- sapply(lam_output$particle_set_B, FUN = length)

  lam_10_results[r, "avg_K_A"] <- mean(tmpKAs)
  lam_10_results[r, "avg_K_B"] <- mean(tmpKBs)
  lam_10_results[r, "wavg_K_A"] <- sum(tmpKAs * lam_output$w)
  lam_10_results[r, "wavg_K_B"] <- sum(tmpKBs * lam_output$w)
   
  Ypred_wtd <- alpha_wtd + beta_wtd * t_new_scaled
  Ypred_top <- lam_output$alpha_particle[,j_max] + lam_output$beta_particle[,j_max] * t_new_scaled
  lam_10_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred_wtd)^2 )) 
  lam_10_results[r, "pred_err_top"] <- sqrt(mean( (Y_outsample - Ypred_top)^2 )) 
  
  rm(lam_10)
  gc()
  
  print(paste("    Found lambda = 10 at", Sys.time()))

  #########

  z <- Sys.time()
  lam_100 <- ensm_givenhyperpar(Y, X, W, L = L, max_iter = miter, 
                         hyperpar = hyperpar, gamma_init_A = gammaA_map, gamma_init_B = gammaB_map, 
                         eta_py = 1.0, sigma_py = 0.0, 
                         lambda = 100, rho = 0.9, A_or_B_first = 0.5)
  time1 <- Sys.time() - z 
  lam_100_results[r, "TIME"] <- time1
  
  lam_output <- lam_100$output

 
  
  alpha_wtd <- rep(0, times = N)
  for(ll in 1:L)  alpha_wtd <- alpha_wtd + lam_output$alpha_particle[, ll ] * lam_output$w[ ll ] 
  beta_wtd <- rep(0, times = N)
  for(ll in 1:L)  beta_wtd <- beta_wtd + lam_output$beta_particle[, ll ] * lam_output$w[ ll ] 
  
  Lstar <- length(lam_output$unik_index)
  lam_100_results[r, "RMSE"] <- sqrt(mean( c((alpha - alpha_wtd)^2, (beta - beta_wtd)^2) )) ## todo
  lam_100_results[r, "L"] <- Lstar
  tmp_K_A <- sapply(lam_output$particle_set_A, FUN = length)
  tmp_K_B <- sapply(lam_output$particle_set_B, FUN = length)
  
  
  j_max <- which.max(lam_output$w)
  gammaA_top <- lam_output$particle_set_A[[j_max]]
  gammaB_top <- lam_output$particle_set_B[[j_max]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_top)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_top)
  
  lam_100_results[r, "RMSE_top"] <- sqrt(mean( c((alpha - lam_output$alpha_particle[,j_max])^2,(beta - lam_output$beta_particle[,j_max])^2) ))
  lam_100_results[r, "K_A_top"] <- length(gammaA_top)
  lam_100_results[r, "K_B_top"] <- length(gammaB_top)
  lam_100_results[r, "LP_top"] <- lam_output$log_posterior[j_max]
  lam_100_results[r, "RAND_A_top"] <- tmprandA$rand
  lam_100_results[r, "RAND_B_top"] <- tmprandB$rand
  lam_100_results[r, "RANDadj_A_top"] <- tmprandA$rand_adj
  lam_100_results[r, "RANDadj_B_top"] <- tmprandB$rand_adj
  lam_100_results[r, "VI_A_top"] <- tmprandA$vi
  lam_100_results[r, "VI_B_top"] <- tmprandB$vi

  tmprandAs <- sapply(lam_output$particle_set_A, rand_adj_vi, N = 400, gamma1 = gamma_alpha)
  tmprandBs <- sapply(lam_output$particle_set_B, rand_adj_vi, N = 400, gamma1 = gamma_beta)
  tmprandAs <- matrix(as.numeric(tmprandAs), ncol = L)
  tmprandBs <- matrix(as.numeric(tmprandBs), ncol = L)
  tmprandAs_avg <- rowMeans(tmprandAs)
  tmprandBs_avg <- rowMeans(tmprandBs)
  tmprandAs_wavg <- as.numeric(tmprandAs %*% matrix(lam_output$w, ncol = 1))
  tmprandBs_wavg <- as.numeric(tmprandBs %*% matrix(lam_output$w, ncol = 1))

  lam_100_results[r, "avg_RAND_A"] <- tmprandAs_avg[1]
  lam_100_results[r, "avg_RAND_B"] <- tmprandBs_avg[1]
  lam_100_results[r, "avg_RANDadj_A"] <- tmprandAs_avg[2]
  lam_100_results[r, "avg_RANDadj_B"] <- tmprandBs_avg[2]
  lam_100_results[r, "avg_VI_A"] <- tmprandAs_avg[3]
  lam_100_results[r, "avg_VI_B"] <- tmprandBs_avg[3]

  lam_100_results[r, "wavg_RAND_A"] <- tmprandAs_wavg[1]
  lam_100_results[r, "wavg_RAND_B"] <- tmprandBs_wavg[1]
  lam_100_results[r, "wavg_RANDadj_A"] <- tmprandAs_wavg[2]
  lam_100_results[r, "wavg_RANDadj_B"] <- tmprandBs_wavg[2]
  lam_100_results[r, "wavg_VI_A"] <- tmprandAs_wavg[3]
  lam_100_results[r, "wavg_VI_B"] <- tmprandBs_wavg[3]

  lam_100_results[r, "p_trueA"] <- mean(tmprandAs[1,] == 1)
  lam_100_results[r, "p_trueB"] <- mean(tmprandBs[1,] == 1)
  lam_100_results[r, "p_trueAB"] <- mean((tmprandAs[1,] == 1) & (tmprandBs[1,] == 1))

  tmpKAs <- sapply(lam_output$particle_set_A, FUN = length)
  tmpKBs <- sapply(lam_output$particle_set_B, FUN = length)

  lam_100_results[r, "avg_K_A"] <- mean(tmpKAs)
  lam_100_results[r, "avg_K_B"] <- mean(tmpKBs)
  lam_100_results[r, "wavg_K_A"] <- sum(tmpKAs * lam_output$w)
  lam_100_results[r, "wavg_K_B"] <- sum(tmpKBs * lam_output$w)
   
  Ypred_wtd <- alpha_wtd + beta_wtd * t_new_scaled
  Ypred_top <- lam_output$alpha_particle[,j_max] + lam_output$beta_particle[,j_max] * t_new_scaled
  lam_100_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred_wtd)^2 )) 
  lam_100_results[r, "pred_err_top"] <- sqrt(mean( (Y_outsample - Ypred_top)^2 )) 
  
  rm(lam_100)
  rm(lam_output)
  gc()
  
  print(paste("    Found lambda = 100 at", Sys.time()))


  #########                      

  for(and_i in 1:length(anderson_ks)){
    and_k <- anderson_ks[and_i]

    z <- Sys.time()
    and <- helper_anderson(C, X, W, n_draws = anderson_ndraws, 
                           num.C = as.numeric(substr(and_k, 1,1)), 
                           num.D = as.numeric(substr(and_k, 2,2)))
    time1 <- Sys.time() - z 

    anderson_results[r, "TIME", and_k] <- time1
    anderson_results[r, "RMSE", and_k] <- sqrt(mean( c((alpha - and$Amean)^2,
                                           (beta - and$Bmean)^2) ))

    anderson_results[r, "post_RAND_A", and_k] <- mean(apply(and$Cstores, MARGIN = 1, FUN = function(x) arandi(x, zA, adjust = FALSE)))
    anderson_results[r, "post_RAND_B", and_k] <- mean(apply(and$Dstores, MARGIN = 1, FUN = function(x) arandi(x, zB, adjust = FALSE)))
    anderson_results[r, "post_RANDadj_A", and_k] <- mean(apply(and$Cstores, MARGIN = 1, FUN = function(x) arandi(x, zA)))
    anderson_results[r, "post_RANDadj_B", and_k] <- mean(apply(and$Dstores, MARGIN = 1, FUN = function(x) arandi(x, zB)))
    anderson_results[r, "post_VI_A", and_k] <- mean(apply(and$Cstores, MARGIN = 1, FUN = function(x) vi.dist(x, zA)))
    anderson_results[r, "post_VI_B", and_k] <- mean(apply(and$Dstores, MARGIN = 1, FUN = function(x) vi.dist(x, zB)))
    anderson_results[r, "post_K_A", and_k] <- mean(apply(and$Cstores, MARGIN = 1, FUN = function(x) length(unique(x))))
    anderson_results[r, "post_K_B", and_k] <- mean(apply(and$Dstores, MARGIN = 1, FUN = function(x) length(unique(x))))

    Cstores_mod <- apply(and$Cstores, MARGIN = 1, FUN = Modify_partition, W); Cstores_mod = t(Cstores_mod)
    Dstores_mod <- apply(and$Dstores, MARGIN = 1, FUN = Modify_partition, W); Dstores_mod = t(Dstores_mod)

    # replacing and$Cstores and and$Dstores with Cstores_mod and Dstores_mod
    anderson_results[r, "post_mod_RAND_A", and_k] <- mean(apply(Cstores_mod, MARGIN = 1, FUN = function(x) arandi(x, zA, adjust = FALSE)))
    anderson_results[r, "post_mod_RAND_B", and_k] <- mean(apply(Dstores_mod, MARGIN = 1, FUN = function(x) arandi(x, zB, adjust = FALSE)))
    anderson_results[r, "post_mod_RANDadj_A", and_k] <- mean(apply(Cstores_mod, MARGIN = 1, FUN = function(x) arandi(x, zA)))
    anderson_results[r, "post_mod_RANDadj_B", and_k] <- mean(apply(Dstores_mod, MARGIN = 1, FUN = function(x) arandi(x, zB)))
    anderson_results[r, "post_mod_VI_A", and_k] <- mean(apply(Cstores_mod, MARGIN = 1, FUN = function(x) vi.dist(x, zA)))
    anderson_results[r, "post_mod_VI_B", and_k] <- mean(apply(Dstores_mod, MARGIN = 1, FUN = function(x) vi.dist(x, zB)))
    anderson_results[r, "post_mod_K_A", and_k] <- mean(apply(Cstores_mod, MARGIN = 1, FUN = function(x) length(unique(x))))
    anderson_results[r, "post_mod_K_B", and_k] <- mean(apply(Dstores_mod, MARGIN = 1, FUN = function(x) length(unique(x))))
    
    tmprandA <- rand_adj_vi(N, gamma_alpha, and$gamma_alpha)
    tmprandB <- rand_adj_vi(N, gamma_beta, and$gamma_beta)

    anderson_results[r, "est_RAND_A", and_k] <- tmprandA$rand
    anderson_results[r, "est_RAND_B", and_k] <- tmprandB$rand
    anderson_results[r, "est_RANDadj_A", and_k] <- tmprandA$rand_adj
    anderson_results[r, "est_RANDadj_B", and_k] <- tmprandB$rand_adj
    anderson_results[r, "est_VI_A", and_k] <- tmprandA$vi
    anderson_results[r, "est_VI_B", and_k] <- tmprandB$vi
    anderson_results[r, "est_K_A", and_k] <- length(and$gamma_alpha)
    anderson_results[r, "est_K_B", and_k] <- length(and$gamma_beta)

    z_alpha_and <- partition_list2z(and$gamma_alpha)
    z_beta_and <- partition_list2z(and$gamma_beta)
    z_alpha_mod <- Modify_partition(z_alpha_and, W)
    z_beta_mod <- Modify_partition(z_beta_and, W)
    gamma_alpha_mod <- partition_z2list(z_alpha_mod)
    gamma_beta_mod <- partition_z2list(z_beta_mod)
    tmprandA_mod <- rand_adj_vi(N, gamma_alpha, gamma_alpha_mod)
    tmprandB_mod <- rand_adj_vi(N, gamma_beta, gamma_beta_mod)

    anderson_results[r, "est_mod_RAND_A", and_k] <- tmprandA_mod$rand
    anderson_results[r, "est_mod_RAND_B", and_k] <- tmprandB_mod$rand
    anderson_results[r, "est_mod_RANDadj_A", and_k] <- tmprandA_mod$rand_adj
    anderson_results[r, "est_mod_RANDadj_B", and_k] <- tmprandB_mod$rand_adj
    anderson_results[r, "est_mod_VI_A", and_k] <- tmprandA_mod$vi
    anderson_results[r, "est_mod_VI_B", and_k] <- tmprandB_mod$vi
    anderson_results[r, "est_mod_K_A", and_k] <- length(gamma_alpha_mod)
    anderson_results[r, "est_mod_K_B", and_k] <- length(gamma_beta_mod)

    Ypred <- and$Amean + and$Bmean * t_new_scaled
    anderson_results[r, "pred_err", and_k] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 

    
    tmp <- summary_givenhyperpar(Y, X, W, 
                                hyperpar = hyperpar,
                                gamma_init_A = and$gamma_alpha, gamma_init_B = and$gamma_beta,
                                eta_input = 1.0, rho_input = 0.9)
    anderson_results[r, "LP", and_k] <- tmp$output$log_post
  }
  print(paste("    Finished anderson2017 at", Sys.time()))


  #########

  fixed_0 <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gamma_alpha, gamma_init_B = gamma_beta,
                              eta_input = 1.0, rho_input = 0.9)
  fixed_output <- fixed_0$output

  gammaA_fixed <- fixed_output$particle_set_A[[1]]
  gammaB_fixed <- fixed_output$particle_set_B[[1]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_fixed)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_fixed)
  
  fixed_0_results[r,"RMSE"] <- sqrt(mean( c((alpha - fixed_output$alpha)^2,
                                            (beta - fixed_output$beta)^2) ))
  fixed_0_results[r, "RAND_A"] <- tmprandA$rand
  fixed_0_results[r, "RAND_B"] <- tmprandB$rand
  fixed_0_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  fixed_0_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  fixed_0_results[r, "VI_A"] <- tmprandA$vi
  fixed_0_results[r, "VI_B"] <- tmprandB$vi
  fixed_0_results[r, "K_A"] <- length(gammaA_fixed)
  fixed_0_results[r, "K_B"] <- length(gammaB_fixed)
  fixed_0_results[r, "LP"] <- fixed_output$log_post

  Ypred <- fixed_output$alpha + fixed_output$beta * t_new_scaled
  fixed_0_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 

  print(paste("    Finished fixed_0 at", Sys.time()))

  

  fixed_1 <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gamma_1, gamma_init_B = gamma_1,
                              eta_input = 1.0, rho_input = 0.9)
  fixed_output <- fixed_1$output

  gammaA_fixed <- fixed_output$particle_set_A[[1]]
  gammaB_fixed <- fixed_output$particle_set_B[[1]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_fixed)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_fixed)
  
  fixed_1_results[r,"RMSE"] <- sqrt(mean( c((alpha - fixed_output$alpha)^2,
                                            (beta - fixed_output$beta)^2) ))
  fixed_1_results[r, "RAND_A"] <- tmprandA$rand
  fixed_1_results[r, "RAND_B"] <- tmprandB$rand
  fixed_1_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  fixed_1_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  fixed_1_results[r, "VI_A"] <- tmprandA$vi
  fixed_1_results[r, "VI_B"] <- tmprandB$vi
  fixed_1_results[r, "K_A"] <- length(gammaA_fixed)
  fixed_1_results[r, "K_B"] <- length(gammaB_fixed)
  fixed_1_results[r, "LP"] <- fixed_output$log_post

  Ypred <- fixed_output$alpha + fixed_output$beta * t_new_scaled
  fixed_1_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 
  print(paste("    Finished fixed_1 at", Sys.time()))

  
  
  fixed_n <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gamma_n, gamma_init_B = gamma_n,
                              eta_input = 1.0, rho_input = 0.9)
  fixed_output <- fixed_n$output

  gammaA_fixed <- fixed_output$particle_set_A[[1]]
  gammaB_fixed <- fixed_output$particle_set_B[[1]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_fixed)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_fixed)
  
  fixed_n_results[r,"RMSE"] <- sqrt(mean( c((alpha - fixed_output$alpha)^2,
                                            (beta - fixed_output$beta)^2) ))
  fixed_n_results[r, "RAND_A"] <- tmprandA$rand
  fixed_n_results[r, "RAND_B"] <- tmprandB$rand
  fixed_n_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  fixed_n_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  fixed_n_results[r, "VI_A"] <- tmprandA$vi
  fixed_n_results[r, "VI_B"] <- tmprandB$vi
  fixed_n_results[r, "K_A"] <- length(gammaA_fixed)
  fixed_n_results[r, "K_B"] <- length(gammaB_fixed)
  fixed_n_results[r, "LP"] <- fixed_output$log_post

  Ypred <- fixed_output$alpha + fixed_output$beta * t_new_scaled
  fixed_n_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 
  print(paste("    Finished fixed_n at", Sys.time()))
  

  ### should we compute the silhouette on the "modified" partition or on the original one???? 
  z <- Sys.time()
  km_silhouetteA <- numeric(KMsplits)
  km_silhouetteB <- numeric(KMsplits)
  dA <- dist(alphas_mle)
  dB <- dist(alphas_mle)
  for(k in 2:KMsplits){
    km <- kmeans(alphas_mle, centers = k, nstart = 25)
    tmp_sil <- silhouette(x = km$cluster, dist = dA)
    km_silhouetteA[k] <- mean(tmp_sil[, 3])

    km <- kmeans(betas_mle, centers = k, nstart = 25)
    tmp_sil <- silhouette(x = km$cluster, dist = dB)
    km_silhouetteB[k] <- mean(tmp_sil[, 3])
  }
  kA <- which.max(km_silhouetteA)
  kB <- which.max(km_silhouetteB)

  kmA <- kmeans(alphas_mle, centers = kA, nstart = 25)
  kmB <- kmeans(betas_mle, centers = kB, nstart = 25)
  time1 <- Sys.time() - z 
  km_results[r, "TIME"] <- time1

  zA <- Modify_partition(kmA$cluster, W)
  zB <- Modify_partition(kmB$cluster, W)

  gammaA_top <- list(); count <- 1
  for(k in unique(zA)){
    gammaA_top[[count]] <- which(zA == k)
    count <- count + 1
  }
  gammaB_top <- list(); count <- 1
  for(k in unique(zB)){
    gammaB_top[[count]] <- which(zB == k)
    count <- count + 1
  }

  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_top)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_top)
  
  km_results[r, "RMSE_mle"] <- sqrt(mean( c((alpha - alphas_mle)^2,(beta - betas_mle)^2) )) # kmeans does not provide any shrinkage
  km_results[r, "K_A"] <- length(gammaA_top)
  km_results[r, "K_B"] <- length(gammaB_top)
  km_results[r, "RAND_A"] <- tmprandA$rand
  km_results[r, "RAND_B"] <- tmprandB$rand
  km_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  km_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  km_results[r, "VI_A"] <- tmprandA$vi
  km_results[r, "VI_B"] <- tmprandB$vi
  
  tmp_km <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gammaA_top, gamma_init_B = gammaB_top,
                              eta_input = 1.0, rho_input = 0.9)
  km_results[r, "RMSE"] <- sqrt(mean( c((alpha - tmp_km$output$alpha)^2,(beta - tmp_km$output$beta)^2) )) # kmeans does not provide any shrinkage
  km_results[r, "LP"] <- tmp_km$output$log_post

  Ypred <- tmp_km$output$alpha + tmp_km$output$beta * t_new_scaled
  km_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 
  
  gc()

  #########

  ### should we compute the silhouette on the "modified" partition or on the original one???? 
  z <- Sys.time()
  sc_silhouetteA <- numeric(SCsplits)
  sc_silhouetteB <- numeric(SCsplits)
  for(k in 2:SCsplits){
    sc <- spectral_clustering(alphas_mle, W, k)
    tmp_sil <- silhouette(x = sc$cluster, dist = dA)
    sc_silhouetteA[k] <- mean(tmp_sil[, 3])

    sc <- spectral_clustering(betas_mle, W, k)
    tmp_sil <- silhouette(x = sc$cluster, dist = dB)
    sc_silhouetteB[k] <- mean(tmp_sil[, 3])
  }
  kA <- which.max(sc_silhouetteA)
  kB <- which.max(sc_silhouetteB)

  scA <- spectral_clustering(alphas_mle, W, kA)
  scB <- spectral_clustering(betas_mle, W, kB)
  time1 <- Sys.time() - z 
  sc_results[r, "TIME"] <- time1

  zA <- Modify_partition(scA$cluster, W)
  zB <- Modify_partition(scB$cluster, W)

  gammaA_top <- list(); count <- 1
  for(k in unique(zA)){
    gammaA_top[[count]] <- which(zA == k)
    count <- count + 1
  }
  gammaB_top <- list(); count <- 1
  for(k in unique(zB)){
    gammaB_top[[count]] <- which(zB == k)
    count <- count + 1
  }

  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_top)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_top)
  
  sc_results[r, "RMSE_mle"] <- sqrt(mean( c((alpha - alphas_mle)^2,(beta - betas_mle)^2) )) # kmeans does not provide any shrinkage
  sc_results[r, "K_A"] <- length(gammaA_top)
  sc_results[r, "K_B"] <- length(gammaB_top)
  sc_results[r, "RAND_A"] <- tmprandA$rand
  sc_results[r, "RAND_B"] <- tmprandB$rand
  sc_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  sc_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  sc_results[r, "VI_A"] <- tmprandA$vi
  sc_results[r, "VI_B"] <- tmprandB$vi
  
  tmp_sc <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gammaA_top, gamma_init_B = gammaB_top,
                              eta_input = 1.0, rho_input = 0.9)
  sc_results[r, "RMSE"] <- sqrt(mean( c((alpha - tmp_sc$output$alpha)^2,(beta - tmp_sc$output$beta)^2) )) # kmeans does not provide any shrinkage
  sc_results[r, "LP"] <- tmp_sc$output$log_post

  Ypred <- tmp_sc$output$alpha + tmp_sc$output$beta * t_new_scaled
  sc_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 

  gc()

  print(paste("    Finished sc at", Sys.time()))

}

assign(paste0("map_results_sim", sim_number, "_", batch), map_results)
assign(paste0("lam_1_results_sim", sim_number, "_", batch), lam_1_results)
assign(paste0("lam_10_results_sim", sim_number, "_", batch), lam_10_results)
assign(paste0("lam_100_results_sim", sim_number, "_", batch), lam_100_results)
assign(paste0("anderson_results_sim", sim_number, "_", batch), anderson_results)
assign(paste0("fixed_0_results_sim", sim_number, "_", batch), fixed_0_results)
assign(paste0("fixed_1_results_sim", sim_number, "_", batch), fixed_1_results)
assign(paste0("fixed_n_results_sim", sim_number, "_", batch), fixed_n_results)
assign(paste0("km_results_sim", sim_number, "_",batch), km_results)
assign(paste0("sc_results_sim", sim_number, "_", batch), sc_results)


save_list <- paste0(c("lam_1_", "lam_10_", "lam_100_",
                      "fixed_0_", "fixed_1_", "fixed_n_",
                      "km_", "sc_", "map_",
                      "anderson_"), "results_sim", sim_number, "_", batch)
save(list = save_list, file = paste0("results/new2newX_sim", sim_number, "_", batch, ".RData"))
