library(Rcpp)
library(RcppArmadillo)
library(mcclust)

source("src/partition_functions.R")
source("src/equalpart/ensm_hyperpar.R")
source("src/equalpart/ensm_givenhyperpar.R")
source("src/equalpart/summary_givenhyperpar.R")

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

batch_size <- 4

particle_metrics <- c("RAND_A_top", "RAND_B_top","RAND_new_top","RANDadj_A_top", "RANDadj_B_top","RANDadj_new_top", "VI_A_top", "VI_B_top", "VI_new_top", "LP_top", 
                      "avg_RAND_A", "avg_RAND_B","avg_RAND_new","avg_RANDadj_A", "avg_RANDadj_B","avg_RANDadj_new","avg_VI_A", "avg_VI_B","avg_VI_new",
                      "wavg_RAND_A", "wavg_RAND_B","wavg_RAND_new","wavg_RANDadj_A", "wavg_RANDadj_B","wavg_RANDadj_new","wavg_VI_A", "wavg_VI_B","wavg_VI_new",
                      "RMSE","RMSE_top", "pred_err", "pred_err_top",
                      "K_A_top","K_B_top","avg_K_A","avg_K_B","wavg_K_A","wavg_K_B",
                      "p_trueA","p_trueB","p_trueAB","p_true_new",
                      "L", "TIME")
lam_1_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
lam_10_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
lam_100_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))

fixed_metrics <- c("RMSE", "RAND_A","RAND_B","RAND_new","RANDadj_A", "RANDadj_B", "RANDadj_new", "VI_A", "VI_B","VI_new", "K_A","K_B","pred_err", "LP")
fixed_0_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_1_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_n_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
fixed_new_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
map_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics) + 1, dimnames = list(c(), c(fixed_metrics, "TIME")))

ihs_transform <- function(x){
  return( log(x + sqrt(x^2 + 1)) - log(2))
}

t <- 12
t_new <- 13
tmp_time <- 1:t - mean(1:t)
t_new_scaled <- (t_new - mean(1:t))/sd(tmp_time)

zA <- numeric(N)
for(i in 1:length(gamma_alpha)){
  zA[ gamma_alpha[[i]] ] <- i
}
zB <- numeric(N)
for(i in 1:length(gamma_beta)){
  zB[ gamma_beta[[i]] ] <- i
}

gamma_new <- list()
tmp <- 1
for(i in 1:length(gamma_alpha)){
  for(j in 1:length(gamma_beta)){
    subset_part <- intersect(gamma_alpha[[i]], gamma_beta[[j]])
    if(length(subset_part)>0){
      gamma_new[[tmp]] <- subset_part
      tmp <- tmp + 1
    }
  }
}

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

  map_orig <- map
  map <- map$adjusted

  gamma_map <- map$particle_set_A[[1]]
  # gammaB_map <- map$particle_set_B[[1]]

  map_results[r, "RMSE"] <- sqrt(mean( c((alpha - map$alpha_particle[,1])^2,
                                         (beta - map$beta_particle[,1])^2) ))
  tmprandA <- rand_adj_vi(N, gamma_alpha, gamma_map)
  tmprandB <- rand_adj_vi(N, gamma_beta, gamma_map)
  tmprand_new <- rand_adj_vi(N, gamma_new, gamma_map)

  map_results[r, "RAND_A"] <- tmprandA$rand
  map_results[r, "RAND_B"] <- tmprandB$rand
  map_results[r, "RAND_new"] <- tmprand_new$rand
  map_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  map_results[r, "RANDadj_B"] <- tmprandB$rand_adj 
  map_results[r, "RANDadj_new"] <- tmprand_new$rand_adj 
  map_results[r, "VI_A"] <- tmprandA$vi
  map_results[r, "VI_B"] <- tmprandB$vi
  map_results[r, "VI_new"] <- tmprand_new$vi
  map_results[r, "K_A"] <- length(gamma_map)
  map_results[r, "K_B"] <- length(gamma_map)
  map_results[r, "LP"] <- map$log_posterior[1]

  Ypred_top <- map$alpha_particle[,1] + map$beta_particle[,1] * t_new_scaled
  map_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred_top)^2 )) 

  rm(map)
  gc()

  print(paste("    Finished MAP at", Sys.time()))
  #########

  z <- Sys.time()
  lam_1 <- ensm_givenhyperpar(Y, X, W, L = L, max_iter = miter, 
                         hyperpar = hyperpar, gamma_init = gamma_map, 
                         eta_py = 1.0, sigma_py = 0.0, 
                         lambda = 1, rho = 0.9)
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
  tmprand_new <- rand_adj_vi(N, gamma_new, gammaB_top)
  
  lam_1_results[r, "RMSE_top"] <- sqrt(mean( c((alpha - lam_output$alpha_particle[,j_max])^2,(beta - lam_output$beta_particle[,j_max])^2) ))
  lam_1_results[r, "K_A_top"] <- length(gammaA_top)
  lam_1_results[r, "K_B_top"] <- length(gammaB_top)
  lam_1_results[r, "LP_top"] <- lam_output$log_posterior[j_max]
  lam_1_results[r, "RAND_A_top"] <- tmprandA$rand
  lam_1_results[r, "RAND_B_top"] <- tmprandB$rand
  lam_1_results[r, "RAND_new_top"] <- tmprand_new$rand
  lam_1_results[r, "RANDadj_A_top"] <- tmprandA$rand_adj
  lam_1_results[r, "RANDadj_B_top"] <- tmprandB$rand_adj
  lam_1_results[r, "RANDadj_new_top"] <- tmprand_new$rand_adj
  lam_1_results[r, "VI_A_top"] <- tmprandA$vi
  lam_1_results[r, "VI_B_top"] <- tmprandB$vi
  lam_1_results[r, "VI_new_top"] <- tmprand_new$vi

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

  tmprand_news <- sapply(lam_output$particle_set_B, rand_adj_vi, N = 400, gamma1 = gamma_new)
  tmprand_news <- matrix(as.numeric(tmprand_news), ncol = L)
  tmprand_news_avg <- rowMeans(tmprand_news)
  tmprand_news_wavg <- as.numeric(tmprand_news %*% matrix(lam_output$w, ncol = 1))
  lam_1_results[r, "avg_RAND_new"] <- tmprand_news_avg[1]
  lam_1_results[r, "avg_RANDadj_new"] <- tmprand_news_avg[2]
  lam_1_results[r, "avg_VI_new"] <- tmprand_news_avg[3]
  lam_1_results[r, "wavg_RAND_new"] <- tmprand_news_wavg[1]
  lam_1_results[r, "wavg_RANDadj_new"] <- tmprand_news_wavg[2]
  lam_1_results[r, "wavg_VI_new"] <- tmprand_news_wavg[3]

  lam_1_results[r, "p_trueA"] <- mean(tmprandAs[1,] == 1)
  lam_1_results[r, "p_trueB"] <- mean(tmprandBs[1,] == 1)
  lam_1_results[r, "p_trueAB"] <- mean((tmprandAs[1,] == 1) & (tmprandBs[1,] == 1))
  lam_1_results[r, "p_true_new"] <- mean(tmprand_news[1,] == 1)

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
                         hyperpar = hyperpar, gamma_init = gamma_map, 
                         eta_py = 1.0, sigma_py = 0.0, 
                         lambda = 10, rho = 0.9)
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
  tmprand_new <- rand_adj_vi(N, gamma_new, gammaB_top)
  
  lam_10_results[r, "RMSE_top"] <- sqrt(mean( c((alpha - lam_output$alpha_particle[,j_max])^2,(beta - lam_output$beta_particle[,j_max])^2) ))
  lam_10_results[r, "K_A_top"] <- length(gammaA_top)
  lam_10_results[r, "K_B_top"] <- length(gammaB_top)
  lam_10_results[r, "LP_top"] <- lam_output$log_posterior[j_max]
  lam_10_results[r, "RAND_A_top"] <- tmprandA$rand
  lam_10_results[r, "RAND_B_top"] <- tmprandB$rand
  lam_10_results[r, "RAND_new_top"] <- tmprand_new$rand
  lam_10_results[r, "RANDadj_A_top"] <- tmprandA$rand_adj
  lam_10_results[r, "RANDadj_B_top"] <- tmprandB$rand_adj
  lam_10_results[r, "RANDadj_new_top"] <- tmprand_new$rand_adj
  lam_10_results[r, "VI_A_top"] <- tmprandA$vi
  lam_10_results[r, "VI_B_top"] <- tmprandB$vi
  lam_10_results[r, "VI_new_top"] <- tmprand_new$vi

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

  tmprand_news <- sapply(lam_output$particle_set_B, rand_adj_vi, N = 400, gamma1 = gamma_new)
  tmprand_news <- matrix(as.numeric(tmprand_news), ncol = L)
  tmprand_news_avg <- rowMeans(tmprand_news)
  tmprand_news_wavg <- as.numeric(tmprand_news %*% matrix(lam_output$w, ncol = 1))
  lam_10_results[r, "avg_RAND_new"] <- tmprand_news_avg[1]
  lam_10_results[r, "avg_RANDadj_new"] <- tmprand_news_avg[2]
  lam_10_results[r, "avg_VI_new"] <- tmprand_news_avg[3]
  lam_10_results[r, "wavg_RAND_new"] <- tmprand_news_wavg[1]
  lam_10_results[r, "wavg_RANDadj_new"] <- tmprand_news_wavg[2]
  lam_10_results[r, "wavg_VI_new"] <- tmprand_news_wavg[3]

  lam_10_results[r, "p_trueA"] <- mean(tmprandAs[1,] == 1)
  lam_10_results[r, "p_trueB"] <- mean(tmprandBs[1,] == 1)
  lam_10_results[r, "p_trueAB"] <- mean((tmprandAs[1,] == 1) & (tmprandBs[1,] == 1))
  lam_10_results[r, "p_true_new"] <- mean(tmprand_news[1,] == 1)

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
                         hyperpar = hyperpar, gamma_init = gamma_map, 
                         eta_py = 1.0, sigma_py = 0.0, 
                         lambda = 100, rho = 0.9)
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
  tmprand_new <- rand_adj_vi(N, gamma_new, gammaB_top)
  
  lam_100_results[r, "RMSE_top"] <- sqrt(mean( c((alpha - lam_output$alpha_particle[,j_max])^2,(beta - lam_output$beta_particle[,j_max])^2) ))
  lam_100_results[r, "K_A_top"] <- length(gammaA_top)
  lam_100_results[r, "K_B_top"] <- length(gammaB_top)
  lam_100_results[r, "LP_top"] <- lam_output$log_posterior[j_max]
  lam_100_results[r, "RAND_A_top"] <- tmprandA$rand
  lam_100_results[r, "RAND_B_top"] <- tmprandB$rand
  lam_100_results[r, "RAND_new_top"] <- tmprand_new$rand
  lam_100_results[r, "RANDadj_A_top"] <- tmprandA$rand_adj
  lam_100_results[r, "RANDadj_B_top"] <- tmprandB$rand_adj
  lam_100_results[r, "RANDadj_new_top"] <- tmprand_new$rand_adj
  lam_100_results[r, "VI_A_top"] <- tmprandA$vi
  lam_100_results[r, "VI_B_top"] <- tmprandB$vi
  lam_100_results[r, "VI_new_top"] <- tmprand_new$vi

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

  tmprand_news <- sapply(lam_output$particle_set_B, rand_adj_vi, N = 400, gamma1 = gamma_new)
  tmprand_news <- matrix(as.numeric(tmprand_news), ncol = L)
  tmprand_news_avg <- rowMeans(tmprand_news)
  tmprand_news_wavg <- as.numeric(tmprand_news %*% matrix(lam_output$w, ncol = 1))
  lam_100_results[r, "avg_RAND_new"] <- tmprand_news_avg[1]
  lam_100_results[r, "avg_RANDadj_new"] <- tmprand_news_avg[2]
  lam_100_results[r, "avg_VI_new"] <- tmprand_news_avg[3]
  lam_100_results[r, "wavg_RAND_new"] <- tmprand_news_wavg[1]
  lam_100_results[r, "wavg_RANDadj_new"] <- tmprand_news_wavg[2]
  lam_100_results[r, "wavg_VI_new"] <- tmprand_news_wavg[3]

  lam_100_results[r, "p_trueA"] <- mean(tmprandAs[1,] == 1)
  lam_100_results[r, "p_trueB"] <- mean(tmprandBs[1,] == 1)
  lam_100_results[r, "p_trueAB"] <- mean((tmprandAs[1,] == 1) & (tmprandBs[1,] == 1))
  lam_100_results[r, "p_true_new"] <- mean(tmprand_news[1,] == 1)

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

  #######

  fixed_0 <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gamma_alpha, gamma_init_B = gamma_beta,
                              eta_input = 1.0, rho_input = 0.9)
  fixed_output <- fixed_0$output

  gammaA_fixed <- fixed_output$particle_set_A[[1]]
  gammaB_fixed <- fixed_output$particle_set_B[[1]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_fixed)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_fixed)
  tmprand_new <- rand_adj_vi(N, gamma_new, gammaA_fixed)
  
  fixed_0_results[r,"RMSE"] <- sqrt(mean( c((alpha - fixed_output$alpha)^2,
                                            (beta - fixed_output$beta)^2) ))
  fixed_0_results[r, "RAND_A"] <- tmprandA$rand
  fixed_0_results[r, "RAND_B"] <- tmprandB$rand
  fixed_0_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  fixed_0_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  fixed_0_results[r, "VI_A"] <- tmprandA$vi
  fixed_0_results[r, "VI_B"] <- tmprandB$vi
  fixed_0_results[r, "RAND_new"] <- tmprand_new$rand
  fixed_0_results[r, "RANDadj_new"] <- tmprand_new$rand_adj
  fixed_0_results[r, "VI_new"] <- tmprand_new$vi
  fixed_0_results[r, "K_A"] <- length(gammaA_fixed)
  fixed_0_results[r, "K_B"] <- length(gammaB_fixed)
  fixed_0_results[r, "LP"] <- fixed_output$log_post

  Ypred <- fixed_output$alpha + fixed_output$beta * t_new_scaled
  fixed_0_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 

  print(paste("    Finished fixed_0 at", Sys.time()))


  fixed_new <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gamma_new, gamma_init_B = gamma_new,
                              eta_input = 1.0, rho_input = 0.9)
  fixed_output <- fixed_new$output

  gammaA_fixed <- fixed_output$particle_set_A[[1]]
  gammaB_fixed <- fixed_output$particle_set_B[[1]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_fixed)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_fixed)
  tmprand_new <- rand_adj_vi(N, gamma_new, gammaA_fixed)
  
  fixed_new_results[r,"RMSE"] <- sqrt(mean( c((alpha - fixed_output$alpha)^2,
                                            (beta - fixed_output$beta)^2) ))
  fixed_new_results[r, "RAND_A"] <- tmprandA$rand
  fixed_new_results[r, "RAND_B"] <- tmprandB$rand
  fixed_new_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  fixed_new_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  fixed_new_results[r, "VI_A"] <- tmprandA$vi
  fixed_new_results[r, "VI_B"] <- tmprandB$vi
  fixed_new_results[r, "RAND_new"] <- tmprand_new$rand
  fixed_new_results[r, "RANDadj_new"] <- tmprand_new$rand_adj
  fixed_new_results[r, "VI_new"] <- tmprand_new$vi
  fixed_new_results[r, "K_A"] <- length(gammaA_fixed)
  fixed_new_results[r, "K_B"] <- length(gammaB_fixed)
  fixed_new_results[r, "LP"] <- fixed_output$log_post

  Ypred <- fixed_output$alpha + fixed_output$beta * t_new_scaled
  fixed_new_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 

  print(paste("    Finished fixed_new at", Sys.time()))

  

  fixed_1 <- summary_givenhyperpar(Y, X, W, 
                              hyperpar = hyperpar,
                              gamma_init_A = gamma_1, gamma_init_B = gamma_1,
                              eta_input = 1.0, rho_input = 0.9)
  fixed_output <- fixed_1$output

  gammaA_fixed <- fixed_output$particle_set_A[[1]]
  gammaB_fixed <- fixed_output$particle_set_B[[1]]
  tmprandA <- rand_adj_vi(N, gamma_alpha, gammaA_fixed)
  tmprandB <- rand_adj_vi(N, gamma_beta, gammaB_fixed)
  tmprand_new <- rand_adj_vi(N, gamma_new, gammaA_fixed)
  
  fixed_1_results[r,"RMSE"] <- sqrt(mean( c((alpha - fixed_output$alpha)^2,
                                            (beta - fixed_output$beta)^2) ))
  fixed_1_results[r, "RAND_A"] <- tmprandA$rand
  fixed_1_results[r, "RAND_B"] <- tmprandB$rand
  fixed_1_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  fixed_1_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  fixed_1_results[r, "VI_A"] <- tmprandA$vi
  fixed_1_results[r, "VI_B"] <- tmprandB$vi
  fixed_1_results[r, "RAND_new"] <- tmprand_new$rand
  fixed_1_results[r, "RANDadj_new"] <- tmprand_new$rand_adj
  fixed_1_results[r, "VI_new"] <- tmprand_new$vi
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
  tmprand_new <- rand_adj_vi(N, gamma_new, gammaA_fixed)
  
  fixed_n_results[r,"RMSE"] <- sqrt(mean( c((alpha - fixed_output$alpha)^2,
                                            (beta - fixed_output$beta)^2) ))
  fixed_n_results[r, "RAND_A"] <- tmprandA$rand
  fixed_n_results[r, "RAND_B"] <- tmprandB$rand
  fixed_n_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  fixed_n_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  fixed_n_results[r, "VI_A"] <- tmprandA$vi
  fixed_n_results[r, "VI_B"] <- tmprandB$vi
  fixed_n_results[r, "RAND_new"] <- tmprand_new$rand
  fixed_n_results[r, "RANDadj_new"] <- tmprand_new$rand_adj
  fixed_n_results[r, "VI_new"] <- tmprand_new$vi
  fixed_n_results[r, "K_A"] <- length(gammaA_fixed)
  fixed_n_results[r, "K_B"] <- length(gammaB_fixed)
  fixed_n_results[r, "LP"] <- fixed_output$log_post

  Ypred <- fixed_output$alpha + fixed_output$beta * t_new_scaled
  fixed_n_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 
  print(paste("    Finished fixed_n at", Sys.time()))

}

assign(paste0("map_results_sim", sim_number, "_", batch), map_results)
assign(paste0("lam_1_results_sim", sim_number, "_", batch), lam_1_results)
assign(paste0("lam_10_results_sim", sim_number, "_", batch), lam_10_results)
assign(paste0("lam_100_results_sim", sim_number, "_", batch), lam_100_results)
assign(paste0("fixed_0_results_sim", sim_number, "_", batch), fixed_0_results)
assign(paste0("fixed_1_results_sim", sim_number, "_", batch), fixed_1_results)
assign(paste0("fixed_n_results_sim", sim_number, "_", batch), fixed_n_results)
assign(paste0("fixed_new_results_sim", sim_number, "_", batch), fixed_new_results)

save_list <- paste0(c("lam_1_", "lam_10_", "lam_100_",
                      "fixed_0_", "fixed_1_", "fixed_n_", "fixed_new_",
                      "map_"), "results_sim", sim_number, "_", batch)
save(list = save_list, file = paste0("results/equalpart_sim", sim_number, "_", batch, ".RData"))
