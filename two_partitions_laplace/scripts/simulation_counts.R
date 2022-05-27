library(Rcpp)
library(RcppArmadillo)
library(mcclust)

source("src/partition_functions.R")
source("src/ensm_hyperpar.R")
source("src/ensm_givenhyperpar.R")
# source("src/summary_givenhyperpar.R")


## need to work on the data generation process
load("data/alphabetas.RData")
load("data/partitions.RData")

N <- 100
T <- 12
batch_size <- 4

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

sim_number <- 3 # 1, 2, 3
batch <- 1

ALPHA <- get(paste0("alpha_", sim_number))
BETA <- get(paste0("beta_", sim_number))

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

fixed_metrics <- c("RMSE", "RAND_A","RAND_B","RANDadj_A", "RANDadj_B", "VI_A", "VI_B", "K_A","K_B","pred_err", "LP")
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

miter <- 100
L <- 10
r = 2

print(paste("  Starting r = ", r, "at", Sys.time()))
set.seed(seed_seq[sim_number] * 100 + 10 * (batch - 1) + r)
alpha <- ALPHA[, r + (batch-1)*batch_size]
beta <- BETA[, r + (batch-1)*batch_size]
mu <- matrix(nrow = N, ncol = T)
C <- matrix(nrow = N, ncol = T)
for(i in 1:N) mu[i,] <- alpha[i] + beta[i] * X[i,]
lambda = exp(mu)
for(i in 1:N) C[i,] <- rpois(n = T, lambda = lambda[i,])
# Y <- ihs_transform(C)
 
# todo
C_outsample <- rpois(N, lambda = exp(alpha + beta * t_new_scaled))
# Y_outsample <- ihs_transform(C_outsample)

if(any(C < 0)){
  warning("Some counts were negative: specifically in locations: ",which(C<0), " with values ", C[which(C<0)])
  C[which(C<0)] <- 0
}
if(any(C_outsample < 0)){
  warning("Some (outsample) counts were negative: specifically in locations: ",which(C_outsample<0), " with values ", C_outsample[which(C_outsample<0)])
  C_outsample[which(C_outsample<0)] <- 0
}

# betas_mle <- numeric(N)
# for(i in 1:N)
#   betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
# alphas_mle <- rowMeans(Y)

## let's run MAP first, and decide the hyper-parameters.
z <- Sys.time()
map <- ensm_hyperpar(C, X, W, L = 1, max_iter = miter, 
                       eta_py = 1.0, sigma_py = 0.0, 
                       lambda = 1, rho = 0.9)
time1 <- Sys.time() - z 
map_results[r, "TIME"] <- time1
beepr::beep()

## let's set the hyper-parameters for the next ones:
hyperpar <- map$hyperpar # a1, a2, b1, b2, alpha_sigma, beta_sigma + estimated sigma2_laplace

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

Cpred_top <- exp(map$alpha_particle[,1] + map$beta_particle[,1] * t_new_scaled)
map_results[r, "pred_err"] <- sqrt(mean( (C_outsample - Cpred_top)^2 ))

# rm(map)
# gc()

print(paste("    Finished MAP at", Sys.time()))
#########

z <- Sys.time()
lam_1 <- ensm_givenhyperpar(C, X, W, L = L, max_iter = miter, 
                       hyperpar = hyperpar, gamma_init_A = gammaA_map, gamma_init_B = gammaB_map, 
                       eta_py = 1.0, sigma_py = 0.0, 
                       lambda = 1, rho = 0.9, A_or_B_first = 0.5)
time1 <- Sys.time() - z 
lam_1_results[r, "TIME"] <- time1
beepr::beep()

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

Cpred_wtd <- exp(alpha_wtd + beta_wtd * t_new_scaled)
Cpred_top <- exp(lam_output$alpha_particle[,j_max] + lam_output$beta_particle[,j_max] * t_new_scaled)
lam_1_results[r, "pred_err"] <- sqrt(mean( (C_outsample - Cpred_wtd)^2 ))
lam_1_results[r, "pred_err_top"] <- sqrt(mean( (C_outsample - Cpred_top)^2 ))

# rm(lam_1)
# gc()

print(paste("    Found lambda = 1 at", Sys.time()))

########

z <- Sys.time()
lam_10 <- ensm_givenhyperpar(C, X, W, L = L, max_iter = miter, 
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
 
Cpred_wtd <- exp(alpha_wtd + beta_wtd * t_new_scaled)
Cpred_top <- exp(lam_output$alpha_particle[,j_max] + lam_output$beta_particle[,j_max] * t_new_scaled)
lam_10_results[r, "pred_err"] <- sqrt(mean( (C_outsample - Cpred_wtd)^2 )) 
lam_10_results[r, "pred_err_top"] <- sqrt(mean( (C_outsample - Cpred_top)^2 )) 

# rm(lam_10)
# gc()

print(paste("    Found lambda = 10 at", Sys.time()))

#########

z <- Sys.time()
lam_100 <- ensm_givenhyperpar(C, X, W, L = L, max_iter = miter, 
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
 
Cpred_wtd <- exp(alpha_wtd + beta_wtd * t_new_scaled)
Cpred_top <- exp(lam_output$alpha_particle[,j_max] + lam_output$beta_particle[,j_max] * t_new_scaled)
lam_100_results[r, "pred_err"] <- sqrt(mean( (C_outsample - Cpred_wtd)^2 )) 
lam_100_results[r, "pred_err_top"] <- sqrt(mean( (C_outsample - Cpred_top)^2 )) 

# rm(lam_100)
# rm(lam_output)
# gc()

print(paste("    Found lambda = 100 at", Sys.time()))



assign(paste0("map_results_sim", sim_number, "_", batch), map_results)
assign(paste0("lam_1_results_sim", sim_number, "_", batch), lam_1_results)
assign(paste0("lam_10_results_sim", sim_number, "_", batch), lam_10_results)
assign(paste0("lam_100_results_sim", sim_number, "_", batch), lam_100_results)


save_list <- paste0(c("lam_1_", "lam_10_", "lam_100_",
                      "map_"), "results_sim", sim_number, "_", batch)
save_list <- c(save_list, "map", "lam_1", "lam_10", "lam_100")
save(list = save_list, file = paste0("results/counts_newX_sim", sim_number, "_", batch, ".RData"))
