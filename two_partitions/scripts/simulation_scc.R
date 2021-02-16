library(glmnet)
library(igraph)
source("src/partition_functions.R")

load("data/alphabetas.RData")
load("data/partitions.RData")
load("data/ssc_design.RData")

N <- 400
T <- 12

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

args <- commandArgs(TRUE)
sim_number <- as.numeric(args[1])
batch <- as.numeric(args[2])
batch_size <- 4

ALPHA <- get(paste0("alpha_", sim_number))
BETA <- get(paste0("beta_", sim_number))


ssc_metrics <- c("RMSE","K_A", "K_B", "RAND_A", "RAND_B", "RANDadj_A", "RANDadj_B", "VI_A", "VI_B", "pred_err")
ssc_results <- matrix(nrow = batch_size, ncol = length(ssc_metrics), dimnames = list(c(), ssc_metrics))

tmp_time <- 1:T
t_new <- T+1
t_new_scaled <- (t_new - mean(tmp_time))/sd(tmp_time)

for(r in 1:batch_size){
  print(paste("  Starting r = ", r, "at", Sys.time()))
  set.seed(seed_seq[sim_number] * 100 + 10 * (batch - 1) + r)
  alpha <- ALPHA[, r + (batch-1)*batch_size]
  beta <- BETA[, r + (batch-1)*batch_size]
  
  Y_vec <- rep(NA, times = N*T)
  for(i in 1:N){
    Y_vec[(1 + (i-1)*T):(i*T)] <- alpha[i] + beta[i] * tmp_time + sqrt(sigma2) * rnorm(T, 0, 1)
  }
  Y_outsample <- rnorm(N, mean = alpha + beta * t_new_scaled, sd = sqrt(sigma2))
  
  cv_fit <- cv.glmnet(x = X_tilde, y = Y_vec, intercept = FALSE, 
                      penalty.factor = penalty_factor)
  # even though we specify intercept = FALSE
  # glmnet still returns an intercept
  theta_hat <- coef(cv_fit, s = "lambda.1se")[-1]
  
  alpha_hat <- H_tilde_inv %*% theta_hat[1:N]
  beta_hat <- H_tilde_inv %*% theta_hat[(N+1):(2*N)]
  
  Y_hat <- alpha_hat + t_new_scaled * beta_hat
  
  
  gamma_hat_alpha <- get_ssc_partition(mst_edges, H, alpha_hat)
  gamma_hat_beta <- get_ssc_partition(mst_edges, H, beta_hat)

  tmprandA <- rand_adj_vi(N, gamma_alpha, gamma_hat_alpha)
  tmprandB <- rand_adj_vi(N, gamma_beta, gamma_hat_beta)

  
  ssc_results[r, "RMSE"] <- sqrt( mean( c( (alpha - alpha_hat)^2, (beta - beta_hat)^2)))
  ssc_results[r, "RAND_A"] <- tmprandA$rand
  ssc_results[r, "RAND_B"] <- tmprandB$rand
  ssc_results[r, "RANDadj_A"] <- tmprandA$rand_adj
  ssc_results[r, "RANDadj_B"] <- tmprandB$rand_adj
  ssc_results[r, "VI_A"] <- tmprandA$vi
  ssc_results[r, "VI_B"] <- tmprandB$vi
  ssc_results[r, "K_A"] <- length(gamma_hat_alpha)
  ssc_results[r, "K_B"] <- length(gamma_hat_beta)
  ssc_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Y_hat)^2 )) 
  
}

assign(paste0("ssc_results_sim", sim_number, "_", batch), ssc_results)
save(list = paste0("ssc_results_sim", sim_number, "_", batch), 
     file = paste0("results/ssc_sim", sim_number, "_", batch, ".RData"))