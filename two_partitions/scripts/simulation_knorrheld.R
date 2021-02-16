# setwd("../review/anderson/")

library(Rcpp)
library(RcppArmadillo)
library(mcclust)
library(cluster)
library(spdep)
library(gbdcd)

source("src/partition_functions.R")

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

batch_size <- 10

gbdcd_metrics <- c("RMSE", "pred_err", "TIME",
                    "post_RAND_A","post_RAND_B", "post_RANDadj_A","post_RANDadj_B", "post_VI_A","post_VI_B", "post_K_A","post_K_B",
                    "est_RAND_A","est_RAND_B", "est_RANDadj_A","est_RANDadj_B", "est_VI_A","est_VI_B", "est_K_A","est_K_B")
gbdcd_results <- matrix(nrow = batch_size, ncol = length(gbdcd_metrics), dimnames = list(c(), gbdcd_metrics))


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

tmp <- listw2sn(mat2listw(W))
neighbors <- array(as.numeric(unlist(tmp[,1:2])), dim = c(nrow(tmp), 2))

## todo:
miter <- 100
anderson_ndraws <- 1000
gBDCD_ndraws <- 100000
gBDCD_burnin <- gBDCD_ndraws/2
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
  
  z <- Sys.time()
  out_A <- gaussianBDCD(y = alphas_mle, 
                    neigh = neighbors, 
                    c = 0.35,                 # parameter indicating the a priori number of clusters
                    n_iterations = gBDCD_ndraws, 
                    burn_in = gBDCD_burnin,   # burnin is subtracted by n_iterations ( total reported is n_iterations-burn_in )
                    mu0 = 0,                  # a priori mean
                    sigma0 = sqrt(10))        # a priori standard deviation
  n_cluster <- as.numeric(names(which.max(table(out_A$k.MCMC))))
  zA_gbdcd <- cutree(out_A$cluster.info, k = n_cluster)

  out_B <- gaussianBDCD(y = betas_mle, 
                    neigh = neighbors, 
                    c = 0.35,                 # parameter indicating the a priori number of clusters
                    n_iterations = gBDCD_ndraws, 
                    burn_in = gBDCD_burnin,   # burnin is subtracted by n_iterations ( total reported is n_iterations-burn_in )
                    mu0 = 0,                  # a priori mean
                    sigma0 = sqrt(10))        # a priori standard deviation
  n_cluster <- as.numeric(names(which.max(table(out_B$k.MCMC))))
  zB_gbdcd <- cutree(out_B$cluster.info, k = n_cluster)

  gamma_A_gbdcd <- list(); cont <- 1
  for(x in unique(zA_gbdcd)){
    gamma_A_gbdcd[[cont]] <- which(zA_gbdcd == x)
    cont <- cont + 1
  }
  gamma_B_gbdcd <- list(); cont <- 1
  for(x in unique(zB_gbdcd)){
    gamma_B_gbdcd[[cont]] <- which(zB_gbdcd == x)
    cont <- cont + 1
  }
  
  partitions_A <- array(NA, dim = c(length(out_A$vec.centers), N))
  centers <- strsplit(out_A$vec.centers, split = ";")
  for(i in 1:length(out_A$vec.centers)){
    partitions_A[i,] <- RcppPartition(neighbors, as.numeric(centers[[i]]))
  }
  partitions_B <- array(NA, dim = c(length(out_B$vec.centers), N))
  centers <- strsplit(out_B$vec.centers, split = ";")
  for(i in 1:length(out_B$vec.centers)){
    partitions_B[i,] <- RcppPartition(neighbors, as.numeric(centers[[i]]))
  }

  time1 <- Sys.time() - z 
  gbdcd_results[r, "TIME"] <- time1
  gbdcd_results[r, "RMSE"] <- sqrt(mean( c((alpha - out_A$mean.info[,2])^2,
                                           (beta - out_B$mean.info[,2])^2) ))
  Ypred <- out_A$mean.info[,2] + out_B$mean.info[,2] * t_new_scaled
  gbdcd_results[r, "pred_err"] <- sqrt(mean( (Y_outsample - Ypred)^2 )) 

  tmprandA <- rand_adj_vi(N, gamma_alpha, gamma_A_gbdcd)
  tmprandB <- rand_adj_vi(N, gamma_beta, gamma_B_gbdcd)
  
  gbdcd_results[r, "est_RAND_A"] <- tmprandA$rand
  gbdcd_results[r, "est_RAND_B"] <- tmprandB$rand
  gbdcd_results[r, "est_RANDadj_A"] <- tmprandA$rand_adj
  gbdcd_results[r, "est_RANDadj_B"] <- tmprandB$rand_adj
  gbdcd_results[r, "est_VI_A"] <- tmprandA$vi
  gbdcd_results[r, "est_VI_B"] <- tmprandB$vi
  gbdcd_results[r, "est_K_A"] <- length(gamma_A_gbdcd)
  gbdcd_results[r, "est_K_B"] <- length(gamma_B_gbdcd)

  gbdcd_results[r, "post_RAND_A"] <- mean(apply(partitions_A, MARGIN = 1, FUN = function(x) arandi(x, zA, adjust = FALSE)))
  gbdcd_results[r, "post_RAND_B"] <- mean(apply(partitions_B, MARGIN = 1, FUN = function(x) arandi(x, zB, adjust = FALSE)))
  gbdcd_results[r, "post_RANDadj_A"] <- mean(apply(partitions_A, MARGIN = 1, FUN = function(x) arandi(x, zA)))
  gbdcd_results[r, "post_RANDadj_B"] <- mean(apply(partitions_B, MARGIN = 1, FUN = function(x) arandi(x, zB)))
  gbdcd_results[r, "post_VI_A"] <- mean(apply(partitions_A, MARGIN = 1, FUN = function(x) vi.dist(x, zA)))
  gbdcd_results[r, "post_VI_B"] <- mean(apply(partitions_B, MARGIN = 1, FUN = function(x) vi.dist(x, zB)))
  gbdcd_results[r, "post_K_A"] <- mean(out_A$k.MCMC)
  gbdcd_results[r, "post_K_B"] <- mean(out_B$k.MCMC)

}

assign(paste0("gbdcd_results_sim", sim_number, "_",batch), gbdcd_results)


save_list <- paste0(c("gbdcd_"), "results_sim", sim_number, "_", batch)
# save_list <- c(save_list, save_list1)
save(list = save_list, file = paste0("results/gbdcd_sim", sim_number, "_", batch, ".RData"))
