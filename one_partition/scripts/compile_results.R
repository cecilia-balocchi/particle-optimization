# Compile the results for sims 1 -- 3




batch_size <- 20

particle_metrics <- c("RMSE", "L", "RMSE_top", "RAND_top", "K_top", "LP_top", paste0("RMSE_", 1:10), paste0("wtd_RMSE_", 1:10), paste0("LP_", 1:10), 
                      paste0("K_", 1:10), paste0("avg_K_", 1:10), paste0("RAND_", 1:10), paste0("avg_RAND_", 1:10),"TIME")

km_metrics <- c(paste0("RMSE_", 1:10), paste0("K_", 1:10), paste0("RAND_", 1:10), paste0("LP_", 1:10), "TIME")

sc_metrics <- c(paste0("RMSE_", 1:10), paste0("K_", 1:10), paste0("RAND_", 1:10), paste0("LP_", 1:10), "TIME")

fixed_metrics <- c("RMSE", "RAND", "K", "LP")


for(sim_number in 1:3){
  lam_1_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  lam_10_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  
  lam_100_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  km_results <- matrix(nrow = batch_size, ncol = length(km_metrics), dimnames = list(c(), km_metrics))
  km_ss <- matrix(nrow = batch_size, ncol = 10)
  
  sc_results <- matrix(nrow = batch_size, ncol = length(sc_metrics), dimnames = list(c(), sc_metrics))
  sc_ss <- matrix(nrow = batch_size, ncol  = 10)
  
  
  fixed_0_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  fixed_1_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  fixed_n_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  map_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics) + 1, dimnames = list(c(), c(fixed_metrics, "TIME")))
  
  
  for(batch in 1:5){
    load(paste0("results/sim", sim_number, "_", batch, ".RData"))
    
    lam_1_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("lam_1_results_sim", sim_number, "_",batch))
    rm(list = paste0("lam_1_results_sim", sim_number, "_",batch))
    
    lam_10_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("lam_10_results_sim", sim_number, "_",batch))
    rm(list = paste0("lam_10_results_sim", sim_number, "_",batch))
    
    lam_100_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("lam_100_results_sim", sim_number, "_",batch))
    rm(list = paste0("lam_100_results_sim", sim_number, "_",batch))
    
    km_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("km_results_sim", sim_number, "_", batch))
    rm(list = paste0("km_results_sim", sim_number, "_",batch))
    
    km_ss[(1 + (batch-1)*4):(4*batch),] <- get(paste0("km_ss_sim", sim_number, "_", batch))
    rm(list = paste0("km_ss_sim", sim_number, "_",batch))
    
    sc_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("sc_results_sim", sim_number, "_", batch))
    rm(list = paste0("sc_results_sim", sim_number, "_",batch))
    sc_ss[(1 + (batch-1)*4):(4*batch),] <- get(paste0("sc_ss_sim", sim_number, "_", batch))
    rm(list = paste0("sc_ss_sim", sim_number, "_",batch))
    
    fixed_0_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("fixed_0_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_0_results_sim", sim_number, "_",batch))
    
    fixed_1_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("fixed_1_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_1_results_sim", sim_number, "_",batch))
    
    fixed_n_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("fixed_n_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_n_results_sim", sim_number, "_",batch))
    
    map_results[(1 + (batch-1)*4):(4*batch),] <- get(paste0("map_results_sim", sim_number, "_", batch))
    rm(list = paste0("map_results_sim", sim_number, "_",batch))
  }
  save_list <- c()
  assign(paste0("sim", sim_number, "_lam_1_results"), lam_1_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_lam_1_results"))
  
  assign(paste0("sim", sim_number, "_lam_10_results"), lam_10_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_lam_10_results"))
  
  assign(paste0("sim", sim_number, "_lam_100_results"), lam_100_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_lam_100_results"))
  
  assign(paste0("sim", sim_number, "_km_results"), km_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_km_results"))
  assign(paste0("sim", sim_number, "_km_ss"), km_ss)
  save_list <- c(save_list, paste0("sim", sim_number, "_km_ss"))
  
  assign(paste0("sim", sim_number, "_sc_results"), sc_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_sc_results"))
  assign(paste0("sim", sim_number, "_sc_ss"), sc_ss)
  save_list <- c(save_list, paste0("sim", sim_number, "_sc_ss"))
  
  assign(paste0("sim", sim_number, "_fixed_0_results"), fixed_0_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_0_results"))
  
  assign(paste0("sim", sim_number, "_fixed_1_results"), fixed_1_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_1_results"))
  
  assign(paste0("sim", sim_number, "_fixed_n_results"), fixed_n_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_n_results"))
  
  assign(paste0("sim", sim_number, "_map_results"), map_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_map_results"))
  
  save(list = save_list, file = paste0("results/sim", sim_number, "_results.RData"))
}

for(sim_number in 4:6){
  lam_1_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  lam_10_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  
  lam_100_results <- matrix(nrow = batch_size, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  km_results <- matrix(nrow = batch_size, ncol = length(km_metrics), dimnames = list(c(), km_metrics))
  km_ss <- matrix(nrow = batch_size, ncol = 10)
  
  sc_results <- matrix(nrow = batch_size, ncol = length(sc_metrics), dimnames = list(c(), sc_metrics))
  sc_ss <- matrix(nrow = batch_size, ncol  = 10)
  
  
  fixed_0_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  fixed_1_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  fixed_n_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  map_results <- matrix(nrow = batch_size, ncol = length(fixed_metrics) + 1, dimnames = list(c(), c(fixed_metrics, "TIME")))
  
  
  for(batch in 1:10){
    load(paste0("results/sim", sim_number, "_lam_1_", batch, ".RData"))
    lam_1_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("lam_1_results_sim", sim_number, "_",batch))[,colnames(lam_1_results)]
    rm(list = paste0("lam_1_results_sim", sim_number, "_",batch))
    
    load(paste0("results/sim", sim_number, "_lam_10_", batch, ".RData"))
    lam_10_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("lam_10_results_sim", sim_number, "_",batch))[,colnames(lam_10_results)]
    rm(list = paste0("lam_10_results_sim", sim_number, "_",batch))
    
    load(paste0("results/sim", sim_number, "_lam_100_", batch, ".RData"))
    lam_100_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("lam_100_results_sim", sim_number, "_",batch))[,colnames(lam_100_results)]
    rm(list = paste0("lam_100_results_sim", sim_number, "_",batch))
    
    load(paste0("results/sim", sim_number, "_fixed_results_", batch, ".RData"))
    
    km_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("km_results_sim", sim_number, "_", batch))
    rm(list = paste0("km_results_sim", sim_number, "_",batch))
    
    km_ss[(1 + (batch-1)*2):(2*batch),] <- get(paste0("km_ss_sim", sim_number, "_", batch))
    rm(list = paste0("km_ss_sim", sim_number, "_",batch))
    
    sc_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("sc_results_sim", sim_number, "_", batch))
    rm(list = paste0("sc_results_sim", sim_number, "_",batch))
    sc_ss[(1 + (batch-1)*2):(2*batch),] <- get(paste0("sc_ss_sim", sim_number, "_", batch))
    rm(list = paste0("sc_ss_sim", sim_number, "_",batch))
    
    fixed_0_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("fixed_0_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_0_results_sim", sim_number, "_",batch))
    
    fixed_1_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("fixed_1_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_1_results_sim", sim_number, "_",batch))
    
    fixed_n_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("fixed_n_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_n_results_sim", sim_number, "_",batch))
    
    load(paste0("results/sim", sim_number, "_map_",batch, ".RData"))
    map_results[(1 + (batch-1)*2):(2*batch),] <- get(paste0("map_results_sim", sim_number, "_", batch))
    rm(list = paste0("map_results_sim", sim_number, "_",batch))
  }
  save_list <- c()
  assign(paste0("sim", sim_number, "_lam_1_results"), lam_1_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_lam_1_results"))
  
  assign(paste0("sim", sim_number, "_lam_10_results"), lam_10_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_lam_10_results"))
  
  assign(paste0("sim", sim_number, "_lam_100_results"), lam_100_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_lam_100_results"))
  
  assign(paste0("sim", sim_number, "_km_results"), km_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_km_results"))
  assign(paste0("sim", sim_number, "_km_ss"), km_ss)
  save_list <- c(save_list, paste0("sim", sim_number, "_km_ss"))
  
  assign(paste0("sim", sim_number, "_sc_results"), sc_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_sc_results"))
  assign(paste0("sim", sim_number, "_sc_ss"), sc_ss)
  save_list <- c(save_list, paste0("sim", sim_number, "_sc_ss"))
  
  assign(paste0("sim", sim_number, "_fixed_0_results"), fixed_0_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_0_results"))
  
  assign(paste0("sim", sim_number, "_fixed_1_results"), fixed_1_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_1_results"))
  
  assign(paste0("sim", sim_number, "_fixed_n_results"), fixed_n_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_n_results"))
  
  assign(paste0("sim", sim_number, "_map_results"), map_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_map_results"))
  
  save(list = save_list, file = paste0("results/sim", sim_number, "_results.RData"))
}




