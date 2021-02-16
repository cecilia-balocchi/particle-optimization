# Compile the results for sims 1 -- 3
rm(list = ls())

n_batches <- c(25, 100, 25)
batch_sizes <- c(4, 1, 4)

particle_metrics <- c("RAND_A_top", "RAND_B_top","RANDadj_A_top", "RANDadj_B_top", "VI_A_top", "VI_B_top", "LP_top", 
                      "avg_RAND_A", "avg_RAND_B","avg_RANDadj_A", "avg_RANDadj_B","avg_VI_A", "avg_VI_B",
                      "wavg_RAND_A", "wavg_RAND_B","wavg_RANDadj_A", "wavg_RANDadj_B","wavg_VI_A", "wavg_VI_B",
                      "RMSE","RMSE_top", "pred_err", "pred_err_top",
                      "K_A_top","K_B_top","avg_K_A","avg_K_B","wavg_K_A","wavg_K_B",
                      "p_trueA","p_trueB","p_trueAB",
                      "L", "TIME")

km_sc_metrics <- c("RMSE","RMSE_mle", "K_A", "K_B",
                   "RAND_A", "RAND_B","RANDadj_A", "RANDadj_B", "VI_A", "VI_B", "pred_err","LP","TIME")
km_metrics <- km_sc_metrics
sc_metrics <- km_sc_metrics

fixed_metrics <- c("RMSE", "RAND_A","RAND_B","RANDadj_A", "RANDadj_B", "VI_A", "VI_B", "K_A","K_B","pred_err", "LP")
anderson_metrics <- c("RMSE", 
                      "post_RAND_A","post_RAND_B", "post_RANDadj_A","post_RANDadj_B", "post_VI_A","post_VI_B", "post_K_A","post_K_B",
                      "post_mod_RAND_A","post_mod_RAND_B", "post_mod_RANDadj_A","post_mod_RANDadj_B", "post_mod_VI_A","post_mod_VI_B", "post_mod_K_A","post_mod_K_B",
                      "est_RAND_A","est_RAND_B", "est_RANDadj_A","est_RANDadj_B", "est_VI_A","est_VI_B", "est_K_A","est_K_B",
                      "est_mod_RAND_A","est_mod_RAND_B", "est_mod_RANDadj_A","est_mod_RANDadj_B", "est_mod_VI_A","est_mod_VI_B", "est_mod_K_A","est_mod_K_B",
                      "pred_err", "LP", "TIME")
anderson_ks <- c("33","35","53","55")
# (1 + (batch-1)*batch_size) : (batch_size*batch)
# batch = 1    1:batch_size
# batch = 2    (batch_size + 1):(2*batch_size)
# batch = 3    (2*batch_size + 1):(3*batch_size)

for(sim_number in 1:3){
  n_batch <- n_batches[sim_number]
  batch_size <- batch_sizes[sim_number]
  tot_batches <- n_batch * batch_size
  
  lam_1_results <- matrix(nrow = tot_batches, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  lam_10_results <- matrix(nrow = tot_batches, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  lam_100_results <- matrix(nrow = tot_batches, ncol = length(particle_metrics), dimnames = list(c(), particle_metrics))
  
  km_results <- matrix(nrow = tot_batches, ncol = length(km_metrics), dimnames = list(c(), km_metrics))
  sc_results <- matrix(nrow = tot_batches, ncol = length(sc_metrics), dimnames = list(c(), sc_metrics))
  
  anderson_results <- array(NA, dim = c(tot_batches,length(anderson_metrics),length(anderson_ks)), dimnames = list(c(), anderson_metrics, anderson_ks))
  
  fixed_0_results <- matrix(nrow = tot_batches, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  fixed_1_results <- matrix(nrow = tot_batches, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  fixed_n_results <- matrix(nrow = tot_batches, ncol = length(fixed_metrics), dimnames = list(c(), fixed_metrics))
  map_results <- matrix(nrow = tot_batches, ncol = length(fixed_metrics) + 1, dimnames = list(c(), c(fixed_metrics, "TIME")))
  
  
  for(batch in 1:n_batch){
    str = paste0("results/new2newX_sim", sim_number, "_", batch, ".RData")
    if(file.exists(str)){
      load(str)
    } else {
      cat(sim_number, " ", batch, "\n")
      next
    }
    
    lam_1_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("lam_1_results_sim", sim_number, "_",batch))
    rm(list = paste0("lam_1_results_sim", sim_number, "_",batch))
    
    lam_10_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("lam_10_results_sim", sim_number, "_",batch))
    rm(list = paste0("lam_10_results_sim", sim_number, "_",batch))
    
    lam_100_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("lam_100_results_sim", sim_number, "_",batch))
    rm(list = paste0("lam_100_results_sim", sim_number, "_",batch))
    
    km_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("km_results_sim", sim_number, "_", batch))
    rm(list = paste0("km_results_sim", sim_number, "_",batch))
    
    sc_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("sc_results_sim", sim_number, "_", batch))
    rm(list = paste0("sc_results_sim", sim_number, "_",batch))

    anderson_results[(1 + (batch-1)*batch_size) : (batch_size*batch),,] <- get(paste0("anderson_results_sim", sim_number, "_",batch))
    rm(list = paste0("anderson_results_sim", sim_number, "_",batch))
    
    fixed_0_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("fixed_0_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_0_results_sim", sim_number, "_",batch))
    
    fixed_1_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("fixed_1_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_1_results_sim", sim_number, "_",batch))
    
    fixed_n_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("fixed_n_results_sim", sim_number, "_", batch))
    rm(list = paste0("fixed_n_results_sim", sim_number, "_",batch))
    
    map_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("map_results_sim", sim_number, "_", batch))
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
  
  assign(paste0("sim", sim_number, "_sc_results"), sc_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_sc_results"))
  
  assign(paste0("sim", sim_number, "_fixed_0_results"), fixed_0_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_0_results"))
  
  assign(paste0("sim", sim_number, "_fixed_1_results"), fixed_1_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_1_results"))
  
  assign(paste0("sim", sim_number, "_fixed_n_results"), fixed_n_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_fixed_n_results"))
  
  assign(paste0("sim", sim_number, "_map_results"), map_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_map_results"))
  
  assign(paste0("sim", sim_number, "_anderson_results"), anderson_results)
  save_list <- c(save_list, paste0("sim", sim_number, "_anderson_results"))
  
  save(list = save_list, file = paste0("results/summaries/newX_sim", sim_number, "_results.RData"))
}

############################################
############################################
#### I actually ran the simulations for Km and Sc in a separate script, so I also compiled them separately
#### and created a separate summary file: paste0("KmSc_sim", sim_number, "_results.RData")
#### If you run the simulations using simulation.R it should be sufficient to compile the results with the code above

# # Compile the results for sims 1 -- 3
# rm(list = ls())

# n_batches <- c(10, 10, 10)
# batch_sizes <- c(10, 10, 10)

# km_sc_metrics <- c("RMSE","RMSE_mle", "K_A", "K_B",
#                    "RAND_A", "RAND_B","RANDadj_A", "RANDadj_B", "VI_A", "VI_B", "pred_err","LP","TIME")
# km_metrics <- km_sc_metrics
# sc_metrics <- km_sc_metrics

# for(sim_number in 1:3){
#   n_batch <- n_batches[sim_number]
#   batch_size <- batch_sizes[sim_number]
#   tot_batches <- n_batch * batch_size
  
#   km_results <- matrix(nrow = tot_batches, ncol = length(km_metrics), dimnames = list(c(), km_metrics))
#   sc_results <- matrix(nrow = tot_batches, ncol = length(sc_metrics), dimnames = list(c(), sc_metrics))
  
#   for(batch in 1:n_batch){
#     str = paste0("results/KmSc_sim", sim_number, "_", batch, ".RData")
#     if(file.exists(str)){
#       load(str)
#     } else {
#       cat(sim_number, " ", batch, "\n")
#       next
#     }
    
#     km_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("km_results_sim", sim_number, "_", batch))
#     rm(list = paste0("km_results_sim", sim_number, "_",batch))
    
#     sc_results[(1 + (batch-1)*batch_size) : (batch_size*batch),] <- get(paste0("sc_results_sim", sim_number, "_", batch))
#     rm(list = paste0("sc_results_sim", sim_number, "_",batch))
    
#   }
#   save_list <- c()
  
#   assign(paste0("sim", sim_number, "_km_results"), km_results)
#   save_list <- c(save_list, paste0("sim", sim_number, "_km_results"))
  
#   assign(paste0("sim", sim_number, "_sc_results"), sc_results)
#   save_list <- c(save_list, paste0("sim", sim_number, "_sc_results"))
  
#   save(list = save_list, file = paste0("results/summaries/KmSc_sim", sim_number, "_results.RData"))
# }
