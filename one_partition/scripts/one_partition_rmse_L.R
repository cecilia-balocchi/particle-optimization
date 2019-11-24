## RMSE as a function of L ##

rmse_lam1 <- matrix(nrow = 6, ncol = 10, dimnames = list(c(), paste0("wtd_RMSE_", 1:10)))
rmse_lam10 <- matrix(nrow = 6, ncol = 10, dimnames = list(c(), paste0("wtd_RMSE_", 1:10)))
rmse_lam100 <- matrix(nrow = 6, ncol = 10, dimnames = list(c(), paste0("wtd_RMSE_", 1:10)))

for(sim_number in 1:6){
  load(paste0("results/sim", sim_number, "_results.RData"))
  rmse_lam1[sim_number,] <- colMeans(get(paste0("sim", sim_number, "_lam_1_results"))[,paste0("wtd_RMSE_",1:10)])
  rmse_lam10[sim_number,] <- colMeans(get(paste0("sim", sim_number, "_lam_10_results"))[,paste0("wtd_RMSE_",1:10)])
  rmse_lam100[sim_number,] <- colMeans(get(paste0("sim", sim_number, "_lam_100_results"))[,paste0("wtd_RMSE_",1:10)])
}


y_lim <- c(0.9 * min(c(rmse_lam1, rmse_lam10, rmse_lam100), na.rm = TRUE), 
           1.1 * max(c(rmse_lam1, rmse_lam10, rmse_lam100), na.rm = TRUE))

load("results/sim1_results.RData")


plot(1, type = "n", xlim = c(0, 10), ylim = range(rmse_lam100[1,], na.rm = TRUE))
points(1:10,rmse_lam1[1,]/, pch = 16)
points(1:10, rmse_lam10[1,]/rmse_lam10[1,10], pch = 15, col = 'blue')
points(1:10, rmse_lam100[1,], pch = 14, col = 'green')


colMeans(sim1_lam_1_results[,paste0("wtd_RMSE_", 1:10)])
