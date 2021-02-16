# change the value of sim_number based on the cluster separation level desidered:
# sim_number = 1 corresponds to high cluster separation, 2 to moderate, 3 to low.
sim_number = 2 # Figure 4 was done with moderate cluster separation

tmp_data <- data.frame(model = NULL, RMSE = NULL, RAND_A = NULL, RAND_B = NULL, pred_err = NULL, K_A = NULL, K_B = NULL, log.post = NULL)
models <- c("lam_100", "km", "sc", "anderson", "ssc", "fixed_0")

tmp <- load(file = paste0("results/summaries/newX_sim", sim_number, "_results.RData"))

i <- 1
tmp1 <- data.frame(model = rep(models[i], 100), 
                   RMSE = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RMSE"],
                   RAND_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "wavg_RANDadj_A"], 
                   RAND_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "wavg_RANDadj_B"],
                   K_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "wavg_K_A"],
                   K_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "wavg_K_B"],
                   pred_err = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "pred_err"],
                   log.post = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "LP_top"])
tmp_data <- rbind(tmp_data, tmp1)

i <- 4
tmp1 <- data.frame(model = rep(models[i], 100), 
                   RMSE = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RMSE",4],
                   RAND_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "post_RANDadj_A",4], 
                   RAND_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "post_RANDadj_B",4],
                   K_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "post_mod_K_A",4],
                   K_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "post_mod_K_B",4],
                   pred_err = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "pred_err",4],
                   log.post = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "LP",4])
tmp_data <- rbind(tmp_data, tmp1)

i <- 6
tmp1 <- data.frame(model = rep(models[i], 100), 
                   RMSE = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RMSE"],
                   RAND_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_A"], 
                   RAND_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_B"],
                   K_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_A"],
                   K_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_B"],
                   pred_err = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "pred_err"],
                   log.post = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "LP"])
tmp_data <- rbind(tmp_data, tmp1)


tmp <- load(file = paste0("results/summaries/KmSc_sim", sim_number, "_results.RData"))

i <- 2
tmp1 <- data.frame(model = rep(models[i], 100), 
                   RMSE = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RMSE"],
                   RAND_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_A"], 
                   RAND_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_B"],
                   K_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_A"],
                   K_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_B"],
                   pred_err = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "pred_err"],
                   log.post = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "LP"])
tmp_data <- rbind(tmp_data, tmp1)

i <- 3
tmp1 <- data.frame(model = rep(models[i], 100), 
                   RMSE = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RMSE"],
                   RAND_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_A"], 
                   RAND_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_B"],
                   K_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_A"],
                   K_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_B"],
                   pred_err = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "pred_err"],
                   log.post = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "LP"])
tmp_data <- rbind(tmp_data, tmp1)

tmp <- load(file = paste0("results/summaries/ssc_sim", sim_number, ".RData"))
assign(paste0("sim",sim_number,"_ssc_results"), get(paste0("ssc_results_sim",sim_number)))

i <- 5
tmp1 <- data.frame(model = rep(models[i], 100), 
                   RMSE = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RMSE"],
                   RAND_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_A"], 
                   RAND_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "RANDadj_B"],
                   K_A = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_A"],
                   K_B = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "K_B"],
                   pred_err = get(paste0("sim",sim_number,"_",models[i],"_results"))[, "pred_err"],
                   log.post = NA)
tmp_data <- rbind(tmp_data, tmp1)


tmp_data$model[which(tmp_data$model == "anderson")] <- "And"
tmp_data$model[which(tmp_data$model == "lam_100")] <- "PartOpt"
tmp_data$model[which(tmp_data$model == "ssc")] <- "SCC"
tmp_data$model[which(tmp_data$model == "km")] <- "KM"
tmp_data$model[which(tmp_data$model == "sc")] <- "SC"
tmp_data$model[which(tmp_data$model == "fixed_0")] <- "true"

index <- which((tmp_data$model != "SCC") & (tmp_data$model != "true"))
index2 <- which(tmp_data$model != "true")


tmp_data1 <- tmp_data[index, ]
tmp_data1$model <- factor(tmp_data1$model , levels=c("PartOpt", "And", "KM", "SC"))
tmp_data2 <- tmp_data[index2, ]
tmp_data2$model <- factor(tmp_data2$model , levels=c("PartOpt", "And", "KM", "SC", "SCC"))

tmp_data$model <- factor(tmp_data$model , levels=c("PartOpt", "true", "And", "KM", "SC", "SCC"))

## this is Figure 4
png(paste0("figures/boxplots_sim",sim_number,".png"), width = 8, height = 8, units = "in", res = 200)
par(mar = c(3,3,2,1), mfrow = c(2,2), cex.axis = 1.2, cex.main = 1.3)
boxplot(tmp_data1$RMSE ~ tmp_data1$model, main = "RMSE", xlab = "model")
boxplot(tmp_data1$pred_err ~ tmp_data1$model, main = "Prediction error", xlab = "model")
boxplot(tmp_data2$RAND_A ~ tmp_data2$model, main = expression(bold("Rand Index for " ~ gamma^{(alpha)})), xlab = "")
boxplot(tmp_data2$RAND_B ~ tmp_data2$model, main = expression(bold("Rand Index for " ~ gamma^{(beta)})), xlab = "")
dev.off()


# this is Figure S3
png(paste0("figures/boxplots2_sim",sim_number,"_1row.png"), width = 9, height = 3, units = "in", res = 300)
par(mar = c(3,3,2,1), mfrow = c(1,3))
boxplot(tmp_data$log.post ~ tmp_data$model, main = "Log-posterior", xlab = "")
boxplot(tmp_data2$K_A ~ tmp_data2$model, main = expression(bold("Number of clusters for " ~ gamma^{(alpha)})), xlab = ""); abline(h = 10, lty = 2)
boxplot(tmp_data2$K_B ~ tmp_data2$model, main = expression(bold("Number of clusters for " ~ gamma^{(beta)})), xlab = ""); abline(h = 4, lty = 2)
dev.off()
