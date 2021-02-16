library(glmnet)
library(igraph)
source("src/partition_functions.R")

load("data/ssc_philly.RData")

# Run cv.glmnet
set.seed(129)
cv_fit <- cv.glmnet(x = X_tilde, y = Y_train, intercept = FALSE, 
                    penalty.factor = penalty_factor, nfolds = 10, family = "gaussian")
# even though we specify intercept = FALSE
# glmnet still returns an intercept
theta_hat <- coef(cv_fit, s = "lambda.1se")[-1]


tmp_time <- 1:T
t_new <- T+1
t_new_scaled <- (t_new - mean(tmp_time))/sd(tmp_time)

alpha_hat <- H_tilde_inv %*% theta_hat[1:N]
beta_hat <- H_tilde_inv %*% theta_hat[(N+1):(2*N)]

Y_hat <- alpha_hat + t_new_scaled * beta_hat

rmse <- sqrt(mean( (Y_test - Y_hat)^2))
gamma_hat_alpha <- get_ssc_partition(mst_edges, H, alpha_hat)
gamma_hat_beta <- get_ssc_partition(mst_edges, H, beta_hat)


save(Y_hat, rmse, gamma_hat_alpha, gamma_hat_beta, file = "results/ssc_philly_results.RData")
