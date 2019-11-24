load("data/crime_tract2019.rdata")
Y <- read.csv("data/Y_tracts.csv", header = F)
X <- read.csv("data/X_tracts.csv", header = F)
Y <- as.matrix(Y)
X <- as.matrix(X)
x_mean <- mean(X[1,])
X <- X - rowMeans(X)
n_tr <- nrow(Y)

n_year <- 2018
index_da <- crime_tract2019$year!= n_year 

outofsample = crime_tract2019$tr.violent[crime_tract2019$year== 2018]
t_new <- 13 - x_mean

strings <- c("results/EPPrior5_0_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata",
             "results/UnifPrior_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata",
             "results/UnifEPPrior5_0_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata",
             "results/EPUnifPrior5_0_tracts2_L10KMInitializeIS_opts02newhyper2_K6lambda100islandFIXED_kmrep.rdata")

RMSE <- array(NA, dim = c(4,2))

# Get MLE values
betas.mle <- numeric(n_tr)
for(i in 1:n_tr)
  betas.mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
alphas.mle <- rowMeans(Y)
predicted.mle = alphas.mle + t_new*betas.mle
RMSE.mle <- sqrt(mean((outofsample-predicted.mle)^2))


for(i in 1:length(strings)){
  load(strings[i]) # the main object saved is called tmp
  # Top particle
  i1 <- which.max(tmp$log_posterior)
  alpha.est <- tmp$alpha_particle[,i1]
  beta.est <- tmp$beta_particle[,i1]
  predicted.est = alpha.est + t_new*beta.est
  RMSE[i,1] <- sqrt(mean((outofsample-predicted.est)^2))
  # BMA
  alpha.estL <- alpha.est*0
  beta.estL <- beta.est*0
  for(j in 1:length(tmp$w)){
    alpha.estL <- alpha.estL + tmp$alpha_particle[,j] * tmp$w[j]
    beta.estL <- beta.estL + tmp$beta_particle[,j] * tmp$w[j]
  }
  predicted.estL = alpha.estL + t_new*beta.estL
  RMSE[i,2] <- sqrt(mean((outofsample-predicted.estL)^2))
}

