load("data/shapefile/phillytracts")
load("data/phillygooglemap.rdata")
load("data/crime_tract2019.rdata")

str <- "_density.csv"
Y <- read.csv(paste0("data/Y",str), header = F)
X <- read.csv(paste0("data/X",str), header = F)
E <- read.csv(paste0("data/E",str), header = F)
Y <- as.matrix(Y)
X <- as.matrix(X)
x_mean <- mean(X[1,])
x_sd <- sd(X[1,])
X <- X - rowMeans(X)
X <- X/apply(X, 1, sd)
E <- as.matrix(E)
colnames(E) <- c()

n_tr <- nrow(Y)

C <- sinh(Y+log(2))*E
C <- round(C) # this is the data used for running Anderson code

n_year <- 2018

# ihs_transform <- function(x){
#   return( log(x + sqrt(x^2 + 1)) - log(2))
# }

outofsample = crime_tract2019$tr.density[crime_tract2019$year== 2018]
outofsample_violent = crime_tract2019$violent[crime_tract2019$year== 2018]
outofsample_area = crime_tract2019$ALAND10_sqmile[crime_tract2019$year== 2018]

C_outofsample <- sinh(outofsample+ log(2))*outofsample_area
t_new <- (13 - x_mean)/x_sd

####
## MAP & PartOpt
####

string <- "results/density_L20priorsep_ep_eta1_lambda100task1.rdata"
load(string)

tmp <- lam_100$output

# Top particle (MAP)
i1 <- which.max(tmp$log_posterior)
alpha.est <- tmp$alpha_particle[,i1]
beta.est <- tmp$beta_particle[,i1]
predicted.est = alpha.est + t_new*beta.est
RMSE.map <- sqrt(mean((outofsample-predicted.est)^2))

# BMA (PartOpt)
alpha.estL <- alpha.est*0
beta.estL <- beta.est*0
for(j in 1:length(tmp$w)){
  alpha.estL <- alpha.estL + tmp$alpha_particle[,j] * tmp$w[j]
  beta.estL <- beta.estL + tmp$beta_particle[,j] * tmp$w[j]
}
predicted.estL = alpha.estL + t_new*beta.estL
RMSE.PartOpt <- sqrt(mean((outofsample-predicted.estL)^2))

####
## MLE
####

betas.mle <- numeric(n_tr)
for(i in 1:n_tr)
  betas.mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
alphas.mle <- rowMeans(Y)
predicted.mle = alphas.mle + t_new*betas.mle
predicted.exp.mle = exp(predicted.mle)
RMSE.mle <- sqrt(mean((outofsample-predicted.mle)^2))

####
## Anderson et al (2017)
####

anderson_ks <- c("31","33","35","51","53","55","81","83","85","151","153","155")

RMSE <- numeric(length(anderson_ks))
R2 <- numeric(length(anderson_ks))

for(i in 1:length(anderson_ks)){
  string <- paste0("results/anderson2017/density_anderson_niter20000ks",anderson_ks[i],".rdata")
  if(!file.exists(string)){
    cat("File missing: ",i, "\n")
  } else {
    load(string)
    predicted <- and$Amean + t_new*and$Bmean
    predicted.exp = exp(and$Amean + t_new*and$Bmean)
    RMSE[i] <- sqrt(mean((outofsample-predicted)^2))
    R2[i] <- 1- mean((outofsample-predicted)^2)/mean((outofsample-mean(outofsample))^2)
  }
}

which.min(RMSE); min(RMSE)
RMSE.And <- min(RMSE)

####
## KM & SC
####

string <- paste0("results/density_KmSc.rdata")
tmp <- load(string)

predicted_km <- tmp_km$output$alpha + t_new*tmp_km$output$beta
RMSE.km <- sqrt(mean((outofsample-predicted_km)^2))

predicted_sc <- tmp_sc$output$alpha + t_new*tmp_sc$output$beta
RMSE.sc <- sqrt(mean((outofsample-predicted_sc)^2))

####
## SCC
####

# See scrips/tracts_scc.R
load("esults/ssc_philly_results.RData")
RMSE.scc <- rmse

####################
####### table ######
####################

RMSEs <- c(RMSE.map, RMSE.PartOpt, RMSE.mle, RMSE.And, RMSE.km, RMSE.sc, RMSE.scc)
names(RMSEs) <- c("MAP", "PartOpt", "MLE", "And", "KM", "SC", "SCC")

library(xtable)
xtable(RMSEs, digits = 4)