library(Rcpp)
library(RcppArmadillo)
library(mcclust)
library(cluster)

source("src/partition_functions.R")
source("src/ensm_hyperpar.R")
source("src/summary_givenhyperpar.R")

str <- "_density.csv"
A_block <- read.csv(paste0("data/A_block",str), header = F)
Y <- read.csv(paste0("data/Y",str), header = F)
X <- read.csv(paste0("data/X",str), header = F)

A_block <- as.matrix(A_block)
colnames(A_block) <- c()
Y <- as.matrix(Y)
colnames(Y) <- c()
X <- as.matrix(X)
colnames(X) <- c()
X <- X - rowMeans(X)
X <- X/apply(X, 1, sd)
W <- A_block

N <- nrow(Y)
T <- ncol(Y)

KMsplits <- 30
SCsplits <- 20

miter <- 500
eta_py <- 1
sigma_py <- 0
lambda <- 100

sim_number <- 1

betas_mle <- numeric(N)
for(i in 1:N)
	betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
alphas_mle <- rowMeans(Y)

set.seed(1110+sim_number)

z <- Sys.time()
map <- ensm_hyperpar(Y, X, W, L = 1, max_iter = miter, 
                     eta_py = eta_py, sigma_py = sigma_py, 
                     lambda = lambda, rho = 0.9, priorA = "ep", priorB = "ep")
time1 <- Sys.time() - z 

## let's set the hyper-parameters for the next ones:
hyperpar <- map$hyperpar 

varlist <- c("map", "hyperpar")
#########

z <- Sys.time()
km_silhouetteA <- numeric(KMsplits)
km_silhouetteB <- numeric(KMsplits)
dA <- dist(alphas_mle)
dB <- dist(alphas_mle)
for(k in 2:KMsplits){
	km <- kmeans(alphas_mle, centers = k, nstart = 25)
	tmp_sil <- silhouette(x = km$cluster, dist = dA)
	km_silhouetteA[k] <- mean(tmp_sil[, 3])

	km <- kmeans(betas_mle, centers = k, nstart = 25)
	tmp_sil <- silhouette(x = km$cluster, dist = dB)
	km_silhouetteB[k] <- mean(tmp_sil[, 3])
}
kA <- which.max(km_silhouetteA)
kB <- which.max(km_silhouetteB)

kmA <- kmeans(alphas_mle, centers = kA, nstart = 25)
kmB <- kmeans(betas_mle, centers = kB, nstart = 25)
time1 <- Sys.time() - z 

varlist <- c(varlist, "kmA", "kmB")

zA <- Modify_partition(kmA$cluster, W)
zB <- Modify_partition(kmB$cluster, W)

gammaA_top_km <- list(); count <- 1
for(k in unique(zA)){
	gammaA_top_km[[count]] <- which(zA == k)
	count <- count + 1
}
gammaB_top_km <- list(); count <- 1
for(k in unique(zB)){
	gammaB_top_km[[count]] <- which(zB == k)
	count <- count + 1
}

tmp_km <- summary_givenhyperpar(Y, X, W, 
                          hyperpar = hyperpar,
                          gamma_init_A = gammaA_top_km, gamma_init_B = gammaB_top_km,
                          eta_input = 1.0, rho_input = 0.9)

varlist <- c(varlist, "gammaA_top_km", "gammaB_top_km", "tmp_km")
#########

z <- Sys.time()
sc_silhouetteA <- numeric(SCsplits)
sc_silhouetteB <- numeric(SCsplits)
for(k in 2:SCsplits){
	sc <- spectral_clustering(alphas_mle, W, k)
	tmp_sil <- silhouette(x = sc$cluster, dist = dA)
	sc_silhouetteA[k] <- mean(tmp_sil[, 3])

	sc <- spectral_clustering(betas_mle, W, k)
	tmp_sil <- silhouette(x = sc$cluster, dist = dB)
	sc_silhouetteB[k] <- mean(tmp_sil[, 3])
}
kA <- which.max(sc_silhouetteA)
kB <- which.max(sc_silhouetteB)

scA <- spectral_clustering(alphas_mle, W, kA)
scB <- spectral_clustering(betas_mle, W, kB)
time2 <- Sys.time() - z 

varlist <- c(varlist, "scA", "scB")

zA <- Modify_partition(scA$cluster, W)
zB <- Modify_partition(scB$cluster, W)

gammaA_top_sc <- list(); count <- 1
for(k in unique(zA)){
	gammaA_top_sc[[count]] <- which(zA == k)
	count <- count + 1
}
gammaB_top_sc <- list(); count <- 1
for(k in unique(zB)){
	gammaB_top_sc[[count]] <- which(zB == k)
	count <- count + 1
}

tmp_sc <- summary_givenhyperpar(Y, X, W, 
                          hyperpar = hyperpar,
                          gamma_init_A = gammaA_top_sc, gamma_init_B = gammaB_top_sc,
                          eta_input = 1.0, rho_input = 0.9)
  
varlist <- c(varlist, "gammaA_top_sc", "gammaB_top_sc", "tmp_sc")

filename <- paste0("results/density_KmSc.rdata")
varlist <- c(varlist, "time1", "time2")
save(list = varlist, file = filename)
