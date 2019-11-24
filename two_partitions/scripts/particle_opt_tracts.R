str <- "_tracts.csv"
A_block <- read.csv(paste0("data/A_block",str), header = F)
A_block <- as.matrix(A_block)
colnames(A_block) <- c()
Y <- read.csv(paste0("data/Y",str), header = F)
Y <- as.matrix(Y)
colnames(Y) <- c()
X <- read.csv(paste0("data/X",str), header = F)
X <- as.matrix(X)
colnames(X) <- c()
X <- X - rowMeans(X)

n_tr <- dim(X)[1]
# find hyper-parameters
betas_mle <- numeric(n_tr)
for(i in 1:n_tr)
  betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
alphas_mle <- rowMeans(Y)
sigmas <- numeric(n_tr)
for(i in 1:n_tr)
  sigmas[i] <- sd(lm(Y[i,]~X[i,])$residuals)
sigma2 <- mean(sigmas^2)
mu <- mean(sigmas^2)
v <- var(sigmas^2)
alpha_sigma <- mu^2/v + 2
beta_sigma <- mu*(alpha_sigma-1)
K <- 6 
rho <- 0.8
tmp <- (max(alphas_mle)-min(alphas_mle))/(K+1)/2
a1 <- tmp^2/sigma2*(1-0.8)
a2 <- (max(abs(alphas_mle))/2)^2/sigma2 - (a1/(1-rho))
tmp <- (max(betas_mle)-min(betas_mle))/(K+1)/2
b1 <- tmp^2/sigma2*(1-0.8)
b2 <- (max(abs(betas_mle))/2)^2/sigma2 - (b1/(1-rho))

Rcpp::sourceCpp('src/main.cpp')
l <- 10
eta <- 5
sigma_py <- 0 

lambda <- 100
miter <- 500
om <- 0
oy <- 2

# to change the prior, set priorA_input, priorA_input to 0 for EP and 3 for Uniform
tmp <- ensm_cluster_mean(Y, X, A_block, L = l, max_iter_input = miter, 
						 a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
						 alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
						 eta_input = eta, sigma_py_input = sigma_py, priorA_input = 0, priorB_input = 0, 
                         lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = 0.8, last_islands = TRUE,
                         Kmeans_initialize = TRUE)
str(tmp)
## This is an example of my very complicated way to save results
filename <- paste0("results/EPPrior",eta,"_",sigma_py,"_tracts2_L",l,"KMInitializeIS_opts",om,oy,"newhyper2_K",K,"lambda",lambda,".rdata")
varlist <- c("tmp", "a1", "b1", "a2", "b2","K", "eta", "sigma_py")
save(list = varlist, file = filename)
