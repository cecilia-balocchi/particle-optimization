str <- "_tracts.csv"
A_block <- read.csv(paste0("data/A",str), header = F)
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
### find hyper-parameters
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
###

Rcpp::sourceCpp('src/main.cpp')
l <- 10 			# 10 particles
eta <- 5 			# eta is the parameter of the Ewens-Pitman prior
sigma_py <- 0 		# if we wanted to use a Pitman Yor prior we could set sigma_py \in (0,1)

lambda <- 100 		# strength of the penalty, also tempering parameter
miter <- 500 		# maximum number of iterations
## parameters that slightly change the algorithm
# this parameter decides how the parameters are estimated
om <- 0 			# for orthogonal predictors, set om = 0, it's the fastest.
					# for non-orthogoal pred, set om = 1, it's slower but correct
					# om = 2 gives a faster but approximate estimation
# this parameter decides which likelihood function should be used
oy <- 2				# for orthogonal predictors, set oy = 2 
					# for non-orthogoal pred, set oy = 1
					# oy = 0 is deprecated and not correct.

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
