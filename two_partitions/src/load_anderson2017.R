# A helper file to load the code from Anderson et al. (2017)
library(MCMCpack)
library(Rcpp)
library(truncdist)
library(mclust)

###Load the R function and Rcpp function to run the model
sourceCpp("src/anderson2017_code/clustfunctions.cpp")
source("src/anderson2017_code/clust_mcmc.r")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

helper_anderson <- function(Y, X, W, n_draws, num.C, num.D,
							E = NULL,
							burnin = NULL){
	
	theta.C <- 1
	theta.D <- 1
	rho <- 0.99
	tau.start <- 0.1
	lambda <- 0.99
	sigma.start <- 0.1
	normal.prior.var <- 10
	gamma.prior.scale <- 0.001
	gamma.prior.shape <- 0.001

	time <- X[1,]
	if(is.null(E)){
		E <- matrix(1, nrow = nrow(Y), ncol = ncol(Y)) # offset is 1 for all values
	}

	tmp_fit <- MCMCfunc(Y = Y, E = E, W = W, n.rep = n_draws, num.C = num.C, num.D = num.D, time = time, 
	                 rho = rho, tau = tau.start, theta.C=theta.C, lambda=lambda, 
	                 sigma=sigma.start, theta.D=theta.D, 
	                 normal.prior.var=normal.prior.var,
	                 gamma.prior.shape=gamma.prior.shape,
	                 gamma.prior.scale=gamma.prior.scale)
	tmp_fit2 <- MCMCfunc(Y = Y, E = E, W = W, n.rep = n_draws, num.C = num.C, num.D = num.D, time = time, 
	                    rho = rho, tau = tau.start, theta.C=theta.C, lambda=lambda, 
	                    sigma=sigma.start, theta.D=theta.D, 
	                    normal.prior.var=normal.prior.var,
	                    gamma.prior.shape=gamma.prior.shape,
	                    gamma.prior.scale=gamma.prior.scale)

	if(is.null(burnin)){
		burnin <- floor(n_draws * 0.5)
	}
	n_tot <- (n_draws - burnin)*2
	index <- (burnin+1):n_draws

	alpha_samples <- matrix(nrow = n_tot, ncol = N)
	beta_samples <- matrix(nrow = n_tot, ncol = N)
	alpha_bar_samples <- matrix(nrow = n_tot, ncol = N)
	beta_bar_samples <- matrix(nrow = n_tot, ncol = N)
	# gamma_alpha_samples <- list()
	# gamma_beta_samples <- list()


	Cstores <- rbind(tmp_fit$C.store[index,],
	                 tmp_fit2$C.store[index,])
	Dstores <- rbind(tmp_fit$D.store[index,],
	                 tmp_fit2$D.store[index,])

	Astores <- rbind(tmp_fit$alpha.store[index,],
	                 tmp_fit2$alpha.store[index,])
	Bstores <- rbind(tmp_fit$beta.store[index,],
	                 tmp_fit2$beta.store[index,])

	Alpha <- colMeans(Astores)
	Beta <- colMeans(Bstores)

	Phi <- colMeans(rbind(tmp_fit$phi.store[index,],tmp_fit2$phi.store[index,]))
	Delta <- colMeans(rbind(tmp_fit$delta.store[index,],tmp_fit2$delta.store[index,]))

	### using the posterior mode of C_i and D_i

	z_C <- apply(Cstores, MARGIN = 2, getmode)
	z_D <- apply(Dstores, MARGIN = 2, getmode)

	unik_z_C <- unique(z_C)
	unik_z_D <- unique(z_D)

	tmp_gamma_alpha <- list()
	for(k in 1:length(unik_z_C)){
	  tmp_gamma_alpha[[k]] <- which(z_C == unik_z_C[k])
	}
	tmp_gamma_beta <- list()
	for(k in 1:length(unik_z_D)){
	  tmp_gamma_beta[[k]] <- which(z_D == unik_z_D[k])
	}

	for(k in 1:length(unik_z_C)){
	  z_C[tmp_gamma_alpha[[k]]] <- as.numeric(names(sort(table(Cstores[,tmp_gamma_alpha[[k]]]), decreasing = TRUE)[1]))
	}
	for(k in 1:length(unik_z_D)){
	  z_D[tmp_gamma_beta[[k]]] <- as.numeric(names(sort(table(Dstores[,tmp_gamma_beta[[k]]]), decreasing = TRUE)[1]))
	}

	Amean <- Alpha[z_C] + Phi
	Bmean <- Beta[z_D] + Delta

	return(list( gamma_alpha = tmp_gamma_alpha, 
               gamma_beta = tmp_gamma_beta, 
               Amean = Amean,
               Bmean = Bmean, 
               Cstores = Cstores,
               Dstores = Dstores,
               fit = tmp_fit, 
               fit2 = tmp_fit2))
}