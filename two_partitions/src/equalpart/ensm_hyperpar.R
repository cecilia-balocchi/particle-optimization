error_safe <- function(expr){
	tryCatch(expr, error = function(e){
			message("An error occurred:\n", e)
			NA})
}

ensm_hyperpar <- function(Y, X, A_block, L, max_iter, 
						eta_py, sigma_py, lambda, rho = 0.9,
						opt_method = 0, opt_Y = 2, 
						last_islands = TRUE, Kmeans_initialize = TRUE,
						prior = "ep"){
	if(prior == "ep"){
		prior_num = 0
	} else if(prior == "unif"){
		prior_num = 3
	} else {
		warning("prior should be a string `ep` for Ewens-Pitman or `unif` for the uniform.")
	}
  	source("src/fun_likelihood.R")
	X <- X - rowMeans(X)

	N <- dim(Y)[1]
	t <- dim(Y)[2]

	n_tr <- dim(X)[1]
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

	K <- round(log(n_tr))
	
	tmp <- (max(alphas_mle)-min(alphas_mle))/(K+1)/2
	a1 <- tmp^2/sigma2*(1-0.8)
	a2 <- (max(abs(alphas_mle))/2)^2/sigma2 - (a1/(1-rho))

	tmp <- (max(betas_mle)-min(betas_mle))/(K+1)/2
	b1 <- tmp^2/sigma2*(1-0.8)
	b2 <- (max(abs(betas_mle))/2)^2/sigma2 - (b1/(1-rho))

	om <- opt_method
	oy <- opt_Y
	
	l <- 1

	Rcpp::sourceCpp('src/main.cpp')
	tmp <- ensm_cluster_mean(Y_input = Y, X_input = X, A_block_input = A_block, L = l, max_iter_input = max_iter, 
             a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
             alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
             eta_input = eta_py, sigma_py_input = sigma_py, 
             prior_input = prior_num, 
             lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = rho, 
             last_islands = last_islands, Kmeans_initialize = Kmeans_initialize)

	## you need to change the logpost function for when the prior is uniform

	j <- which.max(tmp$w)
	part <- tmp$particle_set_A[[j]]

	log_post2 <- function(par){
	  log_post(par, prior_num, prior_num, part, part, 
	  	a2, b2, rho, Y,X,A_block, alpha_sigma, beta_sigma,
	  	eta_py, sigma_py)
	}
	
	tmp_nm <- error_safe(optim(par = c(a1,b1), fn = log_post2, method = "Nelder-Mead"))
	if(!is.na(tmp_nm[1])){
		a1_new <- tmp_nm$par[1]
		b1_new <- tmp_nm$par[2]
	} else {
		a1_new <- a1
		b1_new <- b1
	}
	

	tmp_new <- ensm_cluster_mean(Y, X, A_block, L = L, max_iter_input = max_iter, 
             a1_input = a1_new, b1_input = b1_new, a2_input = a2, b2_input = b2, 
             alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
             eta_input = eta_py, sigma_py_input = sigma_py, 
             prior_input = prior_num, 
             lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = rho, 
             last_islands = last_islands, Kmeans_initialize = Kmeans_initialize,
             gamma_init = part)

	final_hyperpar <- c(a1_new, a2, b1_new, b2, alpha_sigma, beta_sigma)

	# return(tmp_new)
	return(list(first = tmp, adjusted = tmp_new, optim = tmp_nm, hyperpar = final_hyperpar))
}