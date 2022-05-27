library(optimization)
error_safe <- function(expr){
  tryCatch(expr,
           error = function(e){
             message("An error occurred:\n", e)
             NA
           })
}

ensm_givenhyperpar <- function(Y, X, A_block, L, max_iter, 
						hyperpar, gamma_init_A = NULL, gamma_init_B = NULL,
						eta_py, sigma_py = 0.0, lambda = 1.0, rho = 0.9,
						opt_method = 0, opt_Y = 2, 
						last_islands = TRUE, Kmeans_initialize = TRUE,
						priorA = "ep", priorB = "ep",
						A_or_B_first = 1.0){
	if(priorA == "ep"){
		priorA_num = 0
	} else if(priorA == "unif"){
		priorA_num = 3
	} else {
		warning("priorA should be a string `ep` for Ewens-Pitman or `unif` for the uniform.")
	}
	if(priorB == "ep"){
		priorB_num = 0
	} else if(priorB == "unif"){
		priorB_num = 3
	} else {
		warning("priorB should be a string `ep` for Ewens-Pitman or `unif` for the uniform.")
	}
  	source("src/fun_likelihood.R")
	Xorig <- X
	Yorig <- Y
	Xmeans <- rowMeans(X)
	X <- X - Xmeans
	betas_mle <- numeric(N)
	for(i in 1:N)
	betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
	Y <- Y - betas_mle*Xmeans

	om <- opt_method
	oy <- opt_Y

	N <- dim(Y)[1]
	t <- dim(Y)[2]

	K <- round(log(N))
	l <- 1

	Rcpp::sourceCpp('src/main.cpp')
	
	partA <- gamma_init_A
	partB <- gamma_init_B

	a1 <- hyperpar[1]
	a2 <- hyperpar[2]
	b1 <- hyperpar[3]
	b2 <- hyperpar[4]
	alpha_sigma <- hyperpar[5]
	beta_sigma <- hyperpar[6]

	tmp <- ensm_cluster_mean(Y, X, A_block, L = L, 
             lambda_input = lambda, eta_input = eta_py, sigma_py_input = sigma_py, max_iter_input = max_iter, 
             a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
             alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
             opt_method_input = om, opt_Y_input = oy, 
             priorA_input = priorA_num, priorB_input = priorB_num, 
             Kmeans_initialize = Kmeans_initialize, 
             rho_input = rho, A_or_B_first_input = A_or_B_first,
             gamma_init_A = partA, gamma_init_B = partB,
             last_islands = last_islands)

	return(list(output = tmp, hyperpar = hyperpar))
}