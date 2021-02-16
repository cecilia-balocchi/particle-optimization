library(optimization)
error_safe <- function(expr){
  tryCatch(expr,
           error = function(e){
             message("An error occurred:\n", e)
             NA
           })
}

ensm_givenhyperpar <- function(Y, X, A_block, L, max_iter, 
						hyperpar, gamma_init = NULL, 
						eta_py, sigma_py = 0.0, lambda = 1.0, rho = 0.9,
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

	om <- opt_method
	oy <- opt_Y

	N <- dim(Y)[1]
	t <- dim(Y)[2]

	Rcpp::sourceCpp('src/main.cpp')
	
	part <- gamma_init

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
             prior_input = prior_num, 
             Kmeans_initialize = Kmeans_initialize, 
             rho_input = rho, 
             gamma_init = part, 
             last_islands = last_islands)

	return(list(output = tmp, hyperpar = hyperpar))
}