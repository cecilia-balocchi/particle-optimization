#library(optimization)
error_safe <- function(expr){
  tryCatch(expr,
           error = function(e){
             message("An error occurred:\n", e)
             NA
           })
}

PARTOPT <- function(Y, X, A, L = 10, lambda = 10.0, rho = 0.9,
                    max_iter = 25, eta_py = 1, sigma_py = 0.0, 
                    hyperpar = NULL, hyperpar_sigma = NULL, 
                    gamma_init_A = NULL, gamma_init_B = NULL, 
                    last_islands = TRUE, Kmeans_initialize = TRUE,
                    priorA = "ep", priorB = "ep", verbose = TRUE)
{
  opt_method <- 0
  opt_Y <- 2
  if(!priorA %in% c("ep", "unif")){
    print("[PARTOPT]: argument priorA must be equal to one of `ep' or `unif'")
    stop("Invalid argument for priorA")
  } else{
    priorA_num <- ifelse(priorA == "ep", 0, 3)
  }
  if(!priorB %in% c("ep", "unif")){
    print("[PARTOPT]: argument priorB must be equal to one of `ep' or `unif'")
    stop("Invalid argument for priorB")
  } else{
    priorB_num <- ifelse(priorB == "ep", 0, 3)
  }  
  
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
  
  initial_gammas_okay <- FALSE
  if((!is.null(gamma_init_A)) & (!is.null(gamma_init_B))){
    if(all.equal(1:N, sort(unlist(gamma_init_A))) == TRUE){
      if(all.equal(1:N, sort(unlist(gamma_init_B))) == TRUE){
        initial_gammas_okay <- TRUE
      } else {
        warning("[PARTOPT]: invalid gamma_init_B.")
      }
    } else {
      warning("[PARTOPT]: invalid gamma_init_A.")
    }
  }
  if(is.null(gamma_init_A) & (!is.null(gamma_init_B))){
    warning("[PARTOPT]: only one initial partition provided.")
  }
  if((!is.null(gamma_init_A)) & is.null(gamma_init_B)){
    warning("[PARTOPT]: only one initial partition provided.")
  }
  
  hyperpar_okay <- FALSE
  if(!is.null(hyperpar)){
    if(length(hyperpar) == 4){
      if(all(hyperpar > 0)){
        hyperpar_okay <- TRUE
      } else {
        warning("[PARTOPT]: hyperpar should be a vector of positive reals. Using default hyper-parameters.")
      }
    } else {
      warning("[PARTOPT]: hyperpar should be a vector of length 4. Using default hyper-parameters.")
    }
  }
  hyperpar_sigma_okay <- FALSE
  if(!is.null(hyperpar_sigma)){
    if(length(hyperpar_sigma) == 2){
      if(all(hyperpar_sigma > 0)){
        hyperpar_sigma_okay <- TRUE
      } else {
        warning("[PARTOPT]: hyperpar_sigma should be a vector of positive reals. Using default hyper-parameters.")
      }
    } else {
      warning("[PARTOPT]: hyperpar_sigma should be a vector of length 2. Using default hyper-parameters.")
    }
  }
  
  if(hyperpar_okay){
    a1 <- a1_new <- hyperpar[1]
    a2 <- hyperpar[2]
    b1 <- b1_new <- hyperpar[3]
    b2 <- hyperpar[4]

    if(hyperpar_sigma_okay){
      nu_sigma <- hyperpar_sigma[1]
      lambda_sigma <- hyperpar_sigma[2]
      alpha_sigma = nu_sigma/2
      beta_sigma = nu_sigma*lambda_sigma/2
      message_string <- "[PARTOPT]: hyperpar and hyperpar_sigma provided."
    } else {
      sigmas <- numeric(N)
      for(i in 1:N)
        sigmas[i] <- sd(lm(Y[i,]~X[i,])$residuals)
      sigma2 <- mean(sigmas^2)
      
      mu <- mean(sigmas^2)
      v <- var(sigmas^2)
      alpha_sigma <- mu^2/v + 2
      beta_sigma <- mu*(alpha_sigma-1)
      message_string <- "[PARTOPT]: hyperpar provided, hyperpar_sigma not provided (computing default)."
    }
    if(verbose) print(message_string)
    
    tmp_nm = NULL
    l <- 1
    if(initial_gammas_okay){
      # Still search for MAP in case initial partition is not optimal
      if(verbose) print("[PARTOPT]: Computing MAP, starting from the provided initial partitions.")
      tmp <- .ensm_cluster_mean(Y, X, A, L = l, max_iter_input = max_iter, 
                                a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
                                alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
                                eta_input = eta_py, sigma_py_input = sigma_py, 
                                priorA_input = priorA_num, priorB_input = priorB_num,
                                lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = rho, 
                                last_islands = last_islands, Kmeans_initialize = Kmeans_initialize,
                                gamma_init_A = gamma_init_A, gamma_init_B = gamma_init_B, verbose = FALSE)
    } else {
      if(verbose) print("[PARTOPT]: Computing MAP, initial partitions not provided.")
      tmp <- .ensm_cluster_mean(Y, X, A, L = l, max_iter_input = max_iter, 
                                a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
                                alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
                                eta_input = eta_py, sigma_py_input = sigma_py, 
                                priorA_input = priorA_num, priorB_input = priorB_num,
                                lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = rho, 
                                last_islands = last_islands, Kmeans_initialize = Kmeans_initialize, verbose = FALSE)
    }
    j <- 1
    partA <- tmp$particle_set_A[[j]]
    partB <- tmp$particle_set_B[[j]]
    if(verbose) print("[PARTOPT]: Found the MAP. Using that as starting point.")
  } else {
    betas_mle <- numeric(N)
    for(i in 1:N)
      betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
    alphas_mle <- rowMeans(Y) - betas_mle * rowMeans(X) 
    sigmas <- numeric(N)
    for(i in 1:N)
      sigmas[i] <- sd(lm(Y[i,]~X[i,])$residuals)
    sigma2 <- mean(sigmas^2)
    
    if(hyperpar_sigma_okay){
      message_string <- "[PARTOPT]: hyperpar not provided (computing default), hyperpar_sigma provided."
      nu_sigma <- hyperpar_sigma[1]
      lambda_sigma <- hyperpar_sigma[2]
      alpha_sigma = nu_sigma/2
      beta_sigma = nu_sigma*lambda_sigma/2
    } else {
      message_string <- "[PARTOPT]: hyperpar and hyperpar_sigma not provided (computing default)."
      mu <- mean(sigmas^2)
      v <- var(sigmas^2)
      alpha_sigma <- mu^2/v + 2
      beta_sigma <- mu*(alpha_sigma-1)
    }
    if(verbose) print(message_string)
    
    K <- round(log(N))
    
    tmp <- (max(alphas_mle)-min(alphas_mle))/(K+1)/2
    a1 <- tmp^2/sigma2*(1-0.8)
    a2 <- (max(abs(alphas_mle))/2)^2/sigma2 - (a1/(1-rho))
    
    tmp <- (max(betas_mle)-min(betas_mle))/(K+1)/2
    b1 <- tmp^2/sigma2*(1-0.8)
    b2 <- (max(abs(betas_mle))/2)^2/sigma2 - (b1/(1-rho))
    
    l <- 1
    if(initial_gammas_okay){
      # if they specify both gamma_init_A and gamma_init_B
      if(verbose) print("[PARTOPT]: Computing MAP, starting from the provided initial partitions.")
      tmp <- .ensm_cluster_mean(Y, X, A, L = l, max_iter_input = max_iter, 
                                a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
                                alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
                                eta_input = eta_py, sigma_py_input = sigma_py, 
                                priorA_input = priorA_num, priorB_input = priorB_num,
                                lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = rho, 
                                last_islands = last_islands, Kmeans_initialize = Kmeans_initialize,
                                gamma_init_A = gamma_init_A, gamma_init_B = gamma_init_B, verbose = FALSE)
    } else {
      if(verbose) print("[PARTOPT]: Computing MAP, initial partitions not provided.")
      tmp <- .ensm_cluster_mean(Y, X, A, L = l, max_iter_input = max_iter, 
                                a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
                                alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
                                eta_input = eta_py, sigma_py_input = sigma_py, 
                                priorA_input = priorA_num, priorB_input = priorB_num,
                                lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = rho, 
                                last_islands = last_islands, Kmeans_initialize = Kmeans_initialize, verbose = FALSE)
    }
    j <- which.max(tmp$w)
    partA <- tmp$particle_set_A[[j]]
    partB <- tmp$particle_set_B[[j]]
    if(verbose) print("[PARTOPT]: Found the MAP. Using that as starting point.")
    
    log_post2 <- function(par){
      log_post(par, priorA_num, priorB_num, partA, partB, 
               a2, b2, rho, Y,X,A, alpha_sigma, beta_sigma,
               eta_py, sigma_py)
    }
    
    tmp_nm <- error_safe(optim(par = c(a1,b1), fn = log_post2, method = "Nelder-Mead"))
    if(!any(is.na(tmp_nm))){
      a1_new <- tmp_nm$par[1]
      b1_new <- tmp_nm$par[2]
    } else {
      a1_new <- a1
      b1_new <- b1
    }
    if(verbose) print("[PARTOPT]: Found hyperparameters.")
  }
  
  if(verbose) print("[PARTOPT]: Starting particle optimization procedure now.")
  tmp_new <- .ensm_cluster_mean(Y, X, A, L = L, max_iter_input = max_iter, 
                                a1_input = a1_new, b1_input = b1_new, a2_input = a2, b2_input = b2, 
                                alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
                                eta_input = eta_py, sigma_py_input = sigma_py, 
                                priorA_input = priorA_num, priorB_input = priorB_num, 
                                lambda_input = lambda, opt_method_input = om, opt_Y_input = oy, rho_input = rho, 
                                last_islands = last_islands, Kmeans_initialize = FALSE, # we have a valid starting point, no need for KM anymore
                                gamma_init_A = partA, gamma_init_B = partB, verbose = verbose)
  
  nu_sigma = 2*alpha_sigma
  lambda_sigma = beta_sigma/alpha_sigma
  final_hyperpar <- c(a1_new, a2, b1_new, b2, nu_sigma, lambda_sigma)
  names(final_hyperpar) <- c("a1","a2","b1","b2","nu_sigma", "lambda_sigma")

  return(list(final = tmp_new, MAP = tmp, hyperpar = final_hyperpar))
}
