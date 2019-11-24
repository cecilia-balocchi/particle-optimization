source("scripts/partition_functions.R")
get_tot_ss <- function(ybar, gamma = list(1:length(ybar))){
  K <- length(gamma)
  tot_ss <- 0
  for(k in 1:K) tot_ss <- tot_ss + sum( (ybar[gamma[[k]]] - mean(ybar[gamma[[k]]]))^2)
  return(tot_ss)
}