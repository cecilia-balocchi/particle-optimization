### this script takes in command line arguments.
### In particular, args takes values from 1 to 12.

library(Rcpp)
library(RcppArmadillo)
library(mcclust)

source("src/partition_functions.R")
source("src/load_anderson2017.R")

str <- "_density.csv"
A_block <- read.csv(paste0("data/A_block",str), header = F)
Y <- read.csv(paste0("data/Y",str), header = F)
X <- read.csv(paste0("data/X",str), header = F)
E <- read.csv(paste0("data/E",str), header = F)

A_block <- as.matrix(A_block)
colnames(A_block) <- c()
Y <- as.matrix(Y)
colnames(Y) <- c()
X <- as.matrix(X)
colnames(X) <- c()
X <- X - rowMeans(X)
X <- X/apply(X, 1, sd)
E <- as.matrix(E)
colnames(E) <- c()
W <- A_block

C <- sinh(Y+log(2))*E
C <- round(C) # this is the data for running Anderson code

N <- nrow(Y)
T <- ncol(Y)

args <- commandArgs(TRUE)
sim_number <- 1
and_i <- as.numeric(args[1])
anderson_ks <- c("31","33","35","51","53","55","81","83","85","151","153","155")


anderson_ndraws <- 20000


set.seed(1110+sim_number)

and_k <- anderson_ks[and_i]

z <- Sys.time()
if(and_i > 9){
	and <- helper_anderson(C, X, W, n_draws = anderson_ndraws, 
                       num.C = as.numeric(substr(and_k, 1,2)), 
                       num.D = as.numeric(substr(and_k, 3,3)),
                       E = E)
} else {
	and <- helper_anderson(C, X, W, n_draws = anderson_ndraws, 
                   num.C = as.numeric(substr(and_k, 1,1)), 
                   num.D = as.numeric(substr(and_k, 2,2)),
                    E = E)

}
time1 <- Sys.time() - z 


filename <- paste0("results/density_anderson_niter",anderson_ndraws,"ks",and_k,".rdata")
varlist <- c("and", "anderson_ndraws","time1")
save(list = varlist, file = filename)
