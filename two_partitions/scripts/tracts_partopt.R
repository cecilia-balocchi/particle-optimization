### this script takes in command line arguments.
### In particular, args was set to c(1, 1, 1)
### If you want to run this same script with "unif-unif prior", set args to c(1, 2, 2)

library(Rcpp)
library(RcppArmadillo)
library(mcclust)

source("src/partition_functions.R")
source("src/ensm_hyperpar.R")
source("src/ensm_givenhyperpar.R")

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

miter <- 500 		# maximum number of iterations
L <- 10 			# 10 particles
eta_py <- 5 		# eta is the parameter of the Ewens-Pitman prior
sigma_py <- 0 		# if we wanted to use a Pitman Yor prior we could set sigma_py \in (0,1)
lambda <- 100 		# strength of the penalty, also tempering parameter

# args <- commandArgs(TRUE)
args <- c(1, 1, 1)
sim_number <- as.numeric(args[1])
priors_strings <- c("ep", "unif") 	# choose between truncated EwensPitman and Uniform (on the space of Spatial Partitions)
priorsAB <- as.numeric(args[2:3])
priors_stringsAB <- priors_strings[priorsAB]


set.seed(1110+sim_number)

z <- Sys.time()
map <- ensm_hyperpar(Y, X, W, L = 1, max_iter = miter, 
                     eta_py = eta_py, sigma_py = sigma_py, 
                     lambda = lambda, rho = 0.9, priorA = priors_stringsAB[1], priorB = priors_stringsAB[2])
time1 <- Sys.time() - z 

## let's set the hyper-parameters for the next ones:
hyperpar <- map$hyperpar # a1, a2, b1, b2, alpha_sigma, beta_sigma

old_map <- map
map <- map$adjusted
gammaA_map <- map$particle_set_A[[1]]
gammaB_map <- map$particle_set_B[[1]]

z <- Sys.time()
lam_100 <- ensm_givenhyperpar(Y, X, W, L = L, max_iter = miter, 
                     hyperpar = hyperpar, gamma_init_A = gammaA_map, gamma_init_B = gammaB_map, 
                     eta_py = eta_py, sigma_py = sigma_py, 
                     lambda = lambda, rho = 0.9, A_or_B_first = 0.5, priorA = priors_stringsAB[1], priorB = priors_stringsAB[2])
time2 <- Sys.time() - z 


filename <- paste0("results/density_L",L,"priors",priors_stringsAB[1],"_",priors_stringsAB[2],"_eta",eta_py,"_lambda",lambda,"task",sim_number,".rdata")
varlist <- c("old_map", "lam_100", "eta_py", "sigma_py", "miter","time1", "time2")
save(list = varlist, file = filename)
