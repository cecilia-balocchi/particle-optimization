\name{PARTOPT}
\alias{PARTOPT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Particle optimization for spatially clustered regression
}
\description{
Implements the particle optimization method introduced in Balocchi et al. (2022+) for fitting simple linear models with spatially clustered intercepts and spatially clustered slopes.
Full details of the implementation are available at https://arxiv.org/abs/1912.00111
}
\usage{
PARTOPT(Y, X, A, L, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{An n x T matrix of outcomes. Y[i,t] is the t-th observed outcome in spatial unit i. Missing values are not supported.}
  \item{X}{An n x T matrix containing the observed covariates. X[i,t] is the t-th observed covariate in spatial unit i. Missing values are not supported.}
  \item{A}{An n x n symmetric, binary adjaceny matrix for the spatial units. A[i,j] is equal to 1 if spatial units i and j are geographically adjacent. Note, A is assumed to be symmetric and dense.}
  \item{L}{Total number of pairs of partitions (hereafter particles). Default is 10.}
  \item{lambda}{Penalty parameter. Default is 10.}
  \item{rho}{Spatial correlation hyperparameter for the Leroux-style conditionally auto-regressive prior. Must be between 0 and 1. Default is 0.9.}
  \item{max_iter}{Maximum number of iterations. Default is 25.}
  \item{eta_py}{Hyperparameter for the Pitman-Yor prior on partitions. Default is 1.}
  \item{sigma_py}{Hyperparameter for the Pitman-Yor prior on partitions. Default is 0, correpsonding to the Ewens-Pitman distribution.}
  \item{hyperpar}{Vector of length 4 containing the (positive) scale hyperparameters a1, a2, b1, and b2. If NULL, \code{PARTOPT} estimates hyperparameters according to the strategy described in Section S2 of the Supplementary Materials accompanying Balocchi et al. (2022+). Default is NULL}
  \item{hyperpar_sigma}{Vector of length 2 containing the (positive) hyperparameters nu and lambda for the Inverse Gamma(nu/2, nu*lambda/2) prior on the residual variance sigma^2. If NULL, \code{PARTOPT} sets these hyperparameters according to the strategy described in Section 3.2 of Balocchi et al. (2022+). Default is NULL.}
  \item{gamma_init_A}{Optional initial partition of the spatial units according to their corresponding intercepts. If not NULL, gamma_init_A must be a list whose elements are vectors of spatial unit indices. If either gamma_init_A or gamma_init_B is NULL, \code{PARTOPT} will initialize at a maximum a posteriori estimate. Default is NULL.}
  \item{gamma_init_B}{Optional initial partition of the spatial units according to their corresponding slopes. If not NULL, gamma_init_B must be a list whose elements are vectors of spatial unit indices. If either gamma_init_A or gamma_init_B is NULL, \code{PARTOPT} will initialize at a maximum a posteriori estimate. Default is NULL.}
  \item{last_islands}{Logical. If TRUE, after the main loop terminates, \code{PARTOPT} performs a final refinement step that checks whether the objective can be further improved by changing cluster assignment of a single spatial unit. Default is TRUE.}
  \item{Kmeans_initialize}{Logical. If TRUE, initializes the MAP estimation step using a procedure based on K-means clustering. Default is TRUE.}
  \item{prior_A}{Character indicating the type of prior over partitions of the intercepts. Must be one of "ep" (for a Ewens-Pitman style prior) or "unif" (for the uniform prior). Default is "ep".}
  \item{prior_B}{Character indicating the type of prior over partitions of the slopes. Must be one of "ep" (for a Ewens-Pitman style prior) or "unif" (for the uniform prior). Default is "ep".}
  \item{verbose}{Logical. Default is TRUE.}
%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{PARTOPT} finds the a posteriori most likely pairs of partitions as described in Balocchi et al. (2022+). 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{final}{A list containing the estimated partitions for alpha and beta and the importance weights.}
\item{MAP}{A list containing the estimated initial partition corresponding to the maximum a posteriori (MAP) for alpha and beta (the MAP).}
\item{hyperpar}{Vector of scale hyperparameters.}
}
\references{
%% ~put references to the literature/web site here ~
Balocchi C., Deshpande S. K., George E. I., Jensen S. T. (2022+) Crime in Philadelphia: Bayesian Clustering with Particle Optimization. {\emph{arXiv:1912.00111}}.
}



%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
