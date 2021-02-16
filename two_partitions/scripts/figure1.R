source("scripts/plot_functions.R")

# Make MLE plots
str <- "_density.csv"
Y <- read.csv(paste0("data/Y",str), header = F)
X <- read.csv(paste0("data/X",str), header = F)
Y <- as.matrix(Y)
X <- as.matrix(X)
X <- X - rowMeans(X)
X <- X/apply(X, 1, sd)
n_tr <- nrow(Y)

betas.mle <- numeric(n_tr)
for(i in 1:n_tr)
  betas.mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
alphas.mle <- rowMeans(Y)


phillypoly <- tracts
polyfortified <- fortify(phillypoly)
p1 <- plot_borders_ggplot(rep(0, n_tr), phillypoly, w.sym, var = alphas.mle, polyfortified = polyfortified, legend = TRUE, 
                    legend_name = "Mean\nlevel", palette_num = 3, map = googlemap)
p2 <- plot_borders_ggplot(rep(0, n_tr), phillypoly, w.sym, var = betas.mle, polyfortified = polyfortified, legend = TRUE, 
                    legend_name = "Time\ntrend", palette_num = 5, map = googlemap)
ggsave(filename = "alpha_density_mle.png", plot = p1, device = "png", path = "figures/",
         width = 2, height = 2, units = "cm", scale = 10)
ggsave(filename = "beta_density_mle.png", plot = p2, device = "png", path = "figures/",
         width = 2, height = 2, units = "cm", scale = 10)