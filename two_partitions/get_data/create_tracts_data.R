load("data/crime_tract2019.rdata")
load("data/shapefile/phillytracts")
ncols <- 13

# we save the csv in matrix form, with each row corresponding to a neighborhood, each column to a year
Y <- matrix(crime_tract2019$tr.density, ncol = ncols, byrow = T)
X <- matrix(crime_tract2019$year19, ncol = ncols, byrow = T)
E <- matrix(crime_tract2019$ALAND10_sqmile, ncol = ncols, byrow = T)

# we have formatted these variables to include also the year 2018 as the last column.
# however, since we'll use that year for out of sample prediction, we remove it now
Y <- Y[,-ncols]
X <- X[,-ncols]
E <- E[,-ncols]

# Let's now save the adjacency matrix
require(spdep)
# get the adjacency matrix from the shapefile
n_tr <- dim(tracts)[1]
list.poly <- poly2nb(tracts)
w <- matrix(data = 0, nrow = n_tr, ncol = n_tr)
for(i in 1:n_tr){
  w[i, list.poly[[i]]] <- 1
}
# for some arcane reason it's not symmetric
w.sym <- w + t(w) - w*t(w)

# let's now save as CSV in data
write.table(w.sym, row.names = FALSE, col.names = FALSE, sep = ',', file = "data/A_block_density.csv")
write.table(X, row.names = FALSE, col.names = FALSE, sep = ',', file = "data/X_density.csv")
write.table(Y, row.names = FALSE, col.names = FALSE, sep = ',', file = "data/Y_density.csv")
write.table(E, row.names = FALSE, col.names = FALSE, sep = ',', file = "data/E_density.csv")

