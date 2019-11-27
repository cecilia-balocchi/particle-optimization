# Note the crime dataset (CSV) downloaded from the Philly Police Dept 
# (https://www.opendataphilly.org/dataset/crime-incidents)
# should be saved in `data/csv/` as `incidents_part1_part2.csv`
# The sub-diectory `data/shapefile` contains the geographic files 
# for aggregating over the neighborhoods (these shapefiles where downloaded from 
# opendataphilly.org, read into R using 'readOGR' and saved as Rdata objects)

source("get_data/subsetup/clean_crime.R")
# crimes.new is a dataset where each row represents a crime and we save time/location and crimetype info
crimes.new <- clean_crime.from_csv(file_path_name = "data/csv/incidents_part1_part2.csv", 
                                   shapefile_path_name = "data/shapefile/phillytracts")
# crime_agr is the dataframe I use to run the linear models: each line gives the number of crimes for a region and a year.
# use it with variable = "blockgroup" or "block" or "tract" (it will datasets aggregated at different levels).
crime_agr <- aggregate_crime.from_rdata(crimes.new, last_year = 2017)
only_count <- crime_agr[,2:(ncol(crime_agr)-1)]
# let's transform crime with inverse hyperbolic sine transformations
# log(yi+(yi2+1)1/2)
# http://worthwhile.typepad.com/worthwhile_canadian_initi/2011/07/a-rant-on-inverse-hyperbolic-sine-transformations.html

# Note: the function we said we would consider in the paper has a normalizing factor:
# log(y+(y^2+1)^0.5) - log(2)
# Instead of doing the analysis all over again, we can simply consider alpha_i - log(2)
ihst <- function(y){
  log(y+(y^2+1)^0.5)
}
transf_count <- ihst(only_count)
transf_count <- as.matrix(transf_count)
write.table(x = transf_count, file = "data/Y_tracts.csv", col.names = F,
            row.names = F, sep = ',')
X <- matrix(rep(1:ncol(transf_count), each = nrow(transf_count)), nrow = nrow(transf_count))
write.table(x = X, file = "data/X_tracts.csv", col.names = F,
            row.names = F, sep = ',')

# Let's now save the adjacency matrix
require(spdep)
load("data/shapefile/phillytracts")
# get the adjacency matrix from the shapefile
n_tr <- dim(tracts)[1]
list.poly <- poly2nb(tracts)
w <- matrix(data = 0, nrow = n_tr, ncol = n_tr)
for(i in 1:n_tr){
  w[i, list.poly[[i]]] <- 1
}
# for some arcane reason it's not symmetric
w.sym <- w + t(w) - w*t(w)
write.table(x = w.sym, file = "data/A_tracts.csv", col.names = F,
            row.names = F, sep = ',')

# this is useful for out of sample predictions.
# it is not the most efficient way to do it..
crime_tract2019 <- aggregate_crime.complete(crimes.new, last_year = 2018)
crime_tract2019$tr.violent <- ihst(crime_tract2019$violent)
crime_tract2019$log.violent <- log(crime_tract2019$violent+1)
crime_tract2019$year19 <- crime_tract2019$year-2005
save(crime_tract2019, file = "data/crime_tract2019.rdata")
