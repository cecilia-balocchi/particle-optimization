### This file creates the object "crime_tract2019.rdata" from scratch
### For convenience it's been already saved in the folder data/
### So this script can be skipped.

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
# check you have all the years
all.equal(sort(unique(crimes.new$year)),2006:2018)
# crime_agr is the dataframe I use to run the linear models: each line gives the number of crimes for a region and a year.
# use it with variable = "blockgroup" or "block" or "tract" (it will datasets aggregated at different levels).
last_y <- 2018
crime_agr <- aggregate_crime.from_rdata(crimes.new, last_year = last_y)

### only count is in matrix form, with each column representing a year 
### and each row representing a neighborhood. let's transform it so that
### each row is a pair (neighborhood, year) with the number of crimes.
ncols <- last_y - 2006 + 1
c_y <- crime_agr
c_y <- as.data.frame(c_y)
t1 <- c_y[rep(1:dim(c_y)[1],ncols),]
t2 <- arrange(t1, ngbh)
t3 <- mutate(t2,year = rep(2006:(ncols+2005),dim(c_y)[1]))
t4 <- mutate(t3,violent = NA)
for(y in 2006:(ncols+2005)){
  t4[t4$year==y,"violent"] <- c_y[,as.character(y)][[1]]
}
t5 <- dplyr::select(t4, ngbh, GEOID10, year, violent)
t5[,"ngbh"] <- as.factor(t5[,"ngbh"])
t5[,"GEOID10"] <- as.factor(t5[,"GEOID10"])
t6 <- rename(t5, X = ngbh)
crime_tract2019 <- as.data.frame(t6)
crime_tract2019$year19 <- crime_tract2019$year - 2006 + 1
## the name is a bit confusing, because it actually contains only crimes up to 2018 
## the reason for the confusing name is that it was computed in 2019 with the full-year data up to 2018


# import the shapefile for area measurement of each neighborhood
load("data/shapefile/phillytracts")
violent <- crime_tract2019$violent 
tmp <- tracts@data[,c("GEOID10","ALAND10")] 
index <- match(crime_tract2019$GEOID10, tmp$GEOID10)
area_sqmet <- tmp$ALAND10[index] # it's in square meters!
area_sqkm <- area_sqmet/(1000*1000)
area_sqmile <- area_sqkm/2.59
density <- violent/area_sqmile

crime_tract2019$ALAND10 <- area_sqmet
crime_tract2019$ALAND10_sqmile <- area_sqmile
crime_tract2019$density <- density

ihst_l2 <- function(y){
  log(y+(y^2+1)^0.5)-log(2)
}
tr.density <- ihst_l2(crime_tract2019$density)
crime_tract2019$tr.density <- tr.density

save(crime_tract2019, file = "data/crime_tract2019.rdata")
