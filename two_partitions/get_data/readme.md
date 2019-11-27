## Aggregate crime data from the csv downloaded for the Philadelphia Police Department

In `setup_crime.R`, I create the aggregated datasets from the csv crime file, both aggregating over block groups and over tracts. 

Instructions: 
1. Download the crime CSV dataset from https://www.opendataphilly.org/dataset/crime-incidents (download the CSV format) and save the file in `data/csv` (the filename should be `incidents_part1_part2.csv`).
2. Download the shapefile for Philadelphia's census tracts from [open data Philly](https://www.opendataphilly.org/dataset/census-tracts) and save the compressed folder in `data/shapefiles` and unzip it there, so that the shapefiles are in `data/shapefiles/tracts_pc2010/` (for example). Note that the folder will contain a .shp file and other files with the same name but different extensions. Do not delete or change the name to those files.
3. Import the shapefiles in R and save as Rdata file using `get_data/import_shapefiles.R`.
3. Run the script `get_data/setup_crime.R`. This will create `Y_tracts.csv`, `X_tracts.csv` and `A_tracts.csv`. Moreover it will create `crime_tract2019.rdata`, which can be used for out-of-sample predictions.

I use the transformation Inverse Hyperbolic Sine transformation. 
I include years from 2006 to 2017. 
The functions to make all of these changes are contained in `subsetup/clean_crime.R`.

This code can word with different neighborhood levels. You can use other shapefiles, as the ones for census blocks or census block groups.
The shapefiles for Phiadelphia can be downloaded at
https://www.opendataphilly.org/dataset/census-block-groups,
https://www.opendataphilly.org/dataset/census-blocks and https://www.opendataphilly.org/dataset/census-tracts 
(NOTE: I used the shapefiles from 2010).
For the analysis in the paper, census tracts only where used.

Some of this code was created with the help of [Colman Humphrey](https://github.com/ColmanHumphrey/). The scripts here have been adapted from another of my reposotories [Urban-project](https://github.com/cecilia-balocchi/Urban-project).

#### Get shapefiles and import them in R

In 'import_shapefiles.R' we use `readOGR` (from the package 'rgdal') to read and import a shapefile into R. Remember that the folder containing your '.shp' file should also contain a bunch of other files with the same name but different extension (.shx, .prj, .dbf, .cpg).

#### Get google map

Unfortunately I forgot how I created an account with the google map API without having to pay. Once I got an API key I downloaded the googlemap using the R package 'ggmap' and saved it as an .rdata object. The code to download the map using ggmap is in `save googlemap.R`.
