## Aggregate crime data from the csv downloaded for the Philadelphia Police Department

`create_crime_tracts2019.R` creates the aggregated datasets from the downloaded csv crime file, by aggregating over census tracts and then transforming to obtain crime density. It outputs the object `data/crime_tracts2019.rdata`. Detailed instructions are described below.

For convenience I have already saved `crime_tracts2019.rdata` and the Philadelphia shape file in the repository, respectively in `data/crime_tracts2019.rdata` and `data/shapefile/phillytracts`. To create these files, the data was accessed in February 2019.

`create_tracts_data.R` uses `data/crime_tracts2019.rdata` to output the csv files used for running the analysis (`data/Y_density.csv`, `data/X_density.csv`, `data/E_density.csv` and `data/A_block_density.csv`).

#### Instructions to create shapefiles and crime_tracts2019 from raw data: 
1. Download the crime CSV datasets from https://www.opendataphilly.org/dataset/crime-incidents for each year from 2006 to 2018. Combine the CSV files for all the years into one CSV file (make sure you don't copy the header multiple times) named `incidents_part1_part2.csv` and save it in `data/csv`. Note: I was able to download all the incidents from 2006 to 2018 (included) as one csv file from this link: https://phl.carto.com/api/v2/sql?filename=incidents_part1_part2&format=csv&skipfields=cartodb_id,the_geom,the_geom_webmercator&q=SELECT%20*%20,%20ST_Y(the_geom)%20AS%20lat,%20ST_X(the_geom)%20AS%20lng%20FROM%20incidents_part1_part2%20WHERE%20dispatch_date_time%20%3E=%20%272006-01-01%27%20AND%20dispatch_date_time%20%3C%20%272019-01-01%27.
2. Download the shapefile for Philadelphia's census tracts from [open data Philly](https://www.opendataphilly.org/dataset/census-tracts) and save the compressed folder in `data/shapefiles` and unzip it there, so that the shapefiles are in `data/shapefiles/tracts_pc2010/` (for example). Note that the folder will contain a .shp file and other files with the same name but different extensions. Do not delete or change the name to those files.
3. Import the shapefiles in R and save as Rdata file using `get_data/import_shapefiles.R`.
3. Run the script `get_data/setup_crime.R`. This will create `Y_density.csv`, `X_density.csv` and `A_block_density.csv`. Moreover it will create `crime_tract2019.rdata`, which can be used for out-of-sample predictions.

I use the transformation Inverse Hyperbolic Sine transformation. 
I include years from 2006 to 2017. 
The functions to make all of these changes are contained in `subsetup/clean_crime.R`.

This code can word with different neighborhood levels. You could use other shapefiles, as the ones for census blocks or census block groups.
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
