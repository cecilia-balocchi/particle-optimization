library(rgdal)
# adjust the directory and file names accordingly
path_to_downloaded_shapefile <- "data/shapefile/tracts_pc2010/"
path_to_export_shapefile <- "data/shapefile/"
file_name <- "tl_2010_42101_tract10.shp"
new_name <- "phillytracts"
str <- paste0(path_to_downloaded_shapefile, file_name)
if(file.exists(str)){
	obj <- readOGR(str)
	save(obj, file = paste0(path_to_export_shapefile,new_name))
}
