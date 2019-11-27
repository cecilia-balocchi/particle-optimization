clean_crime.from_csv <- function(file_path_name, shapefile_path_name){
  # takes the raw data is csv format (downloaded from the Philadelphia Police Dept)
  # returns a dataframe where each line is a crime, with the location, the crime type, 
  # the date and time, and the neighborhood (depending on the shapefile specified)
  require(sp)
  # the names of the columns are specific for the dataset release by the Philly police
  # import the shapefiles. The code can be changed so that not all of them are used.
  imported_name <- load(shapefile_path_name)
  # depending on the neighborhood level you want (blocks, blockgroups, tracts)
  # the dataset imported have different names. In this way I can just use the new variable 
  # 'neighborhood'
  neighborhoods <- get(imported_name) 
  crime.pre <- read.csv(file_path_name,
                        sep = ",",
                        header = T,
                        stringsAsFactors = F)
  # GPS location of crime:
  temp.pointX <- crime.pre$point_x
  temp.pointY <- crime.pre$point_y
  
  ## actually they have the same NAs, but just in case
  rem.ind.points <- is.na(temp.pointX) | is.na(temp.pointY)
  ## want to make this exact same as old file. Remove crap points, keep order:
  crimes.new <- crime.pre[!rem.ind.points,c('dc_dist', 'psa', 'ucr_general',
                                            'text_general_code', 'point_x', 'point_y')]
  ## crime type oh yes - the order matters (for crimeframe)
  crimetypes <- c('murder', 'rape', 'robbery', 'assault',
                  'burglary', 'theft', 'motortheft',
                  'otherassault', 'arson', 'forgery',
                  'fraud', 'embezzlement', 'receivestolen',
                  'vandalism', 'weaponviolation',
                  'prostitution', 'sexoffense', 'drugviolation',
                  'gambling', 'familyoffense', 'dui', 'liquorlaw',
                  'publicdrunk', 'disorderly', 'vagrancy',
                  'other')
  crimeframe <- data.frame(ucr = sort(unique(crime.pre$ucr_general)),
                           crimetype = crimetypes)
  # the following reorders the crimetype for the ucr in pre.crime, so we can add crimetype to crime.new
  temp.crimetype = crimeframe$crimetype[match(crime.pre$ucr_general, crimeframe$ucr)]
  crimes.new$crimetype = temp.crimetype[!rem.ind.points]
  ## there's a timedate col, but it's in 12h... so I won't use it
  temp.date = as.Date(crime.pre$dispatch_date, tz = 'EST')
  temp.time = crime.pre$dispatch_time
  crimes.new$date = temp.date[!rem.ind.points]
  temp.timedate = as.POSIXct(paste(temp.date, temp.time), tz = 'EST')
  crimes.new$timedate = temp.timedate[!rem.ind.points]
  ## add time with any given Sunday (only so that it has a Sunday zero point
  ## for week plotting)
  temp.time <- format(crimes.new$timedate, format = "%H:%M:%S", tx = 'EST')
  given.sunday <- '2016-05-08'
  crimes.new$time <- as.POSIXct(paste(given.sunday, temp.time),  tz = "America/New_York")
  # the time zone does not work for me (?)
  ## add weekday
  crimes.new$weekday <- weekdays(crimes.new$timedate)

  ##--------- Now the spatial part ----------##
  
  ## spatial situ
  pts <- SpatialPoints(crimes.new[ ,c('point_x', 'point_y')])
  
  ## This is where we check in which region the point falls!
  # over returns a list: each element of the list corresponds to a region
  # and contains a vector with the indices of the crimes happened in that region
  # so the length of that vector is the number of crimes in that region
  crimebyngbhs <- sp::over(neighborhoods, pts, returnList = T)
  
  ## we want to record which block a crime happens in, in the dataframe
  crimes.new$ngbh <- NA
  crimes.new$ngbhID <- NA
  
  # The identifier we'll use is GEOID10 (a factor)
  crimevector <- unlist(crimebyngbhs)
  crimengbhslength <- lapply(crimebyngbhs, length)
  
  levels_region <- levels(neighborhoods@data$GEOID10)[neighborhoods@data$GEOID10]
  
  crimeind.temp <- c(0, cumsum(crimengbhslength))
  for(i in 1:length(crimebyngbhs)){
    temp.index <- (crimeind.temp[i] + 1): (crimeind.temp[i + 1])
    crimes.new$ngbh[crimevector[temp.index]] = i
    crimes.new$ngbhID[crimevector[temp.index]] = levels_region[i]
    # if(i %% 100 == 0) cat(i, " ")
  }
  
  ##------------------------------------##
  
  ## add year
  temp.year = format(crimes.new$date, '%Y')
  crimes.new$year <- as.numeric(temp.year)
  
  ## sort the whole thing
  crimes.new <- crimes.new[order(crimes.new$timedate),]
  
  ## a small fraction (863 or 0.04%) aren't in Philly (according to our shapefile)
  ## remove them
  ## (note: each crime either has both a block and blockgroup, or neither)
  crimes.new <- crimes.new[!is.na(crimes.new$ngbh),]
  return(crimes.new)
}

aggregate_crime.from_rdata <- function(crimes.new, last_year = 2017){
  # takes the dataframe with the crimes (the one formatted with clean_crime.from_csv)
  # returns a dataframe with one row for each region and year and the number of violent and nonviolent crimes for that (region, year)
  # The region used is the ngbh created in clean_crime.from_csv
  variable <- "ngbh"
  ncols <- last_year - 2006 + 1
  require(dplyr)
  violent <- crimes.new[crimes.new$ucr_general <= 400 | crimes.new$ucr_general == 800 | crimes.new$ucr_general == 1700,]
  ## if we wanted to use nonviolent crimes:
  # nonviolent <- crimes.new[crimes.new$ucr_general == 600 | crimes.new$ucr_general == 500,]
  
  v.table <- table(violent[,variable],violent$year)
  # Only up to year 'last_year' (included). Note that you should not use the current year 
  # because it contains less months of crimes.
  total <- v.table[,1:ncols]
  n <- dim(v.table)[1]
  total <- cbind(as.numeric(rownames(total)),total)
  colnames(total)[1] <- variable
  total <- as.data.frame(total)

  # save the GEOID for those regions
  varID <- paste0(variable,"ID")
  total$GEOID10 <- violent[match(total[,1], violent[,variable]),varID]
  return(total)
}


aggregate_crime.complete <- function(crimes.new, last_year = 2018){
  # takes the dataframe with the crimes (ideally the one formatted with clean_crime.from_csv)
  # returns a dataframe with one row for each region and year and the number of violent and nonviolent crimes for that (region, year)
  # The region used depends on 'variable': it can be 'block', 'blockgroup' or 'tract'.
  variable <- "ngbh"
  ncols <- last_year - 2006 + 1
  require(dplyr)
  violent <- crimes.new[crimes.new$ucr_general <= 400 | crimes.new$ucr_general == 800 | crimes.new$ucr_general == 1700,]
  nonviolent <- crimes.new[crimes.new$ucr_general == 600 | crimes.new$ucr_general == 500,]
  
  v.table <- table(violent[,variable],violent$year)
  nv.table <- table(nonviolent[,variable],nonviolent$year)
  
  colnames(v.table) <- paste0("v_",colnames(v.table))
  colnames(nv.table) <- paste0("nv_",colnames(nv.table))
  # Only up to year 'last_year' (included). Note that you should not use the current year 
  # because it contains less months of crimes.
  total <- cbind(v.table[,1:ncols],nv.table[,1:ncols])
  n <- dim(v.table)[1]
  total <- cbind(as.numeric(rownames(total)),total)
  colnames(total)[1] <- variable
  total <- as.data.frame(total)
  
  # save the GEOID for those regions
  varID <- paste0(variable,"ID")
  total$GEOID10 <- violent[match(total[,1], violent[,variable]),varID]
  
  c_y <-tbl_df(total)
  t1 <- c_y[rep(1:dim(c_y)[1],ncols),]
  # sort by variable
  t2 <- arrange_(t1, .dots = variable)
  # add 'year' column
  t3 <- mutate(t2,year = rep(2006:(ncols+2005),dim(c_y)[1]))
  # violent will contain the number of violent crimes for that year and region
  t4 <- mutate(t3,violent = NA)
  for(y in 2006:(ncols+2005)){
    str <- paste0("v_",y)
    t4[t4$year==y,]$violent <- c_y[,str][[1]]
  }
  # same for nonviolent
  t5 <- mutate(t4, nonviolent = NA)
  for(y in 2006:(ncols+2005)){
    str <- paste0("nv_",y)
    t5[t5$year==y,]$nonviolent <- c_y[,str][[1]]
  }
  years <- "year"
  violents <- "violent"
  nonviolents <- "nonviolent"
  geo <- "GEOID10"
  # only keep the relevant variables:
  t5 <- dplyr::select_(t5,.dots = c("variable","geo","years","violents","nonviolents"))
  t5[,variable][[1]] <- as.factor(t5[,variable][[1]])
  t5[,geo][[1]] <- as.factor(t5[,geo][[1]])
  t6 <- mutate_(t5, .dots = setNames(list(variable), "X"))
  t6 <- dplyr::select(t6,X,GEOID10,year,violent,nonviolent)
  d6 <- as.data.frame(t6)
  return(d6)
}
