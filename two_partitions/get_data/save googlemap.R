require(ggmap)
# put google maps API key here (check console.cloud.google.com for the key)
# str <- 
temp.center <- c(-75.12, 39.99)
googlemap <- get_googlemap(center = temp.center,
                           zoom = 11, ## this is important
                           scale = 2,
                           color = "bw",
                           maptype = "roadmap", key = str)
save(googlemap, file = "data/phillygooglemap.rdata")
