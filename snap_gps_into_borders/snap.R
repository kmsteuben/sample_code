## ############################################################################################################
##
## GEOGRAPHY CODEBOOKS SNAP FUNCTION
## Author: Krista Steuben (steuben)
## Purpose: To snap points that are less than 10km outside of the border of the assigned country
##          into the boundary of the country. (note: 10km is the LBD standard)
##
## Sub Functions:
##   set_up
##   format_codebook
##   snap_calculate_points
##   codebook_diagnostics
##   snap_diagnostics
##
## Principle Functions:
##   snap_codebook
##   snap_all_codebooks
##
## ############################################################################################################


## ############################################################################################################
## SUB FUNCTIONS
## ############################################################################################################

# Installs libraries needed for snap function
set_up <- function(){
  rm(list=ls())
  libs <- c("geosphere", "rgdal", "plyr")
  sapply(libs, require, character.only = T)
}

# Prepares the shapefile by creating an ISO3 code that matches codebooks
# NEEDS TO BE UPDATED AS SHAPEFILES ARE UPDATED
prepare_gaul_lookup <- function(admin0){
  temp <- read.csv('<<<< FILE REDACTED >>>>')
  keep <- temp[,c("location_name", "iso3")]
  keep$iso3 <- as.character(keep$iso3)
  gaul_lookup <- merge(admin0@data, keep, by.x = "ADM0_NAME", by.y = "location_name", all.x = T)
  gaul_lookup$iso3 <- as.character(gaul_lookup$iso3)
  gaul_lookup$iso3[is.na(admin0@data$iso3)] <- "TEMP"
  gaul_lookup$iso3[gaul_lookup$ADM0_CODE == 45] <- "CIV"
  gaul_lookup$ADM0_CODE <- as.character(gaul_lookup$ADM0_CODE)
  return(gaul_lookup)
}

# Subsets the codebook into only point data and formats lat/long into numeric
format_codebook <- function(codebook, first_snap){
  point.locations <- subset(codebook, point == 1)
  point.locations <- subset(point.locations, !is.na(long) | !is.na(lat))
  point.locations$iso3 <- as.character(point.locations$iso3)
  # Fills in blank original lat/longs with lat/long columns info
  point.locations$original_lat[is.na(point.locations$original_lat)] <- as.numeric(as.character(point.locations$lat[is.na(point.locations$original_lat)]))
  point.locations$original_long[is.na(point.locations$original_long)] <- as.numeric(as.character(point.locations$long[is.na(point.locations$original_long)]))
  if (first_snap == T) {
    point.locations$lat <- point.locations$original_lat # only used for first snap of codebook/shapefile
    point.locations$long <- point.locations$original_long # only used for first snap of codebook/shapefile
  }
  if (exists('point.locations') && nrow(point.locations) > 0) point.locations$snapped <- NA # only used for first snap of codebook
  return(point.locations)
}

# Finds all points outside borders of assigned countries within a given codebook
# Returns a data.frame of out_of_bounds points, their distances to the country border, and the nearest point on border
snap_calculate_points <- function(codebook, admin0){
  # List countries in codebook
  countries <- sort(unique(codebook$iso3))
  message("Processing ", length(countries), " countries.")
  out_of_bounds <- data.frame()
  # Loop through countries and find points outside of national boundaries
  for (i in 1:length(countries)){
    country <- countries[1]

    #subset and prep survey to country
    sub <- unique(codebook[codebook$iso3 == country, ])

    #subset shapefile to survey country
    gaul_code <- na.exclude(gaul_lookup$GAUL_CODE[gaul_lookup$iso3 == sub$iso3[1]])
    sub.shapefile <- admin0[admin0@data$GAUL_CODE == gaul_code, ]

    # get points outside boundaries
    sub_CRS <-  CRS(proj4string(sub.shapefile))
    outside_pts <- sub[is.na(over(SpatialPoints(sub[, c("long", "lat")], sub_CRS),
                                                as(sub.shapefile, "SpatialPolygons"))),]


    # calculate the distance between point and polygon and note the nearest point on the polygon
    if (nrow(outside_pts) != 0) {
      sp <- SpatialPoints(outside_pts[,c("long", "lat")], sub_CRS)
      nearest_point <- dist2Line(sp, as(sub.shapefile, "SpatialPolygons"), distfun = distHaversine)
      outside_pts$min_distance_km <- nearest_point[,1]/1000
      outside_pts$target_lat <- nearest_point[,3]
      outside_pts$target_long <- nearest_point[,2]

      out_of_bounds <- rbind(out_of_bounds, outside_pts)
    }
  }
  return(out_of_bounds)
}

# Prints the number of rows in a given codebook that have invalid lat/long, invalid ISO3, and pionts >10km from border
codebook_diagnostics <- function(bad_points, bad_iso3, ones_gps, out_of_bounds_over_10){
  message("        Running Codebook Diagnostics")
  if (exists('ones_gps') && nrow(ones_gps) > 0){
    message("         WARNING! THERE ARE (1,1) IN THE LAT/LONG COLUMNS! FIX AND RE-RUN SNAP")
  }
  message("            There are ", nrow(bad_points), " points with invalid lat/longs.")
  message("            There are ", nrow(bad_iso3), " points with an invalid ISO3.")
  message("            There are ", nrow(out_of_bounds_over_10), " points >10km out of country border.")

}

# Prints diagnostics of snap on codebook. Calculates the distance from original lat/long and snapped point.
# Prints the number of points snapped and number of points snapped more than 12km
snap_diagnostics <- function(point.locations){
  message("        Running Snap Diagnostics")
  snapped.points <- subset(point.locations, snapped == "Snapped")
  n <- nrow(snapped.points)
  if (n >0){
    for (i in 1:n){
      snapped.points$dist_from_original[i] <- distm(c(snapped.points$original_long[i],snapped.points$original_lat[i]),
                                                    c(snapped.points$long[i],snapped.points$lat[i]),
                                                    fun = distHaversine)
    }
  }
  snapped.points$dist_from_original <- snapped.points$dist_from_original/1000
  too_far <- subset(snapped.points, dist_from_original > 12) #Using 12km because point should be <10km from border and be snapped 1km into border
  message("            There were ", nrow(snapped.points), " points snapped.")
  message("            There are ", nrow(too_far), " points >12 km from starting point.")
}

## ############################################################################################################
## PRINCIPLE FUNCTIONS
## ############################################################################################################

#Snaps one given codebook
snap_codebook <- function(codebook_filepath, admin0, gaul_lookup, first_snap) {
  message("Processing codebook: ", codebook_filepath)
  codebook <- read.csv(codebook_filepath)
  codebook$ID <- seq.int(nrow(codebook))
  point.locations <- format_codebook(codebook, first_snap)
  not_point.locations <- codebook[!(codebook$ID %in% point.locations$ID),]
  bad_points <- subset(point.locations, lat >90 | lat < -90 | long >180 | long < -180)
  bad_iso3 <- subset(point.locations, !(iso3 %in% gaul_lookup$iso3))
  ones_gps <- subset(point.locations, lat == 1 | long == 1)
  #previously_snapped <- subset(point.locations, snapped %in% c("Not snapped - over 10km away"))
  point.locations <- subset(point.locations,
                            !(lat >90 | lat < -90 | long >180 | long < -180)
                            & (iso3 %in% gaul_lookup$iso3)
                            & !(lat == 1 | long == 1))
  to_snap <- point.locations[,c('ID','iso3', 'lat','long')]
  out_of_bounds_over_10 <- data.frame()

  #repeat snap until no points are out_of_bounds (up to 10 iterations - 10 iterations were sufficient in testing)
  iteration <- 1
  while (nrow(to_snap)>0 & iteration <= 11){
    if (nrow(to_snap) > 0){
      message("Iteration ", iteration)
      out_of_bounds_snapped <- snap_calculate_points(to_snap, admin0)
      message(nrow(out_of_bounds_snapped), " points out of bounds.")
      if (nrow(out_of_bounds_snapped) > 0) {
        if (iteration == 1){
          out_of_bounds_over_10 <- subset(out_of_bounds_snapped, min_distance_km > 10)
        }
        out_of_bounds_snapped <- subset(out_of_bounds_snapped, min_distance_km <= 10)
        if (nrow(out_of_bounds_snapped) > 0) {
          # new lat/long calculated by finding the point 1km inside country boundary along line between original lat/long and target lat/long
          out_of_bounds_snapped$new_lat <- (1/out_of_bounds_snapped$min_distance_km*(out_of_bounds_snapped$target_lat - out_of_bounds_snapped$lat)) + out_of_bounds_snapped$target_lat
          out_of_bounds_snapped$new_long <- (1/out_of_bounds_snapped$min_distance_km*(out_of_bounds_snapped$target_long - out_of_bounds_snapped$long)) + out_of_bounds_snapped$target_long
          # the calculated new lat/longs replace the lat/long in point.locations
          point.locations$lat[match(out_of_bounds_snapped$ID, point.locations$ID)] <- out_of_bounds_snapped$new_lat
          point.locations$long[match(out_of_bounds_snapped$ID, point.locations$ID)] <- out_of_bounds_snapped$new_long
          point.locations$snapped[match(out_of_bounds_snapped$ID, point.locations$ID)] <- "Snapped"
          out_of_bounds_snapped$lat <- out_of_bounds_snapped$new_lat
          out_of_bounds_snapped$long <- out_of_bounds_snapped$new_long
        }
      }
    to_snap <- out_of_bounds_snapped
    iteration <- iteration + 1
    }
  }

  point.locations$snapped <- as.character(point.locations$snapped)
  if (exists('to_snap') && nrow(to_snap) > 0){
    message(paste0("There are still ", nrow(to_snap),"  out of bounds points!"))
    point.locations$snapped[match(to_snap$ID, point.locations$ID)] <- "Flag - still out of bounds"
  }

  if (exists('bad_points') && nrow(bad_points) > 0){
    bad_points$snapped <- "Flag - bad points"
  }
  if (exists('bad_iso3') && nrow(bad_iso3) > 0){
    bad_iso3$snapped <- "Flag - bad Iso3"
  }
  if (exists('ones_gps') && nrow(ones_gps) > 0){
    ones_gps$snapped <- "Flag - ones in lat/long"
  }
  if (exists('out_of_bounds_over_10') && nrow(out_of_bounds_over_10) > 0){
    point.locations$snapped[point.locations$ID %in% out_of_bounds_over_10$ID] <- "Not snapped - over 10km away"
  }

  codebook_diagnostics(bad_points, bad_iso3, ones_gps, out_of_bounds_over_10)
  snap_diagnostics(point.locations)
  final_codebook <- rbind.fill(list(point.locations,not_point.locations, bad_points, bad_iso3, ones_gps))
  final_codebook <- final_codebook[order(final_codebook$ID),]
  final_codebook$ID <- NULL
  return (final_codebook)
}

# Runs snap function on all codebooks in geography codebook directory
snap_all_codebooks <- function(first_snap = FALSE){
  set_up()
  message("Opening shapefile ...")
  admin0 <- readOGR("<<<< FILE REDACTED >>>>")
  gaul_lookup <- prepare_gaul_lookup(admin0)
  all_codebooks <- list.files(paste0("<<<< FILE REDACTED >>>>"), pattern=".csv$", ignore.case = T, full.names = T)
  for (codebook_filepath in all_codebooks){
    final_codebook <- snap_codebook(codebook_filepath, admin0, gaul_lookup, first_snap)
    codebook_name <- grep("/*.csv", codebook_filepath)
    write.csv(final_codebook, paste0("<<<< FILE REDACTED >>>>",basename(codebook_filepath)), na = "", row.names = F)
  }
}
