# 02_createContinentHexRates.R
# Author: Benjamin R. Goldstein
# Date: 6/3/2021

# This file executes code for turning iNaturalist and eBird raw inputs into 
# counts for many species, aggregated to a hexagonal grid. 
# It calls the helper file contHexDetectionRate_fn.R which contains much
# of the code used for achieving the above.

require(rgdal)
require(sp)
require(maptools)
require(raster)
require(tidyverse)
require(taxalight)

# Read in the list of state codes
statecode_df <- read_csv("raw_data/statecodes.csv")

source("helper_code/contHexDetectionRate_fn.R")

# Create the continent-wide hex grid 
continentcode <- statecode_df %>% filter(!Abbreviation %in% c("AK", "HI"))
if (!file.exists("intermediate_data/continent_hexes/continent_hexes.shp")) {
  US_map <- maps::map(database = "usa", regions = "main",
                      fill = T, plot = F)
  IDs <- sapply(strsplit(US_map$names, ":"), function(x) x[1])
  US_poly <- map2SpatialPolygons(US_map, IDs = IDs, 
                                 proj4string=CRS("+proj=longlat +datum=WGS84"))
  US_poly <- spTransform(US_poly, CRS("+init=epsg:2163 +datum=WGS84"))
  
  US_grid <- make_grid(US_poly, cell_diameter = 40000)
  US_grid <- SpatialPolygonsDataFrame(
    US_grid, data = data.frame(gridID = 1:length(US_grid))
  )
  
  statesmap <- maps::map(database = "state", fill = T, plot = F)
  IDs <- sapply(strsplit(statesmap$names, ":"), function(x) x[1])
  statesmap <- map2SpatialPolygons(statesmap, IDs = IDs,
                                   CRS("+proj=longlat +datum=WGS84"))
  statesmap <- spTransform(statesmap, CRS("+init=epsg:2163 +datum=WGS84"))
  
  US_grid$state <- names(statesmap)[US_grid %over% statesmap]
  
  writeOGR(US_grid, dsn = "continent_hexes", 
           layer = "continent_hexes", 
           driver = "ESRI Shapefile", 
           overwrite_layer = TRUE)
} else {
  US_grid <- readOGR("intermediate_data/continent_hexes/continent_hexes.shp", "continent_hexes")
}



# Batch the creation of hex-level files which will eventually be summarized 
# into a single dataset.
# It's critical to wipe before batching so items don't get duplicated!
system("rm intermediate_data/continentHexRates/hex*40k.csv")
continentcode %>% 
  mutate(done = FALSE) %>% 
  write_csv("intermediate_data/continentHexRates/stateProcInfo.csv")
for (i in 1:nrow(continentcode)) {
  # This version of the fn reads to and writes the files internally. It writes
  # many small files to store all this data (1/hex cell)
  writeLines(continentcode$State[i])
  msg <- suppressMessages(suppressWarnings(
    contHexDetectionRate(statename = continentcode$State[i],
                         gridObj = US_grid, 
                         outputPath = "intermediate_data/continentHexRates",
                         outputSuffix = "_40k")
  ))
}



# create summary
hex_files <- list.files("intermediate_data/continentHexRates", pattern = "hex.*40k", full.names = T)

suppressMessages(all_hex_info <- lapply(hex_files, read_csv))
all_hex_info <- do.call(rbind, all_hex_info)

hex_totals <- all_hex_info %>% 
  group_by(hex) %>% 
  summarize(inat_ntotal = sum(inat_n),
            ebird_comp_ntotal = sum(ebird_comp_n))

all_hex_info <- left_join(all_hex_info, hex_totals)

write_csv(all_hex_info, "intermediate_data/continentHexRates/continentHexInfo.csv")

system("rm intermediate_data/continentHexRates/hex*40k.csv")

# Filter for species that were observed at least once in all 3 datasets,
# to ensure that naming conventions are the same
spec_totals <- all_hex_info %>% 
  group_by(species) %>% 
  summarize(inat_ntotal = sum(inat_n),
            ebird_comp_ntotal = sum(ebird_comp_n))

accepted_specs <- spec_totals %>% 
  filter(inat_ntotal > 0 &
         ebird_comp_ntotal > 0)



# Get a list of eBird species and their common names
all_ebird_species <- list()
for (i in 1:length(statecode_df$Abbreviation)) {
  code <- statecode_df$Abbreviation[i]
  all_ebird_species[[i]] <- read_csv(paste0("../eBird_Data/subsets-2020/state_rdbs/", 
                               code, "_Aug2020_species_counts.csv")) %>% 
    distinct(name_clean, SCIENTIFIC.NAME)
}

# Write output
all_ebird_species_df <- do.call(rbind, all_ebird_species) %>% distinct() %>% 
  write_csv("raw_data/ebird_species_tbl.csv")
