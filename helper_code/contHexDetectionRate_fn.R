# contHexDetectionRate_fn.R
# Author: Benjamin R. Goldstein
# Date: 6/3/2021

# This file contains helper code for turning iNaturalist and eBird raw inputs into 
# counts for many species, aggregated to a hexagonal grid. 
# It is called by the main code file 02_createContinentHexRates.


# A manual function for associating species using taxalight, with some hardocoded
# fixes for species that aren't automatically matched for some reason
manual_get_ids <- function(spec_vec, db) {
  result_vec <- character(length(spec_vec))
  duplicate_specs <- c()
  problem_specs <- c()
  
  manual_fixes <- 
    data.frame(source = c(
      "Ammospiza nelsoni", "Dryobates arizonae", "Dryobates albolarvatus",
      "Sporaeginthus subflavus", "Rhabdotorrhinus corrugatus", "Melanitta stejnegeri",
      "Ammospiza caudacuta", "Nesophlox evelynae", "Phonipara canora",
      "Spatula hottentota", "Estrilda coerulescens", "Oressochen jubatus", 
      "Leucogeranus leucogeranus", "Otocichla mupinensis", "Alopecoenas xanthonurus"
      
      
    ), correct = c(
      "Ammodramus nelsoni", "Leuconotopicus arizonae", "Leuconotopicus albolarvatus",
      "Amandava subflava", "Aceros corrugatus", "Melanitta fusca",
      "Ammodramus caudacutus", "Calliphlox evelynae", "Tiaris canorus",
      "Anas hottentota", "Estrilda caerulescens", "Neochen jubata",
      "Grus leucogeranus", "Turdus mupinensis", "Gallicolumba xanthonura"
    ))
  
  for (i in 1:length(spec_vec)) {
    spec_tbl <- taxalight::tl(spec_vec[i], provider = db)
    if (nrow(spec_tbl) == 0) {
      if (spec_vec[i] %in% manual_fixes$source) {
        spec_tbl <- 
          taxalight::tl(manual_fixes$correct[manual_fixes$source == spec_vec[i]], 
                        provider = db)
        
      }
    }
    
    if (nrow(spec_tbl) == 0) {
      problem_specs <- c(problem_specs, spec_vec[i])
      result_vec[i] <- NA
    } else if (nrow(spec_tbl) > 1) {
      duplicate_specs <- c(duplicate_specs, spec_vec[i])
      
      if (length(unique(spec_tbl$acceptedNameUsageID)) > 1) {
        spec_tbl <- spec_tbl %>% filter(scientificName == spec_vec[i])
      }
      
      result_vec[i] <- spec_tbl$acceptedNameUsageID[1]
    } else {
      result_vec[i] <- spec_tbl$acceptedNameUsageID[1]
    }
  }
  
  if (length(problem_specs) > 0) {
    warning("Problem species detected that require manual handling.")
  }
  result_vec
}


# A function for generating a hexagonal grid to cover a given polygon x
# based on either a cell diameter or cell area
make_grid <- function(x, cell_diameter, cell_area, clip = FALSE) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate array of hexagon centers
  g <- spsample(ext, type = "hexagonal", cellsize = cell_diameter, 
                offset = c(0.5, 0.5))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  return(g)
}


# For a given state name, get the iNaturalist and eBird observations for
# that state and process them onto a hexagonal grid
contHexDetectionRate <- function(statename = NULL, statecode = NULL,
                                 gridObj, outputPath, outputSuffix) {
  # Step 1: Prepare data.
  if (is.null(statename) && is.null(statecode)) {
    stop("You must provide either a state name or a 2-letter state code.")
  } else if (is.null(statename)) {
    if (!(statecode %in% statecode_df$Abbreviation)) {
      stop(paste0("State code ", statecode, 
                  " not recognized. Did you provide 2 uppercase letters?"))
    }
    statename <- statecode_df$State[statecode_df$Abbreviation == statecode]
  } else if (is.null(statecode)) {
    if (!(statename %in% statecode_df$State)) {
      stop(paste0("State name ", statename, 
                  " not recognized."))
    }
    statecode <- statecode_df$Abbreviation[statecode_df$State == statename]
  } else {
    if (which(statecode_df$State == statename) !=
        which(statecode_df$Abbreviation == statecode)) {
      stop("You provided a state code and state name that didn't match.")
    }
  }
  
  stateProcInfo <- read_csv("intermediate_data/continentHexRates/stateProcInfo.csv")
  # Can't handle state-by-state overwriting due to summing.
  # Return if this state has already been processed since reset.
  if (stateProcInfo$done[stateProcInfo$State == statename]) {
    return("Already processed.")
  }
  
  
  
  # Read data
  ebird_dat <- read_csv(paste0("raw_data/state_rdbs/", 
                               statecode, "_Aug2020_checklist_info.csv")) %>% 
    dplyr::select(SAMPLING.EVENT.IDENTIFIER,
                  ALL.SPECIES.REPORTED, OBSERVATION.DATE,
                  LATITUDE, LONGITUDE) %>% 
    right_join(read_csv(paste0("raw_data/state_rdbs/", 
                               statecode, "_Aug2020_species_counts.csv"))) %>% 
    filter(lubridate::year(OBSERVATION.DATE) < 2020) %>% 
    dplyr::select(-X1, -OBSERVATION.DATE)
  
  inat_dat <- read_csv(paste0("raw_data/iNat_", 
                              statename, ".csv")) %>% 
    filter(!is.na(species)) %>% 
    filter(class == "Aves") %>% 
    filter(lubridate::year(eventDate) < 2020)
  
  if (file.exists("raw_data/species_codes_tbl.csv")) {
    species_codes_tbl <- read_csv("raw_data/species_codes_tbl.csv")
  } else {
    species_codes_tbl <- data.frame(
      name_clean = character(0),
      SCIENTIFIC.NAME = character(0),
      specID = character(0)
    )
  }
  
  # Get full list of species
  eBird_species <- ebird_dat %>% 
    distinct(name_clean, SCIENTIFIC.NAME) %>% 
    left_join(species_codes_tbl) %>% 
    filter(!grepl("sp\\.", SCIENTIFIC.NAME)) %>% 
    filter(!grepl("\\/", SCIENTIFIC.NAME)) %>% 
    filter(!grepl("\\(", SCIENTIFIC.NAME)) %>% 
    filter(!grepl("\\)", SCIENTIFIC.NAME)) %>% 
    filter(!grepl(" x ", SCIENTIFIC.NAME))
  if (any(is.na(eBird_species$specID))) {
    eBird_species$SCIENTIFIC.NAME[is.na(eBird_species$specID)] <- 
      manual_get_ids(eBird_species$SCIENTIFIC.NAME[is.na(eBird_species$specID)], 
                     "ott")
  }
  
  iNat_species <- inat_dat %>% 
    distinct(species) %>% 
    left_join(species_codes_tbl, by = c("species" = "SCIENTIFIC.NAME"))
  if (any(is.na(iNat_species$specID))) {
    iNat_species$specID[is.na(iNat_species$specID)] <- 
      manual_get_ids(iNat_species$species[is.na(iNat_species$specID)], 
                     "ott")
  }
  
  ebird_dat <- left_join(ebird_dat, eBird_species) %>% 
    filter(!is.na(specID))
  inat_dat <- left_join(inat_dat, iNat_species) %>% 
    filter(!is.na(specID))
  
  # Step 2. Associate observations with hexes.
  # Extract cell IDs to each observation.
  ebird_coords <- ebird_dat[, c("LONGITUDE", "LATITUDE")] %>% distinct() %>% 
    filter(!is.na(LONGITUDE) & !is.na(LATITUDE))
  ebird_pts <- SpatialPoints(ebird_coords, 
                             proj4string = CRS("+proj=longlat +datum=WGS84"))
  ebird_pts <- spTransform(ebird_pts, CRS("+init=epsg:2163 +datum=WGS84"))
  
  ebird_coords$associated_hex <- (ebird_pts %over% gridObj)$gridID
  
  ebird_dat_hexed <- left_join(ebird_dat, ebird_coords, by = c("LATITUDE", "LONGITUDE"))
  ebird_complete_hexed <- ebird_dat_hexed %>% filter(ALL.SPECIES.REPORTED == 1)
  ebird_incomplete_hexed <- ebird_dat_hexed %>% filter(ALL.SPECIES.REPORTED == 0)
  
  
  inat_coords <- distinct(inat_dat[, c("decimalLongitude", "decimalLatitude")])
  inat_pts <- SpatialPoints(inat_coords, 
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  inat_pts <- spTransform(inat_pts, CRS("+init=epsg:2163 +datum=WGS84"))
  
  inat_coords$associated_hex <- (inat_pts %over% gridObj)$gridID
  inat_dat_hexed <- left_join(inat_dat, inat_coords, by = c("decimalLongitude", "decimalLatitude"))
  
  inat_counts_hexed <- inat_dat_hexed %>% 
    mutate(species = specID) %>% 
    count(species, associated_hex)
  ebird_complete_counts_hexed <- ebird_complete_hexed %>% 
    mutate(species = specID) %>% 
    count(species, associated_hex) 

  # Step 3. Save this hexed information.
  this_hexes <- unique(c(inat_counts_hexed$associated_hex,
                         ebird_complete_counts_hexed$associated_hex))
  this_hexes <- this_hexes[!is.na(this_hexes)]
  
  for (hex in this_hexes) {
    this_dat <- lapply(
      list(inat_counts_hexed, 
           ebird_complete_counts_hexed),
      function(x) {
        x %>% filter(associated_hex == hex)
      }
    )
    
    this_specs <- unique(c(this_dat[[1]]$species,
                           this_dat[[2]]$species))
    
    this_file <- file.path(outputPath, paste0("hex_", hex, outputSuffix, ".csv"))
    if (file.exists(this_file)) {
      all_counts_df <- read_csv(this_file)
    } else {
      all_counts_df <- data.frame(species = this_specs,
                                  hex = hex) %>% 
        mutate(species = as.character(species)) %>% 
        mutate(inat_n = 0,
               ebird_comp_n = 0)
    }
    
    if (any(!this_specs %in% all_counts_df$species)) {
      missing_specs_df <- data.frame(
        species = this_specs[!this_specs %in% all_counts_df$species],
        hex = hex,
        inat_n = 0,
        ebird_comp_n = 0
      )
      
      all_counts_df <- bind_rows(all_counts_df, missing_specs_df)
    }
    
    indices <- match(this_dat[[1]]$species, all_counts_df$species)
    all_counts_df$inat_n[indices] <- 
      all_counts_df$inat_n[indices] + 
      this_dat[[1]]$n
    
    indices <- match(this_dat[[2]]$species, all_counts_df$species)
    all_counts_df$ebird_comp_n[indices] <- 
      all_counts_df$ebird_comp_n[indices] + 
      this_dat[[2]]$n
    
    write_csv(all_counts_df, this_file)
    
  }
  
  stateProcInfo$done[stateProcInfo$State == statename] <- TRUE
  write_csv(stateProcInfo, "intermediate_data/continentHexRates/stateProcInfo.csv")
}
