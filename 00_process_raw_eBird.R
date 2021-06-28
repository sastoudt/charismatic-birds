# 00_process_raw_eBird.R
# Author: Benjamin R. Goldstein
# Date: 6/3/2021

# This file transforms the eBird raw dataset into a a series of two-file
# relational database storing observation counts and checklist metadata for each
# US state and DC. To replicate, the user needs to download the eBird full
# release dataset as a .txt and indicate its filepath in the "input file" slot.

# With the full repository, this takes a very long time (days to weeks).
# You may want to run it in pieces.

library(tidyverse)
library(USAboundaries)
verbose <- TRUE # Print progress updates?

input_file <- stop("Put the file path to the eBird Basic Dataset .txt here.")

#### Split the full file into workable chunks ####
dir.create("raw_data/raw_split")
system(paste("split -l 2000000 -a 4", input_file, "raw_split/raw_"))

#### For each of these chunks, separate out by state ####

# Get a list of target files and codes
all_target_files <- list.files("raw_split", full.names = TRUE)
statecodes <- USAboundaries::state_codes %>% 
  filter(jurisdiction_type %in% c("state", "district"))

# Extract the raw eBird column names
master_colnames <- colnames(read_tsv(file = all_target_files[[1]], n_max = 10)) %>% 
  gsub(pattern = " ", replacement = ".")

dir.create("raw_data/state_splits")

# Fn to sum counts for duplicates--if any X, then X, otherwise add them
get_sum <- function(xs) {
  if (any(xs == "X")) "X"
  else as.character(sum(as.numeric(xs)))
}

# Function for converting raw splits to state splits
raw_to_stateraws <- function(filepath) {
  filecode <- substr(filepath, 31, 34)
  
  capture <- capture.output(
    rawdf <- read_tsv(filepath, col_names = master_colnames, quote = "",
                      skip = as.numeric(filepath == all_target_files[[1]]),
                      progress = FALSE, col_types = readr::cols(.default = "c"))
  )
  
  rawdf <- rawdf %>% filter(COUNTRY.CODE == "US")
  capture <- lapply(statecodes$state_abbr, function(x) {
    rawdf %>% 
      filter(STATE.CODE == paste0("US-", x)) %>% 
      write_csv(paste0("subsets-Aug2020/raws_by_state/US-", x, "_", filecode, ".csv"))
  })
  
}

# Do all the states
for (i in 1:length(all_target_files)) {
  if (i %% 4 == 0) cat(i, "/", length(all_target_files), "\n")
  capture <- raw_to_stateraws(all_target_files[[i]])
}


#### Create state-level relational databases of eBird ####
for (state_code in state_codes) {
  
  if (verbose) cat("Starting", state_code, "\n")
  
  this_state_files <- list.files("data/state_splits", 
                                 pattern = state_code, full.names = TRUE)
  if (verbose) cat("...Reading files...\n")
  suppressWarnings(
    suppressMessages(
      state_data_list <- lapply(this_state_files[1:200], function(x) {
        read_csv(x, progress = F, 
                 col_types = list(EFFORT.AREA.HA = col_double(),
                                  SAMPLING.EVENT.IDENTIFIER = col_character(),
                                  OBSERVATION.DATE = col_date(),
                                  TIME.OBSERVATIONS.STARTED = col_time(),
                                  OBSERVER.ID = col_character(),
                                  DURATION.MINUTES = col_double(),
                                  EFFORT.DISTANCE.KM = col_double())) %>% 
          distinct(SAMPLING.EVENT.IDENTIFIER, LATITUDE, LONGITUDE, 
                   OBSERVATION.DATE, TIME.OBSERVATIONS.STARTED,
                   OBSERVER.ID, DURATION.MINUTES, EFFORT.DISTANCE.KM, 
                   PROTOCOL.CODE, EFFORT.AREA.HA, ALL.SPECIES.REPORTED,
                   NUMBER.OBSERVERS)
      })
    ))
  
  if (verbose) cat("...Handling checklist data...\n")
  
  state_data_df <- bind_rows(state_data_list)
  rm(state_data_list)
  
  # Get the list of checklists
  state_checklists <- state_data_df %>% 
    distinct(SAMPLING.EVENT.IDENTIFIER, LATITUDE, LONGITUDE, 
             OBSERVATION.DATE, TIME.OBSERVATIONS.STARTED,
             OBSERVER.ID, DURATION.MINUTES, EFFORT.DISTANCE.KM, 
             PROTOCOL.CODE, EFFORT.AREA.HA, ALL.SPECIES.REPORTED,
             NUMBER.OBSERVERS)
  
  # Check non-distinct cases
  duplicates <- state_checklists %>% 
    count(SAMPLING.EVENT.IDENTIFIER) %>% 
    filter(n > 1)
  # Duplicates occur when information about effort area isn't attributed to all obs.
  # The solution is to delete the row with NA effort area
  state_checklists <- state_checklists %>% 
    filter(!(SAMPLING.EVENT.IDENTIFIER %in% 
               duplicates$SAMPLING.EVENT.IDENTIFIER) |
             !is.na(EFFORT.AREA.HA))
  # Check non-distinct cases AGAIN
  duplicates <- state_checklists %>% 
    count(SAMPLING.EVENT.IDENTIFIER) %>% 
    filter(n > 1)
  
  if (nrow(duplicates) > 0) stop(paste0("Duplicates still found. State code ", state_code))
  
  # Now, we have all the checklist info. We simply have to write out that and also write out
  # the observation columns.
  write_csv(state_checklists, 
            paste0("raw_data/", state_code, "_Aug2020_checklist_info.csv"))
  
  if (verbose) cat("...Finished checklists. Reading files for counts...\n")
  
  suppressWarnings(
    suppressMessages(
      state_data_list <- lapply(this_state_files, function(x) {
        read_csv(x, progress = F,
                 col_types = list(COMMON.NAME = col_character(),
                                  SCIENTIFIC.NAME = col_character(),
                                  OBSERVATION.COUNT = col_character())) %>% 
          select(SAMPLING.EVENT.IDENTIFIER, COMMON.NAME,
                 SCIENTIFIC.NAME, OBSERVATION.COUNT)
      })
    ))
  state_data_df <- bind_rows(state_data_list)
  rm(state_data_list)
  
  if (verbose) cat("...Processing counts...\n")
  
  state_obs_info <- state_data_df %>% 
    select(SAMPLING.EVENT.IDENTIFIER, SCIENTIFIC.NAME, 
           COMMON.NAME, OBSERVATION.COUNT) %>% 
    mutate(name_clean = gsub("['.()]", "", 
                             gsub("[ /]", "_", COMMON.NAME)))
  
  species_to_process <- unique(state_obs_info$name_clean)
  
  data_by_species <- list()
  
  if (verbose) cat("...Handling duplicates...\n")
  
  for (i in 1:length(species_to_process)) { # Handle subspecies duplicates
    # writeLines(as.character(i))
    data_by_species[[i]] <- 
      state_obs_info %>% 
      filter(name_clean == species_to_process[i]) %>% 
      select(-COMMON.NAME) %>% 
      group_by(SAMPLING.EVENT.IDENTIFIER, name_clean, SCIENTIFIC.NAME) %>% 
      summarise(total_count = get_sum(OBSERVATION.COUNT), .groups = "drop")
  }
  
  all_data_unique <- bind_rows(data_by_species)
  
  write.csv(all_data_unique, 
            paste0("raw_data/", state_code, "_Aug2020_species_counts.csv"))

}

if (verbose) cat("...Done.\n")


