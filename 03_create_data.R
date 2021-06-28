# 03_create_data.R
# Author: Benjamin R. Goldstein
# Date: 6/3/2021

# This file gets from the aggregated hexagonal counts for each state to the
# input needed for species-level analysis, including filtering species of
# interest, gathering species-level trait data, and reformatting

schuetz_johnston_input_file <- stop("Point to the supplementary .xlsx file from Schuetz and Johnston (2019) here")
eltontraits_input_file <- stop("Point to the EltonTraits .txt (maybe 'BirdFuncDat.txt'? here.")


library(tidyverse)
library(taxalight)
library(rgdal)

source("contHexDetectionRate_fn.R")

##### Prepare inputs #####
hex_grid <- readOGR(dsn = "intermediate_data/continent_hexes", layer = "continent_hexes")

hex_grid_dat <- cbind(hex_grid@data, coordinates(hex_grid))
colnames(hex_grid_dat) <- c("hex", "state", "lon", "lat")



ebird_specs <- read_csv("raw_data/ebird_species_tbl.csv") %>% 
  filter(!grepl("sp\\.", SCIENTIFIC.NAME)) %>% 
  filter(!grepl("\\/", SCIENTIFIC.NAME)) %>% 
  filter(!grepl("\\(", SCIENTIFIC.NAME)) %>% 
  filter(!grepl("\\)", SCIENTIFIC.NAME)) %>% 
  filter(!grepl(" x ", SCIENTIFIC.NAME))

ebird_specs$specID <- manual_get_ids(ebird_specs$SCIENTIFIC.NAME, db = "ott")

##### Aggregate trait data #####
Schuetz_data <- readxl::read_xlsx(schuetz_johnston_input_file, sheet = 3) %>%
  select(
    scientific.name, Family, Order, log.10.mass,
    resident, crest, max.color.contrast
  )

Schuetz_data <- Schuetz_data %>% 
  mutate(species_code = 
           manual_get_ids(scientific.name, db = "ott")) %>% 
  filter(!is.na(species_code))

EltonTraits <- read_tsv(eltontraits_input_file) %>% 
  mutate(species_code = manual_get_ids(Scientific, db = "ott")) %>% 
  filter(species_code %in% Schuetz_data$species_code)



trait_data <- left_join(Schuetz_data, EltonTraits, by = "species_code") %>% 
  select(species_code, scientific.name, Family, Order, crest,
         max.color.contrast, PelagicSpecialist, `BodyMass-Value`) %>% 
  rename(species = scientific.name, family = Family, order = Order,
         pelagic = PelagicSpecialist, mass = `BodyMass-Value`) %>% 
  filter(!pelagic) %>% 
  filter(species_code %in% ebird_specs$specID)



##### Prepare count data #####
# Read in some important input files
statecode_df <- read_csv("raw_data/statecodes.csv") %>% 
  filter(!(Abbreviation %in% c("AK", "HI", "DC")))

alldat <- read_csv("intermediate_data/continentHexRates/continentHexInfo.csv") %>% 
  filter(species %in% trait_data$species_code)


# Double-check that no hexes don't have data for either species
sum(alldat$inat_n == 0 & alldat$ebird_comp_n == 0)

# Clean the names of species
all_specs <- alldat %>% 
  count(species) %>% 
  filter(!grepl("sp\\.", species)) %>% 
  filter(!grepl("\\/", species)) %>% 
  filter(!grepl("\\(", species)) %>% 
  filter(!grepl("\\)", species)) %>% 
  filter(!grepl(" x ", species))

# Redo totals since we dropped species
new_totals <- 
  alldat %>% 
  group_by(hex) %>% 
  summarize(ebd_total_NEW = sum(ebird_comp_n),
            inat_total_NEW = sum(inat_n))

# Combine all the data into a df
good_dat <- left_join(alldat, new_totals, by = "hex") %>% 
  select(hex, species, ebird_comp_n, inat_n, ebd_total_NEW, inat_total_NEW) %>% 
  rename(ebird_n = ebird_comp_n,
         ebird_total = ebd_total_NEW,
         inat_total = inat_total_NEW) %>% 
  pivot_longer(cols = c(ebird_n, inat_n),
             names_to = "dataset") %>% 
  rename(success = value) %>% 
  mutate(total = ifelse(dataset == "inat_n", inat_total, ebird_total)) %>% 
  select(dataset, species, hex, total, success, ) %>% 
  left_join(hex_grid_dat) %>% 
  left_join(ebird_specs, by = c("species" = "specID")) %>% 
  mutate(specID = species, species = SCIENTIFIC.NAME) %>% 
  select(-SCIENTIFIC.NAME)

# Get list of hexes that have any obs in each dataset
hexes_in_inat <- good_dat %>% 
  filter(total > 0, dataset == "inat_n") %>% 
  distinct(hex)
hexes_in_ebird <- good_dat %>% 
  filter(total > 0, dataset == "ebird_n") %>% 
  distinct(hex)

# Make sure we only have hexes that actually have data (defensive)
good_dat <- good_dat %>% 
  filter(hex %in% hexes_in_inat$hex & hex %in% hexes_in_ebird$hex)


# Only accept species in at least 100 hexes in either dataset:
spec_hexes <- good_dat %>% 
  filter(dataset == "ebird_n") %>% 
  count(species) %>% 
  filter(n >= 100) %>% 
  left_join(ebird_specs[, c("SCIENTIFIC.NAME", "specID")],
            by = c("species" = "SCIENTIFIC.NAME"))

# Filter trait data down to only those species we want
trait_data <- trait_data %>% 
  filter(species_code %in% spec_hexes$specID) %>% 
  filter(!duplicated(species))  # This is probably unnecessary
good_dat <- good_dat %>% 
  filter(species %in% spec_hexes$species)

#### eBird checklist overall rates per species ####

# Number of hexes each species was observed in, now in a separate df
nhex_per_spec <- good_dat %>% 
  filter(dataset == "ebird_n", success > 0) %>% 
  distinct(species, hex) %>% 
  count(species) %>% 
  rename(nhex = n)

# Number of total observations of each species in each dataset
all_state_files <- list.files("../eBird_Data/subsets-2020/state_rdbs/",
                              pattern = "checklist", full.names = T)
checklists_per_state <- lapply(all_state_files, function(x) {
  nrow(read_csv(x))
})
total_cls <- sum(unlist(checklists_per_state))

spec_rates <- good_dat %>% 
  group_by(species) %>%
  summarize(all_success = sum(success)) %>% 
  mutate(rate = all_success / total_cls)

rarity_df <- left_join(spec_rates, nhex_per_spec, by = "species")

#### Write files ####
write_csv(rarity_df, "intermediate_data/species_rarity.csv")
write_csv(trait_data, "intermediate_data/trait_data.csv")
write_csv(good_dat, "intermediate_data/hex_data.csv.gz")
write_csv(ebird_specs, "intermediate_data/ebird_species_tbl.csv")

