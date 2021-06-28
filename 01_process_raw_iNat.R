# 01_process_raw_iNat.R
# Author: Sara Stoudt
# Date: 6/15/2021

library(readr)
library(dplyr)


# DOI: https://doi.org/10.15468/dl.7mt5mn 
# Creation Date: 19:59:48 9 March 2021
# Records included: 24638398 records from 1 published datasets
# Compressed data size: 3.9 GB
# Download format: simple tab-separated values (TSV)
# Filter used:
#
#   {
#     "DatasetKey" : [
#       "is iNaturalist Research-grade Observations"
#       ]
#   }

input_file <- stop("Put the filepath to your raw iNaturalist .csv file here")
output_file_dir <- "raw_data/"

all <- read_tsv(input_file)

us <- subset(all, countryCode == "US")
rm(all)


byState <- split(us, us$stateProvince)

## only birds
for (i in 1:length(byState)) {
  data <- byState[[i]]
  birds <- data %>% filter(class == "Aves")
  write.csv(birds, paste(output_file_dir, "iNat_", names(byState)[i], ".csv", sep = ""), row.names = F)
  ## save by state
  print(i)
}
