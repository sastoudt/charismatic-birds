# 05_predictAllGAM_hiK.R
# Author: Benjamin R. Goldstein, Sara Stoudt
# Date: 6/3/2021

# This file contains the same procedure as 04_predictAllGAM.R but is set up
# to run with higher values of k only on species that haven't yet produced
# adequate results. k can be manually set to iteratively increase it.

library(mgcv)
library(tidyverse)
library(taxadb)
library(parallel)
library(tryCatchLog)
library(futile.logger)
library(R.utils)
library(nimble)

flog.appender(appender.file("warning.log"))
flog.threshold(WARN)

source("helper_code/predictAllGAM_fn.R")

# We create some directories. lpmtx will store some GAM outputs that we won't
# actually use but needed for reproducing the surfaces. ind_results is where
# we'll store the results of interest. Warnings is a folder where warnings and
# errors will be written (if they occur--they shouldn't if unless we're timing
# out). If an individual species throws an error, the error will be saved to a
# file and computation will continue
if (!dir.exists("intermediate_data/warnings")) dir.create("intermediate_data/warnings")
if (!dir.exists("intermediate_data/ind_results")) dir.create("intermediate_data/ind_results")
if (!dir.exists("intermediate_data/lpmtx")) dir.create("intermediate_data/lpmtx")

nCores <- 1 # recommend set to 1, since mgcv operates across threads
timeout <- Inf # Stop after a while for each dataset? (Helps catch super slow species)

# Grab inputs
hex_data <- read_csv("intermediate_data/hex_data.csv.gz")
trait_data <- read_csv("intermediate_data/trait_data.csv") %>% 
  select(-species)
ebird_specs <- read_csv("intermediate_data/ebird_species_tbl.csv")

output_file <- "intermediate_data/GAM_diffs.csv"

target_species <- unique(hex_data$species) # All species

## store p-value


# 1
if (nCores > 1) { # If parallel
  cl <- makeCluster(nCores)
  ignored <- clusterEvalQ(cl = cl, {
    library(mgcv)
    library(tidyverse)
    library(tryCatchLog)
    library(futile.logger)
    library(R.utils)
    flog.appender(appender.file("warning.log"))
    flog.threshold(WARN)
  })
  clusterExport(cl = cl, varlist = c("process_spec_diff", "hex_data", "ebird_specs"))
  
  time.taken <- system.time(
    result_list <-
      parLapply(cl = cl, X = target_species, fun = function(x) {
        tryLog(process_spec_diff(x, hex_data = hex_data, ebird_specs = ebird_specs, timeout = timeout))
      })
  )
} else { # ^ End if parallel | If one thread (recommended) v

  for (this_species in target_species) {
    # Get the common name to make a nice filename
    this_name_clean <- ebird_specs$name_clean[ebird_specs$SCIENTIFIC.NAME == this_species]
    
    # Set up warning printing and output filepath
    flog.appender(appender.file(paste0("intermediate_data/warnings/", this_name_clean, "_warning.log")))
    this_output_file <- paste0("intermediate_data/ind_results/res_", this_name_clean, ".csv")
    
    # Only run if we haven't done this species yet
    if (!file.exists(this_output_file)) {
      cat("Processing", this_species, "-", this_name_clean, "\n")
      
      # We do write a file and use the Rscript command to make sure that we're
      # running each model in a fresh R session. The idea here is to encourage
      # the garbage collector to do its job; before using this method we were
      # experiencing some hanging/memory overflow issues
      cat(paste0(c('
       library(mgcv)
       library(tidyverse)
       library(tryCatchLog)
       library(futile.logger)
       library(R.utils)
       
       flog.appender(appender.file("warning.log"))
       flog.threshold(WARN)
      
       hex_data <- read_csv("intermediate_data/hex_data.csv.gz")
       ebird_specs <- read_csv("intermediate_data/ebird_species_tbl.csv")
       source("predictAllGAM_fn.R")
       tryLog(process_spec_diff("',
        this_species,
       '", hex_data = hex_data, ebird_specs = ebird_specs, timeout =',
        timeout,
       ', maxK = 20))'), collapse = ""), file = "tempfile.R")
    
      sink("captured.txt")
      captured <- system("Rscript tempfile.R")
      sink()
    }
  }
}

# Write the results (not hugely necessary)
result_df <- do.call(rbind, lapply(list.files("intermediate_data/ind_results", full.names = T), read_csv))
write_csv(result_df, output_file)

# How many have converged?
hist(result_df$ebird_GAM_pval)
sum(result_df$ebird_GAM_pval < 0.1)
mean(result_df$ebird_GAM_pval < 0.1)

hist(result_df$inat_GAM_pval)
sum(result_df$inat_GAM_pval < 0.1)
mean(result_df$inat_GAM_pval < 0.1)


# How many succeeded?
files_ind <- list.files("intermediate_data/ind_results", full.names = T)
fitspecs <- substr(files_ind, 35, nchar(files_ind) - 4)

# How many gave errors or warnings?
files_warn <- list.files("intermediate_data/warnings", full.names = T)
warnspecs <- substr(files_warn, 28, nchar(files_warn) - 12)

