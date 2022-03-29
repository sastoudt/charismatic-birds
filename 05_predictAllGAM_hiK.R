# 05_predictAllGAM_hiK.R
# Author: Benjamin R. Goldstein, Sara Stoudt
# Date: 6/3/2021

# This file contains the same procedure as 04_predictAllGAM.R but is set up
# to run with higher values of k only on species that haven't yet produced
# adequate results. k can be manually set to iteratively increase it.

# See 04 for better commenting as these two scripts are fairly equivalent. The
# main differences are (1) setting k and (2) reading in results so far such that
# we only re-run species that haven't converged

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

if (!dir.exists("intermediate_data/warnings30")) dir.create("intermediate_data/warnings30")
if (!dir.exists("intermediate_data/ind_results30")) dir.create("intermediate_data/ind_results30")
if (!dir.exists("intermediate_data/lpmtx")) dir.create("intermediate_data/lpmtx")

result_df <- rbind(
  do.call(rbind, lapply(list.files("intermediate_data/ind_results", full.names = T), read_csv)),
  do.call(rbind, lapply(list.files("intermediate_data/ind_results30", full.names = T), read_csv))
) %>% 
  arrange(-k) %>% 
  filter(!duplicated(common_name))
write_csv(result_df, "intermediate_data/GAM_diffs.csv")

nCores <- 1
timeout <- Inf
fitK <- 35
overwrite <- T

# Grab inputs. We need this to see if we've hit "max k" yet (numerically fixed
# ceiling equal to sqrt(nhex)).
hex_data <- read_csv("intermediate_data/hex_data.csv.gz")
trait_data <- read_csv("intermediate_data/trait_data.csv") %>% 
  select(-species)
ebird_specs <- read_csv("intermediate_data/ebird_species_tbl.csv")

result_df$maxK <- NA
for (i in 1:nrow(result_df)){
  this_spec_dat <- hex_data %>% filter(species == result_df$species[i], dataset == "ebird_n")
  result_df$maxK[i] <- floor(sqrt(nrow(this_spec_dat)))
}

# Only do those species that we still need to run
target_species <- result_df$species[(result_df$ebird_GAM_pval < 0.1 |
                                    result_df$inat_GAM_pval < 0.1) &
                                    result_df$maxK != result_df$k &
                                    result_df$k < fitK] # All species

# Skip erroring specs (these were excluded for practical reasons, see manuscript
# supplement)--feel free to comment this line out to try them
target_species <- 
  target_species[!target_species %in% c("Phalacrocorax auritus",
                                        "Tyrannus vociferans", 
                                        "Setophaga pensylvanica")]

# 1
if (nCores > 1) {
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
        tryLog(process_spec_diff(x, hex_data = hex_data, ebird_specs = ebird_specs, 
                                 timeout = timeout, maxK = fitK, 
                                 output_path = "intermediate_data/ind_results30/"))
      })
  )
} else {
  # time.taken <- system.time(
  #   result_list <-
  #     lapply(sample(target_species), function(x) {
  #       tryLog(process_spec_diff(x, hex_data = hex_data, ebird_specs = ebird_specs, timeout = timeout, verbose = T))
  #     })
  # )
  
  for (this_species in target_species[length(target_species):1]) {
    this_name_clean <- ebird_specs$name_clean[ebird_specs$SCIENTIFIC.NAME == this_species]
    
    flog.appender(appender.file(paste0("intermediate_data/warnings30/", this_name_clean, "_warning.log")))
    this_output_file <- paste0("intermediate_data/ind_results30/res_", this_name_clean, ".csv")
    
    if (!file.exists(this_output_file) | overwrite) {
      cat("Processing", this_species, "-", this_name_clean, "\n")
      
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
                   timeout, ", maxK = ", fitK, ", overwrite = T, output_path = 'intermediate_data/ind_results30/',",
                   " warn_path = 'intermediate_data/warnings30/'",
                   '))'), collapse = ""), file = "tempfile.R")
      
      sink("captured.txt")
      captured <- system("Rscript tempfile.R")
      sink()
    }
  }
}


result_df_30 <- do.call(rbind, lapply(list.files("intermediate_data/ind_results30", full.names = T), read_csv))
sum(result_df_30$inat_GAM_pval < 0.1)
sum(result_df_30$ebird_GAM_pval < 0.1)

