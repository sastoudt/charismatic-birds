# getPhyloTreeMat.R
# Author: Sara Stoudt
# Date: 6/15/2021
# Updated: 10/24/2021


### read in necessary data ####

allData <- read.csv("intermediate_data/GAM_diffs.csv")
dim(allData) ## 482

toRemove <- c("Chestnut-sided_Warbler", "Double-crested_Cormorant", "Swainsons_Thrush", "Red_Crossbill", "Black-throated_Blue_Warbler", "Western_Kingbird", "Nashville_Warbler", "Ash-throated_Flycatcher", "Cassins_Kingbird") ## convergence issues in the first stage

allData <- allData[-which(allData$common_name %in% toRemove), ]
dim(allData) ## 473


numHexes <- read.csv("intermediate_data/species_rarity.csv")
dim(numHexes) # 482

trait_dat <- read_csv("intermediate_data/trait_data.csv")
trait_dat <- trait_dat %>% select(-species)
ebird_specs <- read_csv("intermediate_data/ebird_species_tbl.csv")

s1_dat <- allData %>%
  left_join(ebird_specs, by = c("species" = "SCIENTIFIC.NAME")) %>%
  left_join(trait_dat, by = c("specID" = "species_code")) %>%
  filter(!is.na(common_name), !is.na(family))
dim(s1_dat) ## 472

s1_dat <- merge(s1_dat, numHexes, by.x = "species", by.y = "species")
dim(s1_dat) ## 472

#### get tree - only run manually once ####
## https://birdtree.org/subsets/

if (!file.exists("intermediate_data/phylo_covSmall.rds")) {
  treeNames <- read.csv("intermediate_data/BLIOCPhyloMasterTax.csv")
  ## created from species list on https://birdtree.org/subsets/

  which(treeNames$specID %in% unique(s1_dat$specID)) %>% length() ## 472

  unique(s1_dat$species) %>% length() ## 472

  write.csv(treeNames$Scientific[which(treeNames$specID %in% unique(s1_dat$specID))], "birdsWeNeed.csv", row.names = F) ## copy and paste into https://birdtree.org/subsets/ to get tree
  ## Ericson All Species: a set of 10000 trees with 9993 OTUs each, create 100 trees

  your_tree_results <- "tree-pruner-ed3d7a63-57e3-4fc3-9780-21e8bfe6d1dc" ## replace with your own
  ## put in intermediate_phylo_materials folder
  tree <- read.nexus(paste("intermediate_data/", your_tree_results, "/output.nex", sep = ""))

  # http://blog.phytools.org/2016/03/method-to-compute-consensus-edge.html
  branchLengths <- consensus.edges(tree, if.absent = "ignore") ## will take a little bit

  phylo_cov <- vcv(branchLengths) ## covariance matrix
  dim(phylo_cov) ## 472 x 472
  ## save this

  saveRDS(phylo_cov, "intermediate_data/phylo_covSmall.rds")
} else {
  phylo_cov <- readRDS("intermediate_data/phylo_covSmall.rds")
}