# getPhyloTreeMatTrim.R
# Author: Sara Stoudt
# Date: 6/15/2021

### setup ###
dropped <- gsub("_", " ", allData$common_name[abs(allData$difference) > 10])
cat(sort(dropped), sep = ", ")
allData <- allData %>% filter(difference < 10 & difference > -10)
dim(allData) ## 424

numHexes <- read.csv("intermediate_data/species_rarity.csv") #
dim(numHexes) # 482

trait_dat <- read_csv("intermediate_data/trait_data.csv")
trait_dat <- trait_dat %>% select(-species)
ebird_specs <- read_csv("intermediate_data/ebird_species_tbl.csv")

s1_dat <- allData %>%
  left_join(ebird_specs, by = c("species" = "SCIENTIFIC.NAME")) %>%
  left_join(trait_dat, by = c("specID" = "species_code")) %>%
  filter(!is.na(common_name), !is.na(family))
dim(s1_dat) ## 424

s1_dat <- merge(s1_dat, numHexes, by.x = "species", by.y = "species")
dim(s1_dat) ## 424

#### get tree - only run manually once ####

## https://birdtree.org/subsets/

if (!file.exists("intermediate_phylo_materials/phylo_covSmallDrop.rds")) {
  treeNames <- read.csv("intermediate_phylo_materials/BLIOCPhyloMasterTax.csv")
  ## created from species list on https://birdtree.org/subsets/


  which(treeNames$specID %in% unique(s1_dat$specID)) %>% length() ## 424


  unique(s1_dat$species) %>% length() ## 424

  write.csv(treeNames$Scientific[which(treeNames$specID %in% unique(s1_dat$specID))], "birdsWeNeedDrop.csv", row.names = F) ## copy and paste into online thing to get tree

  your_tree_results <- "tree-pruner-dcef2944-56b5-4250-9db5-5b40567fcfed" ## replace with your own
  ## put in intermediate_phylo_materials folder
  tree <- read.nexus(paste("intermediate_phylo_materials/", your_tree_results, "/output.nex", sep = ""))

  # http://blog.phytools.org/2016/03/method-to-compute-consensus-edge.html
  branchLengths <- consensus.edges(tree, if.absent = "ignore") ## will take a little bit

  phylo_cov <- vcv(branchLengths) ## covariance matrix
  dim(phylo_cov) ## 424 x 424
  ## save this

  saveRDS(phylo_cov, "intermediate_phylo_materials/phylo_covSmallDrop.rds")
} else {
  phylo_cov <- readRDS("intermediate_phylo_materials/phylo_covSmallDrop.rds")
}
