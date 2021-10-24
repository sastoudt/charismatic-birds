# 05_phyloMeta.R
# Author: Sara Stoudt
# Date: 6/15/2021

# This file contains code for running the meta-analysis summary to identify
# associations between species-level traits and the overreporting index with
# uncertainty propagation.

library(tidyverse)
library(taxadb)
library(metafor) ## meta-analysis
library(ape) ## get phylogenetic covariance matrix
library(phytools) ## consensus tree with edges

source("helper_code/getPhyloTreeMat.R") ## has some one time only manual steps, after that automatic
source("helper_code/fullMetaAnalysis.R") ## model with all indices
source("helper_code/fullMetaAnalysis_endanger.R") ## meta-analysis involving endangered status
source("helper_code/fullMetaAnalysis_noPhylo.R") ## full meta-analysis without the phylogenetic error structure
source("helper_code/getPhyloTreeMatTrim.R") ## has some one time only manual steps, after that automatic
source("helper_code/trimmedMetaAnalysis.R") ## model dropping low-information outliers
source("helper_code/trimmedMetaAnalysis_endanger.R") ## trimmed meta-analysis involving endangered status
source("helper_code/trimmedMetaAnalysis_noPhylo.R") ## trimmed meta-analysis without the phylogenetic error structure
source("helper_code/orderEffects.R") ## get the order effects on trimmed data

