#### setup ####
attr(phylo_cov, which = "dimnames")[[1]]

treeNames <- read.csv("intermediate_phylo_materials/BLIOCPhyloMasterTax.csv")

dimName <- gsub("_", " ", attr(phylo_cov, which = "dimnames")[[1]])

testM <- merge(cbind.data.frame(id = 1:nrow(phylo_cov), dimName = dimName), treeNames, by.x = "dimName", by.y = "Scientific")

testM2 <- testM %>% arrange(id)


attr(phylo_cov, which = "dimnames")[[1]] <- testM2$specID
attr(phylo_cov, which = "dimnames")[[2]] <- testM2$specID

setdiff(attr(phylo_cov, which = "dimnames")[[1]], s1_dat$specID) ## none
setdiff(s1_dat$specID, attr(phylo_cov, which = "dimnames")[[1]]) ## none

s1_dat$target_species2 <- as.factor(s1_dat$specID)

setdiff(s1_dat$target_species2, attr(phylo_cov, which = "dimnames")[[1]]) ## none
setdiff(attr(phylo_cov, which = "dimnames")[[1]], s1_dat$target_species2) ## none

#### model ####

id <- diag(nrow(s1_dat))
dimnames(id) <- dimnames(phylo_cov)

s1_dat$target_species3 <- s1_dat$target_species2 ## need a dummy variable for two forms of species effect, identity and phylo cov
class(s1_dat$order)
s1_dat$order <- as.factor(s1_dat$order)


s1_dat$scaled_mass <- scale(log(s1_dat$mass))
s1_dat$scaled_color <- scale(s1_dat$max.color.contrast)
s1_dat$scaled_nhex <- scale(log(s1_dat$nhex))
s1_dat$scaled_rate <- scale(log(s1_dat$rate))


fit <- rma.mv(difference ~ scaled_mass + scaled_color + scaled_nhex + scaled_rate,
  diff_SE^2,
  random = list(~ 1 | target_species2, ~ 1 | target_species3),
  R = list(target_species2 = phylo_cov, target_species3 = id),
  Rscale = "cov0", data = s1_dat
)
summary(fit)
