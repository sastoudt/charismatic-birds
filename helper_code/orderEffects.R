# orderEffects.R
# Author: Sara Stoudt
# Date: 6/15/2021

fit <- rma.mv(difference ~ order - 1,
  diff_SE^2,
  random = list(~ 1 | species),
  R = list(target_species3 = id),
  Rscale = "cov0", data = s1_dat
)

summary(fit)

nice_names <- c(
  "Hawks, Eagles, and Kites", "Ducks, Geese, and Swans", "Nightjars",
  "New World Vultures", "Waders and Gulls", "Storks", "Pigeons and Doves",
  "Kingfishers", "Cuckoos", "Falcons", "Gamefowl", "Loons",
  "Gruiformes", "Songbirds", "Pelicans", "Woodpeckers", "Grebes",
  "Parrots", "Owls", "Cormorants"
)

toSave <- cbind.data.frame(order = row.names(fit$beta), beta = fit$beta[, 1], se = fit$se, common_order_name = nice_names)


write.csv(toSave, "intermediate_phylo_materials/orderMetaAnalysisResults.csv", row.names = F)
