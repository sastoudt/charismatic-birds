# trimmedMetaAnalysis_endanger.R
# Author: Sara Stoudt
# Date: 10/24/2021

s1_dat$endangerment <- relevel(as.factor(s1_dat$endangerment), ref = "Least Concern")
fit_end <- rma.mv(difference ~ endangerment,
                  diff_SE^2, random = list(~ 1 | target_species2, ~ 1 | target_species3),
                  R = list(target_species2 = phylo_cov, target_species3 = id ),
                  Rscale = "cov0", data = s1_dat)
summary(fit_end)