# trimmedMetaAnalysis_noPhylo
# Author: Sara Stoudt
# Date: 10/24/2021

fit_noPhylo <- rma.mv(difference ~ scaled_mass + scaled_color + scaled_nhex + scaled_rate,
                      diff_SE^2,
                      random = list(~ 1 | target_species2, ~ 1 | target_species3),
                      data = s1_dat
)
summary(fit_noPhylo)

