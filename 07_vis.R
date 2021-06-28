# 07_vis.R
# Author: Benjamin R. Goldstein, Sara Stoudt
# Date: 6/14/2021

# This file produces graphic visualizations for the manuscript.

# Manually drop those species excluded due to computational issues:
specs_to_drop <- c(
  "Leiothlypis ruficapilla", "Loxia curvirostra", "Catharus ustulatus", "Tyrannus verticalis", 
  "Myiarchus cinerascens", "Setophaga caerulescens", "Setophaga pensylvanica", "Phalacrocorax auritus",
  "Tyrannus vociferans")


# Create visualizations
library(tidyverse)

negcol <- "#FF9000"
poscol <- "#3888ff"
noscol <- "#858585"

# Load data, calculate CIs and uncorrected p-values
allData <- read.csv("intermediate_data/GAM_diffs.csv") %>% 
  filter(!species %in% specs_to_drop) %>% 
  mutate(lb95CI = difference + qnorm(0.025)*diff_SE,
         ub95CI = difference + qnorm(0.975)*diff_SE,
         pval_uncorr = 2*ifelse(pnorm(difference / diff_SE)>0.5, 1-pnorm(difference / diff_SE), pnorm(difference / diff_SE)),
         )

# Adjust p-value for false detection rate
allData$pval <- p.adjust(allData$pval_uncorr,"fdr")

# Code each species with its significance level according to corrected p-value
allData$sigdir <- ifelse(allData$pval < 0.05,
                         ifelse(allData$difference < 0, "Underreported", "Overreported"),
                         "Non-significant")


# Load supplemental data
species_rarity <- read_csv("intermediate_data/species_rarity.csv")
hex_data <- read_csv("intermediate_data/hex_data.csv.gz")
ebird_specs <- read_csv("intermediate_data/ebird_species_tbl.csv")

trait_dat <- read_csv("intermediate_data/trait_data.csv") %>% 
  select(-species)

# get total observations per hex in iNaturalist
total_inat_obs <- hex_data %>% 
  filter(dataset == "inat_n") %>% 
  distinct(hex, total) %>% 
  .$total %>% 
  sum()

# Add trait data
allData <- allData %>%
  left_join(ebird_specs, by = c("species" = "SCIENTIFIC.NAME")) %>%
  left_join(trait_dat, by = c("specID"="species_code")) %>%
  filter(!is.na(common_name), !is.na(family))

# Get order effects and CIs w/ adjusted p-values
orderEffects <- read_csv("orderMetaAnalysisResults.csv") %>% 
  mutate(order = substr(order, 6, nchar(order)),
         ub = beta + 1.96*se,
         lb = beta - 1.96*se,
         pval_uncorr = 2*ifelse(pnorm(beta / se)>0.5, 1-pnorm(beta / se), pnorm(beta / se)))

orderEffects$pval <- p.adjust(orderEffects$pval_uncorr,"fdr")

orderEffects$sigdir <- ifelse(orderEffects$pval < 0.05,
                         ifelse(orderEffects$beta < 0, "Underreported", "Overreported"),
                         "Non-significant")

# Interpretation: species are often those that fall into the trait patterns
# (e.g. peafowl is one of the largest and most beautiful birds, turkeys are
# large) but there are also specific species relationships that defy trends
# (e.g. burrowing owl). So, while trait-based findings are interesting,
# they certainly don't tell the whole story, and having species-by-species
# estimates is very useful


#### Fig. 2: Species index summary plots ####
# Get counts of data in each significance label
sigcounts <- allData %>% 
  mutate(
    sigdir = factor(sigdir, levels = c(
      "Underreported", "Non-significant", "Overreported"
    ))
  ) %>% 
  count(sigdir) %>% 
  mutate(lab = paste0("n = ", n))

# Plot the counts
(countplot <- allData %>% 
  mutate(
    sigdir = factor(sigdir, levels = c(
      "Underreported", "Non-significant", "Overreported"
    ))
  ) %>% 
  ggplot(aes(sigdir)) +
  xlab("") + ylab("") +
  geom_bar(aes(fill = sigdir), show.legend = F) +
  scale_fill_manual(values = c(negcol, "#888888", poscol)) +
  geom_text(data = sigcounts, aes(sigdir, n, label = lab), vjust = "bottom", nudge_y = 2) +
  theme_minimal() +
  ylim(c(0, 375)) +
  labs(tag = "A") +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()))

# Plot the eight most overreported species
(posplot <- allData %>% 
  left_join(species_rarity) %>% 
  mutate(common_name = gsub("_", " ", common_name)) %>% 
  # filter(rate > 0.01) %>% 
  filter(abs(difference) < 10) %>% 
  arrange(-lb95CI) %>%
  .[1:8,] %>% 
  arrange(lb95CI) %>% 
  mutate(common_name = factor(common_name, levels = common_name)) %>%
  ggplot(aes(common_name, difference)) +
  geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI)) +
  geom_point(col = poscol, cex = 3) +
  coord_flip() +
  labs(tag = "C") +
  xlab("") +
  ylab("") +
  # ylab("Overreporting index (95% CI)") +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()))

# Plot the eight most underreported species
(negplot <- allData %>% 
  left_join(species_rarity) %>% 
  mutate(common_name = gsub("_", " ", common_name)) %>% 
  # filter(rate > 0.01) %>% 
  filter(abs(difference) < 10) %>% 
  arrange(ub95CI) %>%
  mutate(common_name = factor(common_name, levels = common_name[length(common_name):1])) %>% 
  .[1:8,] %>% 
  ggplot(aes(common_name, difference)) +
  geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI)) +
  geom_point(col = negcol, cex = 3) +
  theme_minimal() +
  coord_flip() +
  labs(tag = "B") +
  xlab("") +
  ylab("") +
  # ylab("Overreporting index (95% CI)") +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()))

Fig2 <- gridExtra::grid.arrange(grobs = list(negplot, countplot, posplot), 
                                layout_matrix = matrix(c(
                                  NA, 2, 2, NA,
                                  1, 1, 3, 3,
                                  1, 1, 3, 3), byrow = T, nrow= 3
                                ),
                        bottom = "Overreporting index (95% CI)")
ggsave("output/Fig2.jpg", Fig2, width = 8, height = 4)


#### Fig. 3: trait plots ####
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

(massplot <- allData %>% 
  filter(abs(difference) <= 10) %>%
  mutate(
    scaled_mass = scale(log(mass)),
  ) %>% 
  filter(diff_SE <= 5) %>%
  ggplot(aes(scaled_mass, difference, col = sigdir, alpha = sigdir)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI, alpha = sigdir)) +
  scale_color_manual("", values = c("#888888", poscol, negcol)) +
  scale_alpha_manual("", values = c(0.3, 1, 1)) +
  theme_minimal() +
  labs(tag = "A") +
  xlab("Scaled log mass") +
  ylab("Overreporting index"))
trait_legend <- get_legend(massplot)
massplot <- massplot + theme(legend.position = "none", panel.grid.minor = element_blank())

(colorplot <- allData %>% 
  filter(abs(difference) <= 10) %>% 
  mutate(
    scaled_color = scale(max.color.contrast),
  ) %>% 
  filter(diff_SE <= 5) %>% 
    ggplot(aes(scaled_color, difference, col = sigdir, alpha = sigdir)) + 
    geom_point(show.legend = F) +
    scale_color_manual(values = c("#888888", poscol, negcol)) +
    geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI, alpha = sigdir)) +
    scale_alpha_manual("", values = c(0.3, 1, 1)) +
  theme_minimal() +
    labs(tag = "B") +
  xlab("Scaled color contrast") +
    ylab("")+
    theme(axis.text.y = element_blank())+
    theme(legend.position = "none", panel.grid.minor = element_blank())
  #ylab("Overreporting index")
  )


(prevplot <- allData %>% 
    filter(abs(difference) <= 10) %>% 
  left_join(species_rarity) %>% 
  mutate(
    scaled_rate = scale(log(rate))
  ) %>% 
  filter(diff_SE <= 5) %>% 
    ggplot(aes(scaled_rate, difference, col = sigdir, alpha = sigdir)) + 
    geom_point(show.legend = F) +
    scale_color_manual(values = c("#888888", poscol, negcol)) +
    geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI, alpha = sigdir)) +
    scale_alpha_manual("", values = c(0.3, 1, 1)) +
  theme_minimal() +
    labs(tag = "C") +
  xlab("Scaled log prevalence") +
  ylab("Overreporting index")+
  theme(legend.position = "none", panel.grid.minor = element_blank()))

(rangeplot <- allData %>% 
    filter(abs(difference) <= 10) %>% 
    left_join(species_rarity) %>%
    mutate(
      scaled_nhex = scale(log(nhex))
    ) %>% 
    filter(diff_SE <= 5) %>% 
    ggplot(aes(scaled_nhex, difference, col = sigdir, alpha = sigdir)) + 
    geom_point(show.legend = F) +
    scale_color_manual(values = c("#888888", poscol, negcol)) +
    geom_errorbar(aes(ymin = lb95CI, ymax = ub95CI, alpha = sigdir)) +
    scale_alpha_manual("", values = c(0.3, 1, 1)) +
    theme_minimal() +
    labs(tag = "D") +
    xlab("Scaled log range size") +
    ylab("") +
    theme(axis.text.y = element_blank()) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  )


Fig3 <- gridExtra::grid.arrange(
  massplot, colorplot, prevplot, rangeplot,
  right = trait_legend
)
ggsave("output/Fig3.jpg", Fig3, width = 10, height = 6)




#### Fig. 4: Order plot ####

orderEffects_wN <- allData %>% 
  filter(difference > -10) %>% 
  count(order) %>% 
  right_join(orderEffects) %>% 
  mutate(common_order_name_n = paste0(common_order_name, " (n=", n,")") )

(allData %>% 
  left_join(orderEffects_wN[,c("order", "common_order_name_n")]) %>% 
  filter(abs(difference) < 10) %>% 
  ggplot(aes(common_order_name_n, difference)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(aes(ymin = lb95CI, ymax = ub95CI), outlier.shape = NA) +
  geom_pointrange(data = orderEffects_wN, 
                  aes(common_order_name_n, beta, ymin = lb, ymax = ub, col = sigdir)) +
  scale_color_manual("", values = c("#888888", poscol, negcol)) +
  ylab("Overreporting index") +
  xlab("") +
  theme_minimal() +
  coord_flip() +
  ylim(c(-5.5, 5)) + 
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom")) %>% 
  ggsave(filename = "output/Fig4.jpg", width = 6, height = 6)




#### Supplemental figures ####

# Figure: histogram before and after filtering for difference > -10
diff_unfiltered <- allData %>% 
  filter(difference > -200) %>% 
  ggplot() + 
  geom_histogram(aes(difference), binwidth = 0.5) +
  geom_vline(xintercept = -10, col = "darkred", cex = 0.5) +
  theme_minimal() +
  labs(tag = "A") +
  xlab("Overreporting index") +
  ylab("Count")

diff_filtered <- allData %>% 
  filter(difference > -10) %>% 
  ggplot() + 
  geom_histogram(aes(difference), binwidth = 0.5) +
  theme_minimal() +
  labs(tag = "B") +
  xlab("Overreporting index") +
  ylab("Count")

filterplot <- gridExtra::grid.arrange(
  diff_unfiltered,
  diff_filtered,
  nrow = 2
)

ggsave(filename = "output/FigS1.jpg", filterplot, width = 6.5, height = 4.5)



