# Trevor Pettit
# June 28, 2023
# MANOVA.R

# The purpose of this script is to...

setwd("/Users/trevorpettit/Library/CloudStorage/OneDrive-Personal/Desktop/Lindo Lab/thesis")


# Step 1: Load packages
library(dplyr)
library(ggplot2)
library(soilfoodwebs)
library(corrplot)
library(tidyr)
library(data.table)
library(reshape2)
library(ggpubr)
library(rstatix)
library(ade4)
library(vegan)
library(MASS)
library(car)

# Step 2: Read Data

abundances <- read.csv("Outputs/abundances.csv")

abundances$destructive_time <- as.factor(abundances$destructive_time)
abundances$moisture_tx <- as.factor(abundances$moisture_tx)
abundances$temp_tx <- as.factor(abundances$temp_tx)
abundances$block_effect <- as.factor(abundances$block_effect)

ab.man <- Manova(cbind(total_micro_ab_standardized, 
                       total_meso_ab_standardized, total_ab_standardized) ~ 
                   (destructive_time*moisture_tx*temp_tx) + 
                   block_effect, 
                 data = abundances)
summary.aov(ab.man)

# one-way example: http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance
# two-way example: https://www.youtube.com/watch?v=_m3QQK53QmI&ab_channel=WakjiraTesfahun

test <- cbind(abundances$total_micro_ab_standardized, 
            abundances$total_meso_ab_standardized)

output <- lm(test~(destructive_time*moisture_tx*temp_tx), data = abundances)
Manova(output, type="III")

abundances2 <- vegdist(abundances$total_ab_standardized, method="bray")

# default test by terms

abundances_T1 <- abundances %>% filter(destructive_time == "T1") %>%
  filter(!Sample_ID == "ExRes 5")

abundances_T2 <- abundances %>% filter(destructive_time == "T2") 

abundances.divT1 <- adonis2(abundances_T1$total_ab_standardized ~ 
                            (moisture_tx*temp_tx) + block_effect, data = abundances_T1, permutations = 999, method="bray")

abundances.divT1

abundances.divT2 <- adonis2(abundances_T2$total_ab_standardized ~ 
                              (moisture_tx*temp_tx), data = abundances_T2, permutations = 999, method="bray")

abundances.divT2
