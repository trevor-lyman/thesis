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

# Step 2: Read Data
biom <- as.data.frame(read.csv("Test/biomass_long.csv", header=TRUE, row.names=1)) %>%
  dplyr::select(c(Sample_ID, total_biomass, meso_biomass, micro_biomass))
rownames(biom) <- biom$Sample_ID
biom <- biom %>%
  dplyr::select(-Sample_ID)

Cmin <- read.csv("Test/cmin_test.csv", header=TRUE) %>%
  dplyr::filter(!Sample_ID == "ExRes 5")
rownames(Cmin) <- Cmin$Sample_ID
Cmin <- Cmin %>%
  dplyr::select(-Sample_ID)

fauna <- as.data.frame(read.csv("Test/biom_wide.csv", header=TRUE, row.names=2)) %>%
  dplyr::select(-X)

metadata <- read.csv("Test/sample metadata v3.csv", header=TRUE) %>%
  dplyr::filter(!Sample_ID == "ExRes 5")
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata %>%
  dplyr::select(-Sample_ID)

all.equal(rownames(Cmin), rownames(metadata))  
# do the two datasets match?
all.equal(rownames(biom), rownames(metadata))  
# do the two datasets match?
all.equal(rownames(fauna), rownames(metadata))  

Table <- metadata$moisture_tx  
Temp <- metadata$temp_tx
Blk <- metadata$block_effect
Time <- metadata$destructive_time
TT <- interaction(metadata$moisture_tx, metadata$temp_tx, metadata$destructive_time)

sp.dist3 <- vegdist(Cmin, method = "bray")
sp.dist3 ### Bray-Curtis will vary between 0-1 (so proportional similarity)

# For Cmin
adonis2(sp.dist3 ~ Table, data = metadata)  
### this tests for differences among water table treatments and is a permutation test using 999 permutations
adonis2(sp.dist3 ~ Temp, data = metadata) 
adonis2(sp.dist3 ~ Temp*Table, data = metadata)
### this tests for differences among temperature treatments
adonis2(sp.dist3 ~ Temp*Table*Time, data = metadata)
adonis2(sp.dist3 ~ Temp*Table*Time + Blk, data = metadata)

### vegan uses the command metaMDS() for non-metric multidimensional scaling
### The metaMDS function automatically transforms data and checks solution robustness - this is done through permutation tests.
mds <- metaMDS(Cmin, dist = "bray")
mds  # Assesses goodness of ordination (stress value) fit
### the stress value for this run is approximately 0.25 on 2 dimensions. (because it is a permutation test, each run provide a slightly different result)
### stress values should be low.  This value is of marginal acceptance for the plot.  A good rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation. 
stressplot(mds)  # Assesses goodness of fit in a stress plot

### Okay, where's the NMDS plot that shows community similarity? (Dude, where's my plot?)
### plot site scores as text (sample labels)
mds.fig <- ordiplot(mds, type = "text", xlab = "NMDS 1", ylab = "NMDS 2")  
# this is plot of NMDS with each sample as unique ID.   But we want by significant treatment (i.e. temp). 
# but I use this plot to decide how to best visualise my data - e.g. length of axes

citation("vegan")
ordiplot(mds, display = "sites", type = "text", xlab = "NMDS 1", ylab = "NMDS 2")  
# this is plot of NMDS with each sample as unique ID.   But we want by significant treatment (i.e. temp). 
# but I use this plot to decide how to best visualise my data - e.g. length of axes

sim1<-with(metadata, simper(Cmin, Table))
summary(sim1)
sim2<-with(metadata, simper(Cmin, temp))
summary(sim2)