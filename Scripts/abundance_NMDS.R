library(ade4)
library(vegan)
library(MASS)
library(rmarkdown)
library(dplyr)

T1 <- read.csv("Outputs/ab_NMDS_T1.csv", header=TRUE, row.names=2) %>%
  dplyr::select(-X)
dim(T1)

T2 <- read.csv("Outputs/ab_NMDS_T2.csv", header=TRUE, row.names=2) %>%
  dplyr::select(-X)
dim(T2)

metadata <- read.csv("Data/sample metadata v2.csv", header=TRUE, row.names=1)

metadata_T1 <- metadata %>%
  filter(destructive_time == "T1") %>%
  filter(!row_number() %in% c(2)) # remove ExRes 5

metadata_T2 <- metadata %>%
  filter(destructive_time == "T2")

T1 <- T1[order(row.names(T1)), ]
metadata_T1 <- metadata_T1[order(row.names(metadata_T1)), ]
all.equal(rownames(T1), rownames(metadata_T1))  

T2 <- T2[order(row.names(T2)), ]
metadata_T2 <- metadata_T2[order(row.names(metadata_T2)), ]
all.equal(rownames(T2), rownames(metadata_T2))  

### calculate Bray-Curtis distance among samples
sp.dist3_T1 <- vegdist(T1, method = "bray")
sp.dist3_T1 ### Bray-Curtis will vary between 0-1 (so proportional similarity)

sp.dist3_T2 <- vegdist(T2, method = "bray")
sp.dist3_T2 ### Bray-Curtis will vary between 0-1 (so proportional similarity)

adonis2(sp.dist3_T1 ~ moisture_tx, data = metadata_T1)  
adonis2(sp.dist3_T1 ~ temp_tx, data = metadata_T1)  
adonis2(sp.dist3_T1 ~ block_effect, data = metadata_T1)  

adonis2(sp.dist3_T2 ~ moisture_tx, data = metadata_T2)  
adonis2(sp.dist3_T2 ~ temp_tx, data = metadata_T2)  
# moisture appears to be significant in both?

T1.mds <- metaMDS(T1, dist = "bray")
T1.mds
stressplot(T1.mds)

T2.mds <- metaMDS(T2, dist = "bray")
T2.mds
stressplot(T2.mds)

# T1 NMDS
ordiplot(T1.mds, display = "sites", type = "text", xlab = "NMDS 1", ylab = "NMDS 2")  
mds.fig_T1 <- ordiplot(T1.mds, type = "none", xlab = "NMDS 1", ylab = "NMDS 2", 
                    xlim = c(-0.8, 0.8), ylim = c(-0.8, 0.8)) 
points(mds.fig_T1, "sites", pch = 19, cex = 1.5, col = "green4", 
       select = metadata_T1$moisture_tx == "Ambient")
points(mds.fig_T1, "sites", pch = 19, cex = 1.5, col = "peru", 
       select = metadata_T1$moisture_tx == "High")
ordiellipse(T1.mds, metadata_T1$moisture_tx, conf = 0.95, label = TRUE)

# T2 NMDS
ordiplot(T2.mds, display = "sites", type = "text", xlab = "NMDS 1", ylab = "NMDS 2")  
mds.fig_T2 <- ordiplot(T2.mds, type = "none", xlab = "NMDS 1", ylab = "NMDS 2", 
                       xlim = c(-0.8, 0.8), ylim = c(-0.8, 0.8)) 
points(mds.fig_T2, "sites", pch = 19, cex = 1.5, col = "green4", 
       select = metadata_T2$moisture_tx == "Ambient")
points(mds.fig_T2, "sites", pch = 19, cex = 1.5, col = "peru", 
       select = metadata_T2$moisture_tx == "High")
ordiellipse(T2.mds, metadata_T2$moisture_tx, conf = 0.95, label = TRUE)

# try superimposing elipse from T1 NMDS on T2 NMDS plot
mds.fig_T2 <- ordiplot(T2.mds, type = "none", xlab = "NMDS 1", ylab = "NMDS 2", 
                       xlim = c(-0.8, 0.8), ylim = c(-0.8, 0.8)) 
points(mds.fig_T2, "sites", pch = 19, cex = 1.5, col = "green4", 
       select = metadata_T2$moisture_tx == "Ambient")
points(mds.fig_T2, "sites", pch = 19, cex = 1.5, col = "peru", 
       select = metadata_T2$moisture_tx == "High")
ordiellipse(T2.mds, metadata_T2$moisture_tx, conf = 0.95, label = TRUE)
ordiellipse(T1.mds, metadata_T1$moisture_tx, conf = 0.95, col = "grey", label = TRUE)
