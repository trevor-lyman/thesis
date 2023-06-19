# Trevor Pettit
# May 28, 2023
# webs.R

# The purpose of this script is to...

# Step 1: Load packages
library(dplyr)
library(ggplot2)
library(soilfoodwebs)
library(corrplot)
library(tidyr)
library(data.table)
library(reshape2)

# Step 2: Read Data
temp2 <- read.csv("Outputs/webs.csv")
model_biomasses <- read.csv("Outputs/model_biomasses.csv")
sample_metadata <- read.csv("Data/sample metadata v2.csv")
ExRes_imat <- as.data.frame(read.csv("Data/ExRes_imat_11nodes.csv"))
# rows denote consumers, cols denote resources
# 1 denotes the presence of a trophic interaction, 0 denotes the absence 
ExRes_prop_control <- 
  as.data.frame(read.csv("Data/ExRes_control_params_11nodes.csv")) 
# control model parameters 
# includes death rate, assimilation efficiency, production efficiency, C:N 
# for each node at 12 deg C (ambient temperature tx)
ExRes_prop_warming <- 
  as.data.frame(read.csv("Data/ExRes_warming_params_11nodes.csv")) 
# warming model parameters
# includes death rate, assimilation efficiency, production efficiency, C:N 
# for each node at 12 deg C (ambient temperature tx)

# Step 3: Set up 'prop' tabel to draw model parameters from
prop <- as.data.frame(ExRes_prop_control$ID)
colnames(prop) <- "Node_ID"
prop <- prop %>%
  mutate(CN = ExRes_prop_control$CN) %>%
  mutate(a = ExRes_prop_control$a) %>%
  mutate(d_12C = ExRes_prop_control$d) %>%
  mutate(d_20C = ExRes_prop_warming$d) %>%
  mutate(p_12C = ExRes_prop_control$p) %>%
  mutate(p_20C = ExRes_prop_warming$p)

# Step 4: Calculate Flux
# 4.1: Set up data frame with cols for each model
ExRes_biomasses <- as.data.frame(model_biomasses)

# 4.2: Aggregate T1 sample names
names <- sample_metadata$Sample_ID

# 4.3: Check formatting of each df
ExRes_imat$X = NULL
ExRes_imat = as.matrix(ExRes_imat)

# Step 4: Food Web Figs
# 8.1: Tidy Data
temp3 <- merge(x = temp2, y = sample_metadata, by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5")
temp3$moisture_tx <- as.factor(temp3$moisture_tx)
temp3$temp_tx <- as.factor(temp3$temp_tx)

# 8.2: Calculate mean biomass and SE by node by treatment
mean_biomasses_temp <- temp3 %>% 
  # grab tidied data
  group_by(interaction(moisture_tx, temp_tx, destructive_time)) %>%
  # group by treatment 
  summarize(mesostigs = mean(mesostigs_biomass),
            # summarize by means
            zerconidae = mean(zerconidae_biomass),
            prostig_astig = mean(prostig_astig_biomass),
            juv_oribatids = mean(juv_oribatids_biomass),
            oribatids = mean(oribatids_biomass),
            collembola = mean(collembola_biomass),
            n.predator = mean(n.predator_biomass),
            n.bacterivore = mean(n.bacterivore_biomass),
            n.fungivore = mean(n.fungivore_biomass),
            n.omnivore = mean(n.omnivore_biomass),
            basal_resource = mean(basal_resource)
  )

# 8.3: Transpose Data
mean_biomasses <- as.data.frame(t(mean_biomasses_temp)) %>%
  slice(-1) # transpose means
colnames(mean_biomasses) <- 
  mean_biomasses_temp$`interaction(moisture_tx, temp_tx, destructive_time)`
# fix colnames

# for some reason transposing these data coerces them as characters
mean_biomasses$Ambient.12.T1 <- as.numeric(mean_biomasses$Ambient.12.T1) 
mean_biomasses$High.12.T1 <- as.numeric(mean_biomasses$High.12.T1)
mean_biomasses$Ambient.20.T1 <- as.numeric(mean_biomasses$Ambient.20.T1)
mean_biomasses$High.20.T1 <- as.numeric(mean_biomasses$High.20.T1)
mean_biomasses$Ambient.12.T2 <- as.numeric(mean_biomasses$Ambient.12.T2) 
mean_biomasses$High.12.T2 <- as.numeric(mean_biomasses$High.12.T2)
mean_biomasses$Ambient.20.T2 <- as.numeric(mean_biomasses$Ambient.20.T2)
mean_biomasses$High.20.T2 <- as.numeric(mean_biomasses$High.20.T2)

# 8.4: Make foodwebs
# Ambient.12.T1
Ambient.12.T1_prop <- ExRes_prop_control %>% 
  mutate(B = mean_biomasses$Ambient.12.T1)
Ambient.12.T1_prop <- 
  Ambient.12.T1_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                         "DetritusRecycling", "isDetritus", "isPlant", 
                         "canIMM")]
Ambient.12.T1 <- list(imat = ExRes_imat, prop = Ambient.12.T1_prop)
# Check some basic properties of the community:
colnames(Ambient.12.T1$imat) == rownames(Ambient.12.T1$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
Ambient.12.T1$prop$a = Ambient.12.T1$prop$a/100
Ambient.12.T1$prop$p = Ambient.12.T1$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
Ambient.12.T1_model <- comana(Ambient.12.T1, mkplot=T, whattoplot = "web", 
                              BOX.SIZE = 0.05,
                              BOX.PROP = 0.3, 
                              # Box proportion (height: width)
                              arrowlog = F, 
                              # Keep it on the normal scale
                              arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# Ambient.20.T1
Ambient.20.T1_prop <- ExRes_prop_warming %>% 
  mutate(B = mean_biomasses$Ambient.20.T1)
Ambient.20.T1_prop <- 
  Ambient.20.T1_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                         "DetritusRecycling", "isDetritus", "isPlant", 
                         "canIMM")]
Ambient.20.T1 <- list(imat = ExRes_imat, prop = Ambient.20.T1_prop)
# Check some basic properties of the community:
colnames(Ambient.20.T1$imat) == rownames(Ambient.20.T1$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
Ambient.20.T1$prop$a = Ambient.20.T1$prop$a/100
Ambient.20.T1$prop$p = Ambient.20.T1$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
Ambient.20.T1_model <- comana(Ambient.20.T1, mkplot=T, whattoplot = "web", 
                              BOX.SIZE = 0.05,
                              BOX.PROP = 0.3, 
                              # Box proportion (height: width)
                              arrowlog = F, 
                              # Keep it on the normal scale
                              arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# High.12.T1
High.12.T1_prop <- ExRes_prop_control %>% 
  mutate(B = mean_biomasses$High.12.T1)
High.12.T1_prop <- 
  High.12.T1_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                      "DetritusRecycling", "isDetritus", "isPlant", 
                      "canIMM")]
High.12.T1 <- list(imat = ExRes_imat, prop = High.12.T1_prop)
# Check some basic properties of the community:
colnames(High.12.T1$imat) == rownames(High.12.T1$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
High.12.T1$prop$a = High.12.T1$prop$a/100
High.12.T1$prop$p = High.12.T1$prop$p/100
# pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
High.12.T1_model <- comana(High.12.T1, mkplot=T, whattoplot = "web", 
                           BOX.SIZE = 0.05,
                           BOX.PROP = 0.3, 
                           # Box proportion (height: width)
                           arrowlog = F, 
                           # Keep it on the normal scale
                           arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# High.20.T1
High.20.T1_prop <- ExRes_prop_warming %>% 
  mutate(B = mean_biomasses$High.20.T1)
High.20.T1_prop <- 
  High.20.T1_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                      "DetritusRecycling", "isDetritus", "isPlant", 
                      "canIMM")]
High.20.T1 <- list(imat = ExRes_imat, prop = High.20.T1_prop)
# Check some basic properties of the community:
colnames(High.20.T1$imat) == rownames(High.20.T1$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
High.20.T1$prop$a = High.20.T1$prop$a/100
High.20.T1$prop$p = High.20.T1$prop$p/100
# pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
High.20.T1_model <- comana(High.20.T1, mkplot=T, whattoplot = "web", 
                           BOX.SIZE = 0.05,
                           BOX.PROP = 0.3, 
                           # Box proportion (height: width)
                           arrowlog = F, 
                           # Keep it on the normal scale
                           arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# Ambient.12.T2
Ambient.12.T2_prop <- ExRes_prop_control %>% 
  mutate(B = mean_biomasses$Ambient.12.T2)
Ambient.12.T2_prop <- 
  Ambient.12.T2_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                         "DetritusRecycling", "isDetritus", "isPlant", 
                         "canIMM")]
Ambient.12.T2 <- list(imat = ExRes_imat, prop = Ambient.12.T2_prop)
# Check some basic properties of the community:
colnames(Ambient.12.T2$imat) == rownames(Ambient.12.T2$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
Ambient.12.T2$prop$a = Ambient.12.T2$prop$a/100
Ambient.12.T2$prop$p = Ambient.12.T2$prop$p/100
# pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
Ambient.12.T2_model <- comana(Ambient.12.T2, mkplot=T, whattoplot = "web", 
                              BOX.SIZE = 0.05,
                              BOX.PROP = 0.3, 
                              # Box proportion (height: width)
                              arrowlog = F, 
                              # Keep it on the normal scale
                              arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# Ambient.20.T2
Ambient.20.T2_prop <- ExRes_prop_warming %>% 
  mutate(B = mean_biomasses$Ambient.20.T2)
Ambient.20.T2_prop <- 
  Ambient.20.T2_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                         "DetritusRecycling", "isDetritus", "isPlant", 
                         "canIMM")]
Ambient.20.T2 <- list(imat = ExRes_imat, prop = Ambient.20.T2_prop)
# Check some basic properties of the community:
colnames(Ambient.20.T2$imat) == rownames(Ambient.20.T2$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
Ambient.20.T2$prop$a = Ambient.20.T2$prop$a/100
Ambient.20.T2$prop$p = Ambient.20.T2$prop$p/100
# pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
Ambient.20.T2_model <- comana(Ambient.20.T2, mkplot=T, whattoplot = "web", 
                              BOX.SIZE = 0.05,
                              BOX.PROP = 0.3, 
                              # Box proportion (height: width)
                              arrowlog = F, 
                              # Keep it on the normal scale
                              arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# High.12.T2
High.12.T2_prop <- ExRes_prop_control %>% 
  mutate(B = mean_biomasses$High.12.T2)
High.12.T2_prop <- 
  High.12.T2_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                      "DetritusRecycling", "isDetritus", "isPlant", 
                      "canIMM")]
High.12.T2 <- list(imat = ExRes_imat, prop = High.12.T2_prop)
# Check some basic properties of the community:
colnames(High.12.T2$imat) == rownames(High.12.T2$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
High.12.T2$prop$a = High.12.T2$prop$a/100
High.12.T2$prop$p = High.12.T2$prop$p/100
# pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
High.12.T2_model <- comana(High.12.T2, mkplot=T, whattoplot = "web", 
                           BOX.SIZE = 0.05,
                           BOX.PROP = 0.3, 
                           # Box proportion (height: width)
                           arrowlog = F, 
                           # Keep it on the normal scale
                           arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# High.20.T2
High.20.T2_prop <- ExRes_prop_warming %>% 
  mutate(B = mean_biomasses$High.20.T2)
High.20.T2_prop <- 
  High.20.T2_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                      "DetritusRecycling", "isDetritus", "isPlant", 
                      "canIMM")]
High.20.T2 <- list(imat = ExRes_imat, prop = High.20.T2_prop)
# Check some basic properties of the community:
colnames(High.20.T2$imat) == rownames(High.20.T2$imat) 
# Rescale a and p to be [0,1] instead of [0,100]
High.20.T2$prop$a = High.20.T2$prop$a/100
High.20.T2$prop$p = High.20.T2$prop$p/100
# pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
High.20.T2_model <- comana(High.20.T2, mkplot=T, whattoplot = "web", 
                           BOX.SIZE = 0.05,
                           BOX.PROP = 0.3, 
                           # Box proportion (height: width)
                           arrowlog = F, 
                           # Keep it on the normal scale
                           arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)
