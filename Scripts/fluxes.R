# Trevor Pettit
# May 26, 2023
# fluxes.R

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
sample_metadata <- read.csv("Data/sample metadata v2.csv")
model_biomasses <- read.csv("Outputs/model_biomasses.csv")
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
rownames(ExRes_imat) = ExRes_imat$X
ExRes_imat$X = NULL
ExRes_imat = as.matrix(ExRes_imat)

# 4.4: Run models
# ExRes 3
ExRes_3_prop <- ExRes_prop_control %>% 
  mutate("B" = model_biomasses$ExRes.3)
ExRes_3_prop <- ExRes_3_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", "isPlant", 
                                 "canIMM")]
ExRes_3 <- list(imat = ExRes_imat, prop = ExRes_3_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_3$prop$a = ExRes_3$prop$a/100
ExRes_3$prop$p = ExRes_3$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_3_model <- comana(ExRes_3, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 5
ExRes_5_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.5)
ExRes_5_prop <- ExRes_5_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", "isPlant", 
                                 "canIMM")]
ExRes_5 <- list(imat = ExRes_imat, prop = ExRes_5_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_5$prop$a = ExRes_5$prop$a/100
ExRes_5$prop$p = ExRes_5$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_5_model <- comana(ExRes_5, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 7
ExRes_7_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.7)
ExRes_7_prop <- ExRes_7_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", "isPlant", 
                                 "canIMM")]
ExRes_7 <- list(imat = ExRes_imat, prop = ExRes_7_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_7$prop$a = ExRes_7$prop$a/100
ExRes_7$prop$p = ExRes_7$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_7_model <- comana(ExRes_7, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 8
ExRes_8_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.8)
ExRes_8_prop <- ExRes_8_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", "isPlant", 
                                 "canIMM")]
ExRes_8 <- list(imat = ExRes_imat, prop = ExRes_8_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_8$prop$a = ExRes_8$prop$a/100
ExRes_8$prop$p = ExRes_8$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_8_model <- comana(ExRes_8, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 11
ExRes_11_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.11)
ExRes_11_prop <- ExRes_11_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_11 <- list(imat = ExRes_imat, prop = ExRes_11_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_11$prop$a = ExRes_11$prop$a/100
ExRes_11$prop$p = ExRes_11$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_11_model <- comana(ExRes_11, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 13
ExRes_13_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.13)
ExRes_13_prop <- ExRes_13_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_13 <- list(imat = ExRes_imat, prop = ExRes_13_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_13$prop$a = ExRes_13$prop$a/100
ExRes_13$prop$p = ExRes_13$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_13_model <- comana(ExRes_13, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 16
ExRes_16_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.16)
ExRes_16_prop <- ExRes_16_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_16 <- list(imat = ExRes_imat, prop = ExRes_16_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_16$prop$a = ExRes_16$prop$a/100
ExRes_16$prop$p = ExRes_16$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_16_model <- comana(ExRes_16, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 17
ExRes_17_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.17)
ExRes_17_prop <- ExRes_17_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_17 <- list(imat = ExRes_imat, prop = ExRes_17_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_17$prop$a = ExRes_17$prop$a/100
ExRes_17$prop$p = ExRes_17$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_17_model <- comana(ExRes_17, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 18
ExRes_18_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.18)
ExRes_18_prop <- ExRes_18_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_18 <- list(imat = ExRes_imat, prop = ExRes_18_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_18$prop$a = ExRes_18$prop$a/100
ExRes_18$prop$p = ExRes_18$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_18_model <- comana(ExRes_18, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 19
ExRes_19_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.19)
ExRes_19_prop <- ExRes_19_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_19 <- list(imat = ExRes_imat, prop = ExRes_19_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_19$prop$a = ExRes_19$prop$a/100
ExRes_19$prop$p = ExRes_19$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_19_model <- comana(ExRes_19, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 20
ExRes_20_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.20)
ExRes_20_prop <- ExRes_20_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_20 <- list(imat = ExRes_imat, prop = ExRes_20_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_20$prop$a = ExRes_20$prop$a/100
ExRes_20$prop$p = ExRes_20$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_20_model <- comana(ExRes_20, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 21
ExRes_21_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.21)
ExRes_21_prop <- ExRes_21_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_21 <- list(imat = ExRes_imat, prop = ExRes_21_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_21$prop$a = ExRes_21$prop$a/100
ExRes_21$prop$p = ExRes_21$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_21_model <- comana(ExRes_21, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 24
ExRes_24_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.24)
ExRes_24_prop <- ExRes_24_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_24 <- list(imat = ExRes_imat, prop = ExRes_24_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_24$prop$a = ExRes_24$prop$a/100
ExRes_24$prop$p = ExRes_24$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_24_model <- comana(ExRes_24, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 25
ExRes_25_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.25)
ExRes_25_prop <- ExRes_25_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_25 <- list(imat = ExRes_imat, prop = ExRes_25_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_25$prop$a = ExRes_25$prop$a/100
ExRes_25$prop$p = ExRes_25$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_25_model <- comana(ExRes_25, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 29
ExRes_29_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.29)
ExRes_29_prop <- ExRes_29_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_29 <- list(imat = ExRes_imat, prop = ExRes_29_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_29$prop$a = ExRes_29$prop$a/100
ExRes_29$prop$p = ExRes_29$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_29_model <- comana(ExRes_29, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 31
ExRes_31_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.31)
ExRes_31_prop <- ExRes_31_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_31 <- list(imat = ExRes_imat, prop = ExRes_31_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_31$prop$a = ExRes_31$prop$a/100
ExRes_31$prop$p = ExRes_31$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_31_model <- comana(ExRes_31, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 33
ExRes_33_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.33)
ExRes_33_prop <- ExRes_33_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_33 <- list(imat = ExRes_imat, prop = ExRes_33_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_33$prop$a = ExRes_33$prop$a/100
ExRes_33$prop$p = ExRes_33$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_33_model <- comana(ExRes_33, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 36
ExRes_36_prop <- ExRes_prop_warming %>% mutate(B = ExRes_biomasses$ExRes.36)
ExRes_36_prop <- ExRes_36_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_36 <- list(imat = ExRes_imat, prop = ExRes_36_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_36$prop$a = ExRes_36$prop$a/100
ExRes_36$prop$p = ExRes_36$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_36_model <- comana(ExRes_36, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 37
ExRes_37_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.37)
ExRes_37_prop <- ExRes_37_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_37 <- list(imat = ExRes_imat, prop = ExRes_37_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_37$prop$a = ExRes_37$prop$a/100
ExRes_37$prop$p = ExRes_37$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_37_model <- comana(ExRes_37, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 43
ExRes_43_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.43)
ExRes_43_prop <- ExRes_43_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_43 <- list(imat = ExRes_imat, prop = ExRes_43_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_43$prop$a = ExRes_43$prop$a/100
ExRes_43$prop$p = ExRes_43$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_43_model <- comana(ExRes_43, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 1
ExRes_1_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.1)
ExRes_1_prop <- ExRes_1_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", 
                                 "isPlant", "canIMM")]
ExRes_1 <- list(imat = ExRes_imat, prop = ExRes_1_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_1$prop$a = ExRes_1$prop$a/100
ExRes_1$prop$p = ExRes_1$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_1_model <- comana(ExRes_1, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 2
ExRes_2_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.2)
ExRes_2_prop <- ExRes_2_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", 
                                 "isPlant", "canIMM")]
ExRes_2 <- list(imat = ExRes_imat, prop = ExRes_2_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_2$prop$a = ExRes_2$prop$a/100
ExRes_2$prop$p = ExRes_2$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_2_model <- comana(ExRes_2, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 4
ExRes_4_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.4)
ExRes_4_prop <- ExRes_4_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", 
                                 "isPlant", "canIMM")]
ExRes_4 <- list(imat = ExRes_imat, prop = ExRes_4_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_4$prop$a = ExRes_4$prop$a/100
ExRes_4$prop$p = ExRes_4$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_4_model <- comana(ExRes_4, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 6 
ExRes_6_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.6)
ExRes_6_prop <- ExRes_6_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", 
                                 "isPlant", "canIMM")]
ExRes_6 <- list(imat = ExRes_imat, prop = ExRes_6_prop)

ExRes_6$prop$a = ExRes_6$prop$a/100
ExRes_6$prop$p = ExRes_6$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_6_model <- comana(ExRes_6, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 9 
ExRes_9_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.9)
ExRes_9_prop <- ExRes_9_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", 
                                 "isPlant", "canIMM")]
ExRes_9 <- list(imat = ExRes_imat, prop = ExRes_9_prop)

ExRes_9$prop$a = ExRes_9$prop$a/100
ExRes_9$prop$p = ExRes_9$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_9_model <- comana(ExRes_9, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 10
ExRes_10_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.10)
ExRes_10_prop <- ExRes_10_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_10 <- list(imat = ExRes_imat, prop = ExRes_10_prop)

ExRes_10$prop$a = ExRes_10$prop$a/100
ExRes_10$prop$p = ExRes_10$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_10_model <- comana(ExRes_10, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 12
ExRes_12_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.12)
ExRes_12_prop <- ExRes_12_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_12 <- list(imat = ExRes_imat, prop = ExRes_12_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_12$prop$a = ExRes_12$prop$a/100
ExRes_12$prop$p = ExRes_12$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_12_model <- comana(ExRes_12, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 14
ExRes_14_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.14)
ExRes_14_prop <- ExRes_14_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_14 <- list(imat = ExRes_imat, prop = ExRes_14_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_14$prop$a = ExRes_14$prop$a/100
ExRes_14$prop$p = ExRes_14$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_14_model <- comana(ExRes_14, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 15
ExRes_15_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.15)
ExRes_15_prop <- ExRes_15_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_15 <- list(imat = ExRes_imat, prop = ExRes_15_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_15$prop$a = ExRes_15$prop$a/100
ExRes_15$prop$p = ExRes_15$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_15_model <- comana(ExRes_15, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 26
ExRes_26_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.26)
ExRes_26_prop <- ExRes_26_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_26 <- list(imat = ExRes_imat, prop = ExRes_26_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_26$prop$a = ExRes_26$prop$a/100
ExRes_26$prop$p = ExRes_26$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_26_model <- comana(ExRes_26, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 27
ExRes_27_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.27)
ExRes_27_prop <- ExRes_27_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_27 <- list(imat = ExRes_imat, prop = ExRes_27_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_27$prop$a = ExRes_27$prop$a/100
ExRes_27$prop$p = ExRes_27$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_27_model <- comana(ExRes_27, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 28
ExRes_28_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.28)
ExRes_28_prop <- ExRes_28_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_28 <- list(imat = ExRes_imat, prop = ExRes_28_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_28$prop$a = ExRes_28$prop$a/100
ExRes_28$prop$p = ExRes_28$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_28_model <- comana(ExRes_28, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 30
ExRes_30_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.30)
ExRes_30_prop <- ExRes_30_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_30 <- list(imat = ExRes_imat, prop = ExRes_30_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_30$prop$a = ExRes_30$prop$a/100
ExRes_30$prop$p = ExRes_30$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_30_model <- comana(ExRes_30, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 32
ExRes_32_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.32)
ExRes_32_prop <- ExRes_32_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_32 <- list(imat = ExRes_imat, prop = ExRes_32_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_32$prop$a = ExRes_32$prop$a/100
ExRes_32$prop$p = ExRes_32$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_32_model <- comana(ExRes_32, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 34
ExRes_34_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.34)
ExRes_34_prop <- ExRes_34_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_34 <- list(imat = ExRes_imat, prop = ExRes_34_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_34$prop$a = ExRes_34$prop$a/100
ExRes_34$prop$p = ExRes_34$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_34_model <- comana(ExRes_34, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 35
ExRes_35_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.35)
ExRes_35_prop <- ExRes_35_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_35 <- list(imat = ExRes_imat, prop = ExRes_35_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_35$prop$a = ExRes_35$prop$a/100
ExRes_35$prop$p = ExRes_35$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_35_model <- comana(ExRes_35, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 38
ExRes_38_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.38)
ExRes_38_prop <- ExRes_38_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_38 <- list(imat = ExRes_imat, prop = ExRes_38_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_38$prop$a = ExRes_38$prop$a/100
ExRes_38$prop$p = ExRes_38$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_38_model <- comana(ExRes_38, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 40
ExRes_40_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.40)
ExRes_40_prop <- ExRes_40_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_40 <- list(imat = ExRes_imat, prop = ExRes_40_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_40$prop$a = ExRes_40$prop$a/100
ExRes_40$prop$p = ExRes_40$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_40_model <- comana(ExRes_40, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 41
ExRes_41_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.41)
ExRes_41_prop <- ExRes_41_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_41 <- list(imat = ExRes_imat, prop = ExRes_41_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_41$prop$a = ExRes_41$prop$a/100
ExRes_41$prop$p = ExRes_41$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_41_model <- comana(ExRes_41, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 42
ExRes_42_prop <- ExRes_prop_control %>% mutate(B = ExRes_biomasses$ExRes.42)
ExRes_42_prop <- ExRes_42_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", 
                                   "isPlant", "canIMM")]
ExRes_42 <- list(imat = ExRes_imat, prop = ExRes_42_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_42$prop$a = ExRes_42$prop$a/100
ExRes_42$prop$p = ExRes_42$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_42_model <- comana(ExRes_42, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# Step 5: Aggregate Data
# 5.1: Set up df and initialize vectors for loop outputs
names <- c("ExRes_3_model", "ExRes_5_model", "ExRes_7_model", "ExRes_8_model",
           "ExRes_11_model", "ExRes_13_model", "ExRes_16_model", 
           "ExRes_17_model", "ExRes_18_model", "ExRes_19_model", 
           "ExRes_20_model", "ExRes_21_model", "ExRes_24_model", 
           "ExRes_25_model", "ExRes_29_model", "ExRes_31_model", 
           "ExRes_33_model", "ExRes_36_model", "ExRes_37_model",
           "ExRes_43_model",
           "ExRes_1_model", "ExRes_2_model", "ExRes_4_model", "ExRes_6_model",
           "ExRes_9_model", "ExRes_10_model", "ExRes_12_model", 
           "ExRes_14_model", "ExRes_15_model", "ExRes_26_model",
           "ExRes_27_model", "ExRes_28_model", "ExRes_30_model",
           "ExRes_32_model", "ExRes_34_model", "ExRes_35_model", 
           "ExRes_38_model", "ExRes_40_model", "ExRes_41_model",
           "ExRes_42_model"
) # subset T1 & T2 models

node_consumption_temp <- list() # initialize
node_Cmin_temp <- list()
node_Nmin_temp <- list()
total_consumption <- list()
total_Cmin <- list()
total_Nmin <- list()

# 5.2: Loop through each model to aggregate total consumption, total C min, 
# and total N min
for (i in names){ # loop to get consumption data, indexed from list as df[1]
  input <- paste(i)
  output <- get(input)
  consumption_temp <- output[1]
  consumption <- as.data.frame(consumption_temp)
  total_consumption_temp <- sum(consumption)
  Cmin_temp <- output[2]
  Cmin <- as.data.frame(Cmin_temp)
  total_Cmin_temp <- sum(Cmin)
  Nmin_temp <- output[3]
  Nmin <- as.data.frame(Nmin_temp)
  total_Nmin_temp <- sum(Nmin)
  node_consumption_temp[[i]] = consumption
  node_Cmin_temp[[i]] = Cmin
  node_Nmin_temp[[i]] = Nmin
  total_consumption[i] = as.data.frame(total_consumption_temp)
  total_Cmin[i] = as.data.frame(total_Cmin_temp)
  total_Nmin[i] = as.data.frame(total_Nmin_temp)
}

# 5.3: Tidy Data
total_consumption <- unlist(total_consumption)
total_Cmin <- unlist(total_Cmin)
total_Nmin <- unlist(total_Nmin)

Sample_ID <- c("ExRes 3", "ExRes 5", "ExRes 7", "ExRes 8",
               "ExRes 11", "ExRes 13", "ExRes 16", 
               "ExRes 17", "ExRes 18", "ExRes 19", 
               "ExRes 20", "ExRes 21", "ExRes 24", 
               "ExRes 25", "ExRes 29", "ExRes 31", 
               "ExRes 33", "ExRes 36", "ExRes 37",
               "ExRes 43",
               "ExRes 1", "ExRes 2", "ExRes 4", "ExRes 6", "ExRes 9",
               "ExRes 10", "ExRes 12", "ExRes 14", "ExRes 15", "ExRes 26",
               "ExRes 27", "ExRes 28", "ExRes 30", "ExRes 32", "ExRes 34",
               "ExRes 35", "ExRes 38", "ExRes 40", "ExRes 41", "ExRes 42"
) # subset T1 & T2 models

model_outputs_temp <- as.data.frame(cbind(Sample_ID, total_consumption, total_Cmin, total_Nmin))

# Step 6: Merge with Sample Metadata  
model_outputs <- merge(x = sample_metadata, y = model_outputs_temp, 
                       by = "Sample_ID") 
total_consumption <- model_outputs %>% select(-total_Cmin, -total_Nmin)
total_Cmin <- model_outputs %>% select(-total_consumption, -total_Nmin)
total_Nmin <- model_outputs %>% select(-total_Cmin, -total_consumption)

# 6.1: Export .csv Files
write.csv(model_outputs, "Outputs/model_outputs.csv")
write.csv(total_consumption, "Outputs/total_consumption.csv")
write.csv(total_Cmin, "Outputs/total_Cmin.csv")
write.csv(total_Nmin, "Outputs/total_Nmin.csv")

# Step 7: Tidy Data
model_outputs <- as.data.frame(model_outputs) %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time)) %>%
  filter(!Sample_ID == "ExRes 5")

model_outputs$total_consumption <- as.numeric(model_outputs$total_consumption)
model_outputs$total_Cmin <- as.numeric(model_outputs$total_Cmin)
model_outputs$total_Nmin <- as.numeric(model_outputs$total_Nmin)

model_outputs$destructive_time <- as.factor(model_outputs$destructive_time)
model_outputs$moisture_tx <- as.factor(model_outputs$moisture_tx)
model_outputs$temp_tx <- as.factor(model_outputs$temp_tx)
model_outputs$block_effect <- as.factor(model_outputs$block_effect)

model_outputs$tx <- droplevels(model_outputs$tx)
key_flux <- levels(model_outputs$tx)
levels(model_outputs$tx) <- 
  c("a", "b", "c", "d", "e", "f", "g", "h")
## a = 12 Ambient T1, b = 20 Ambient T1, c = 12 High T1, d = 20 High T1,
## e = 12 Ambient T2, f = 20 Ambient T2, g = 12 High T2, h = 20 High T2

# Step 8: Consumption
# 8.1: Plot
consumption_plot <- 
  ggplot(model_outputs, aes(x=interaction(temp_tx, moisture_tx), y=total_consumption, 
                            fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# 8.2: Summary Table
consumption_summary <- model_outputs %>%
  group_by(tx) %>%
  summarize(mean = mean(total_consumption),
            se = sd(total_consumption)/
              sqrt(length(total_consumption))) %>%
  mutate(key = key_flux)

# 8.3: AOV
consumption_aov <- 
  aov(total_consumption ~ (moisture_tx * temp_tx * destructive_time) 
      + block_effect, data=model_outputs)

consumption_plot; consumption_summary; summary(consumption_aov)

# Step 9: Cmin
# 9.1: Plot
Cmin_plot <- 
  ggplot(model_outputs, aes(x=interaction(temp_tx, moisture_tx), y=total_Cmin, 
                            fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")


# 9.2: Summary Table
Cmin_summary <- model_outputs %>%
  group_by(tx) %>%
  summarize(mean = mean(total_Cmin),
            se = sd(total_Cmin)/
              sqrt(length(total_Cmin))) %>%
  mutate(key = key_flux)

# 9.3: AOV
Cmin_aov <- 
  aov(total_Cmin ~ (moisture_tx * temp_tx * destructive_time) 
      + block_effect, data=model_outputs)

Cmin_plot; Cmin_summary; summary(Cmin_aov)

# Step 10: Nmin
# 10.1: Plot
Nmin_plot <- 
  ggplot(model_outputs, aes(x=interaction(temp_tx, moisture_tx), y=total_Nmin, 
                            fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# 10.2: Summary Table
Nmin_summary <- model_outputs %>%
  group_by(tx) %>%
  summarize(mean = mean(total_Nmin),
            se = sd(total_Nmin)/
              sqrt(length(total_Nmin))) %>%
  mutate(key = key_flux)

# 10.3: AOV
Nmin_aov <- 
  aov(total_Nmin ~ (moisture_tx * temp_tx * destructive_time) 
      + block_effect, data=model_outputs)

Nmin_plot; Nmin_summary; summary(Nmin_aov)

