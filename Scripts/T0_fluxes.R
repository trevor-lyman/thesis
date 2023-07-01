# Trevor Pettit
# July 1, 2023
# T0_fluxes.R

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

# 22, 23, 39, 44, 45

# ExRes 22
ExRes_22_prop <- ExRes_prop_control %>% 
  mutate("B" = model_biomasses$ExRes.22)
ExRes_22_prop <- ExRes_22_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                 "DetritusRecycling", "isDetritus", "isPlant", 
                                 "canIMM")]
ExRes_22 <- list(imat = ExRes_imat, prop = ExRes_22_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_22$prop$a = ExRes_22$prop$a/100
ExRes_22$prop$p = ExRes_22$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_22_model <- comana(ExRes_22, mkplot=F, whattoplot = "web", 
                        BOX.SIZE = 0.05,
                        BOX.PROP = 0.3, # Box proportion (height: width)
                        arrowlog = F, # Keep it on the normal scale
                        arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 23
ExRes_23_prop <- ExRes_prop_control %>% 
  mutate("B" = model_biomasses$ExRes.23)
ExRes_23_prop <- ExRes_23_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", "isPlant", 
                                   "canIMM")]
ExRes_23 <- list(imat = ExRes_imat, prop = ExRes_23_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_23$prop$a = ExRes_23$prop$a/100
ExRes_23$prop$p = ExRes_23$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_23_model <- comana(ExRes_23, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 39
ExRes_39_prop <- ExRes_prop_control %>% 
  mutate("B" = model_biomasses$ExRes.39)
ExRes_39_prop <- ExRes_39_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", "isPlant", 
                                   "canIMM")]
ExRes_39 <- list(imat = ExRes_imat, prop = ExRes_39_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_39$prop$a = ExRes_39$prop$a/100
ExRes_39$prop$p = ExRes_39$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_39_model <- comana(ExRes_39, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 44
ExRes_44_prop <- ExRes_prop_control %>% 
  mutate("B" = model_biomasses$ExRes.44)
ExRes_44_prop <- ExRes_44_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", "isPlant", 
                                   "canIMM")]
ExRes_44 <- list(imat = ExRes_imat, prop = ExRes_44_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_44$prop$a = ExRes_44$prop$a/100
ExRes_44$prop$p = ExRes_44$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_44_model <- comana(ExRes_44, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

# ExRes 45
ExRes_45_prop <- ExRes_prop_control %>% 
  mutate("B" = model_biomasses$ExRes.45)
ExRes_45_prop <- ExRes_45_prop[, c("X", "ID", "d", "a", "p", "B", "CN", 
                                   "DetritusRecycling", "isDetritus", "isPlant", 
                                   "canIMM")]
ExRes_45 <- list(imat = ExRes_imat, prop = ExRes_45_prop)

# Rescale a and p to be [0,1] instead of [0,100]
ExRes_45$prop$a = ExRes_45$prop$a/100
ExRes_45$prop$p = ExRes_45$prop$p/100
#pdf("PFC.pdf", width = 8, height = 6, bg= "white", colormodel = "cmyk") 
# This will put the plot in the working directory
ExRes_45_model <- comana(ExRes_45, mkplot=F, whattoplot = "web", 
                         BOX.SIZE = 0.05,
                         BOX.PROP = 0.3, # Box proportion (height: width)
                         arrowlog = F, # Keep it on the normal scale
                         arrowsizerange = c(0.1, 30)) 
# What range of line sizes do you want: c(min, max)

names <- c("ExRes_22_model", "ExRes_23_model", "ExRes_39_model", "ExRes_44_model",
           "ExRes_45_model"
)

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

Sample_ID <- c("ExRes 22", "ExRes 23", "ExRes 39", "ExRes 44",
               "ExRes 45"
)
model_outputs_temp <- as.data.frame(cbind(Sample_ID, total_consumption, 
                                          total_Cmin, total_Nmin))

model_outputs <- merge(x = sample_metadata, y = model_outputs_temp, 
                       by = "Sample_ID") 

write.csv(model_outputs, "Outputs/T0_models.csv")
