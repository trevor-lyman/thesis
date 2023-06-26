# Trevor Pettit
# June 25, 2023
# respiration.R

# The purpose of this script is to...

# Step 1: Load packages
library(dplyr)
library(ggplot2)
library(soilfoodwebs)
library(corrplot)
library(tidyr)
library(data.table)
library(reshape2)

sample_metadata <- read.csv("Data/sample metadata v2.csv")

# Data folder
file.path <- "~/Library/CloudStorage/OneDrive-Personal/Desktop/Lindo Lab/ExRes/Data/Respiration/"

# Import names from data folder
file.names <- list.files(path = file.path, all.files=FALSE,
                         full.names=FALSE, ignore.case=FALSE, 
                         include.dirs=FALSE)
file.names <- file.names[-22] # remove 'Txt files' folder

write.csv(file.names, "Data/respiration_file.names.csv")

# Create custom read function
read <- function(input = " "){
  data.in <- read.csv(
    file = input, header=T) 
  data.out <- data.in
  return(data.out)
}

# Loop read function through respiration folder using file.path
for (i in file.names){
  data <- read(input = paste(file.path, i, sep=""))
  assign(i, data)
}

W0 <- bind_rows(Mesocosms_W0_1to15.csv, Mesocosms_W0_16to30.csv, 
                Mesocosms_W0_31to45.csv) %>%
  mutate(Sampling_Time = "W0")

W1 <- bind_rows(Mesocosms_W1_12C_T1.csv, Mesocosms_W1_12C_T2.csv, 
                Mesocosms_W1_20C_T1.csv, Mesocosms_W1_20C_T2.csv) %>%
  mutate(Sampling_Time = "W1")

W2 <- bind_rows(Mesocosms_W2_12C_T1.csv, Mesocosms_W2_12C_T2.csv,
                Mesocosms_W2_20C_T1.csv, Mesocosms_W2_20C_T2.csv) %>%
  mutate(Sampling_Time = "W2")

W3 <- bind_rows(Mesocosms_W3_12C_T1.csv, Mesocosms_W3_12C_T2.csv,
                Mesocosms_W3_20C_T1.csv, Mesocosms_W3_20C_T2.csv) %>%
  mutate(Sampling_Time = "W3")

W4 <- bind_rows(Mesocosms_W4_12C_T2.csv, Mesocosms_W4_20C_T2.csv) %>%
  mutate(Sampling_Time = "W4")

W5 <- bind_rows(Mesocosms_W5_12C_T2.csv, Mesocosms_W5_20C_T2.csv) %>%
  mutate(Sampling_Time = "W5")

W6 <- bind_rows(Mesocosms_W6_12C_T2.csv, Mesocosms_W6_20C_T2.csv) %>%
  mutate(Sampling_Time = "W6")

resp.data.temp <- bind_rows(W0, W1, W2, W3, W4, W5, W6)

write.csv(resp.data.temp, "Data/respiration data.csv")

mol_per_umol <- 10E-6
molar_weight_CO2 <- 44.01 # g/mol
s_per_year <- 3.1536E7

resp.data.temp2 <- resp.data.temp %>%
  select(-c(X.Msgs, Obs., Port., X.Raw, IV.Date)) %>%
  mutate(Lin_Flux_gofCO2perm2pers = Lin_Flux * mol_per_umol * molar_weight_CO2) %>%
  mutate(Lin_Flux_gofCO2perm2peryr = Lin_Flux * mol_per_umol * molar_weight_CO2 * s_per_year) %>%
  mutate(Exp_Flux_gofCO2perm2pers = Exp_Flux * mol_per_umol * molar_weight_CO2) %>%
  mutate(Exp_Flux_gofCO2perm2peryr = Exp_Flux * mol_per_umol * molar_weight_CO2 * s_per_year) %>%
  filter(!Label == "Dummy") %>%
  filter(!Label == "ExRes 5") %>%
  mutate("Sample_ID" = Label) %>%
  select(-c(Label, Lin_Flux, Exp_Flux, Lin_Flux_gofCO2perm2pers, Exp_Flux_gofCO2perm2pers))

resp.data.temp2$Lin_Flux_gofCO2perm2peryr[resp.data.temp2$Lin_Flux_gofCO2perm2peryr<0] <- 0

resp.data.temp3 <- resp.data.temp2 %>%
  group_by(Sample_ID) %>%
  summarize(mean_flux = mean(Lin_Flux_gofCO2perm2peryr),
            se_flux = sd(Lin_Flux_gofCO2perm2peryr)/
              sqrt(length(Lin_Flux_gofCO2perm2peryr)))

resp.data <- merge(x = sample_metadata, y = resp.data.temp2, by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  select(!c(Exp_FluxCV, Exp_Flux_gofCO2perm2peryr))

resp.data_means <- merge(x = sample_metadata, y = resp.data.temp3, by = "Sample_ID") %>%
  filter(!destructive_time == "T0")

resp.data_means$moisture_tx <- as.factor(resp.data_means$moisture_tx)
resp.data_means$temp_tx <- as.factor(resp.data_means$temp_tx)
resp.data_means$destructive_time <- as.factor(resp.data_means$destructive_time)

resp.data2 <- resp.data %>%
  select(-Lin_FluxCV) %>%
  pivot_wider(names_from = Sampling_Time, values_from = Lin_Flux_gofCO2perm2peryr) %>%
  group_by(Sample_ID)

resp.data2 <- resp.data2[, c("Sample_ID", "destructive_time", 
                             "temp_tx", "moisture_tx", "block_effect",
                             "W0", "W1", "W2", "W3", "W4", "W5", "W6")]

write.csv(resp.data2, "Outputs/respiration_output.csv")
write.csv(resp.data_means, "Outputs/mean_respiration_output.csv")

resp_plot <- ggplot(resp.data_means, 
       aes(x=interaction(temp_tx, moisture_tx), 
           y=mean_flux, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
resp_summary <- resp.data_means %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(
    mean = mean(mean_flux), 
    se = sd(mean_flux)/sqrt(length(mean_flux)))

# create aov
resp_aov <- 
  aov(mean_flux ~ (moisture_tx * temp_tx * destructive_time) + block_effect, 
      data = resp.data_means)

# call everything
resp_plot; resp_summary; summary(resp_aov)

# check assumptions
model_resp <- lm(mean_flux ~ 
                     (moisture_tx * temp_tx * destructive_time) 
                   + block_effect, data=resp.data_means)

ggqqplot(residuals(model_resp))
shapiro_test(residuals(model_resp))

plot(model_resp, 1)
resp.data_means %>% levene_test(mean_flux ~ 
                                         (moisture_tx * temp_tx * destructive_time))
