# Trevor Pettit
# Sept 14, 2023
# respiration_revised_v2.R

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
  dplyr::select(-c(X.Msgs, Obs., Port., X.Raw, IV.Date)) %>%
  mutate(Lin_Flux_gofCO2perm2pers = Lin_Flux * mol_per_umol * molar_weight_CO2) %>%
  mutate(Lin_Flux_gofCO2perm2peryr = Lin_Flux * mol_per_umol * molar_weight_CO2 * s_per_year) %>%
  mutate(Exp_Flux_gofCO2perm2pers = Exp_Flux * mol_per_umol * molar_weight_CO2) %>%
  mutate(Exp_Flux_gofCO2perm2peryr = Exp_Flux * mol_per_umol * molar_weight_CO2 * s_per_year) %>%
  filter(!Label == "Dummy") %>%
  filter(!Label == "ExRes 5") %>%
  mutate("Sample_ID" = Label) %>%
  dplyr::select(-c(Label, Lin_Flux, Exp_Flux, Lin_Flux_gofCO2perm2pers, Exp_Flux_gofCO2perm2pers))

resp.data.temp2$Lin_Flux_gofCO2perm2peryr[resp.data.temp2$Lin_Flux_gofCO2perm2peryr<0] <- 0

resp.data.temp2$filter <- interaction(resp.data.temp2$Sample_ID, resp.data.temp2$Sampling_Time)

boxplot(resp.data.temp2$Lin_Flux_gofCO2perm2peryr)
max(resp.data.temp2$Lin_Flux_gofCO2perm2peryr)
resp.data.temp3 <- resp.data.temp2%>%filter(!filter == "ExRes 29.W1")
# outlier 1: W1, ExRes 29
boxplot(resp.data.temp3$Lin_Flux_gofCO2perm2peryr)
max(resp.data.temp3$Lin_Flux_gofCO2perm2peryr)
resp.data.temp4 <- resp.data.temp3%>%filter(!filter == "ExRes 13.W1")
# outlier 2: W1, ExRes 13
boxplot(resp.data.temp4$Lin_Flux_gofCO2perm2peryr)

resp.data <- merge(x = sample_metadata, y = resp.data.temp4, by = "Sample_ID") %>%
  dplyr::select(!c(Exp_FluxCV, Exp_Flux_gofCO2perm2peryr))

resp.data2 <- resp.data %>%
  dplyr::select(-c(Lin_FluxCV, filter)) %>%
  pivot_wider(names_from = Sampling_Time, values_from = Lin_Flux_gofCO2perm2peryr) %>%
  group_by(Sample_ID)

resp.data2 <- resp.data2[, c("Sample_ID", "destructive_time", 
                             "temp_tx", "moisture_tx", "block_effect",
                             "W0", "W1", "W2", "W3", "W4", "W5", "W6")]

resp.data2$cum_resp <- log(rowSums(resp.data2[,6:10], na.rm = T)/4)
#W0 to W3

write.csv(resp.data2, "Outputs/respiration_output.csv")

resp.data <- resp.data2 %>% filter(!destructive_time == "T0")

resp.data$destructive_time <- 
  as.factor(resp.data$destructive_time)
resp.data$block_effect <- 
  as.factor(resp.data$block_effect)
resp.data$temp_tx <- 
  as.factor(resp.data$temp_tx)
resp.data$moisture_tx <- 
  as.factor(resp.data$moisture_tx)

resp_plot <- ggplot(resp.data, 
       aes(x=interaction(destructive_time, temp_tx, moisture_tx), 
           y=cum_resp, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
resp_summary <- resp.data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(
    mean = mean(cum_resp), 
    se = sd(cum_resp)/sqrt(length(cum_resp)))

# create LM
resp.data$block_effect <- factor(resp.data$block_effect, ordered = F)

resp_lm <- lm(cum_resp ~ block_effect + (temp_tx * moisture_tx * destructive_time), 
              data = resp.data)

resp_lm2_temp <- lm(cum_resp ~
                    moisture_tx * temp_tx * destructive_time + 
                    moisture_tx * temp_tx +
                    moisture_tx * destructive_time + 
                    temp_tx * destructive_time +
                    moisture_tx +
                    temp_tx +
                    destructive_time +
                    block_effect,
                  data = resp.data)

resp_lm2 <- lm(cum_resp ~
                
                
                
                
                moisture_tx +
                temp_tx +
                destructive_time +
                block_effect,
            data = resp.data)

summary(resp_lm2)

# create aov
# car package -- ANOVA function & type III ANOVAs:
# https://cran.r-project.org/web/packages/car/car.pdf
library(car)
resp_aov <- 
  Anova(resp_lm, type = 'III')

resp_aov

drop1(resp_lm2_temp, ~., test="F")

resp_aov2 <- 
  Anova(resp_lm2, type = 'III')

resp_aov2

# ddredge()? take a look into this for model selection
# currently: backwards model selection
# to add: interpretation of interaction effects
# note: model selection alpha 0.10

# test assumptions
# test assumptions
ggqqplot(residuals(resp_lm2))
shapiro_test(residuals(resp_lm2))

plot(resp_lm2, 1)
as.data.frame(resp.data) %>% 
  levene_test(cum_resp ~ (moisture_tx * temp_tx * destructive_time))


# call everything
resp_plot; resp_summary; summary(resp_lm); resp_aov
