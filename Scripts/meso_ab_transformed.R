# Trevor Pettit
# Sept 14, 2023
# meso_ab_transformed.R

# The purpose of this script is to...

# Step 1: Load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(soilfoodwebs)
library(corrplot)
library(tidyr)
library(data.table)
library(reshape2)
library(car)
library(rstatix)

# Step 2: Read Data
sample_metadata <- read.csv("Data/sample metadata v2.csv")
extraction_moisture_data <- 
  read.csv("Data/extraction moisture data.csv") 
mesofauna_ab_temp <- 
  read.csv("Data/mesofauna abundances.csv") 
# raw abundance data
meso_body_mass <- 
  read.csv("Data/Barreto et al complete dataset .csv") 
# Carlos' body mass data, in ug of C
extraction_moisture_data <- 
  read.csv("Data/extraction moisture data.csv") 
microfauna_ab_temp <- 
  read.csv("Data/microfauna abundances.csv") %>% # raw abundance data
  mutate(mean_abundance = (ab_1 + ab_2)/2) # average of 2 passes
micro_body_mass <- read.csv("Data/Kamath et al dataset.csv") %>% 
  dplyr::select(-prop_ab) # body mass from Dev's thesis, in ug wwt
micro_prop_ab <- read.csv("Data/Kamath et al dataset.csv") %>% 
  dplyr::select(-body_mass) # proportional abundances from Dev's thesis

sample_metadata$block_effect <- ordered(sample_metadata$block_effect, levels = c("low", "high"))

# Step 3:
mesofauna_ab <- merge(x = mesofauna_ab_temp, y = extraction_moisture_data %>% 
                        filter(extraction_type == "Dry"), 
                      by = "Sample_ID")
microfauna_ab <- merge(x = microfauna_ab_temp, 
                       y = extraction_moisture_data %>% 
                         filter(extraction_type == "Wet"), 
                       by = "Sample_ID") 

# Step 4:
wet_extraction <- extraction_moisture_data %>% 
  filter(extraction_type == "Wet") %>%
  dplyr::select(Sample_ID, extraction_type, dry_weight)
dry_extraction <- extraction_moisture_data %>% 
  filter(extraction_type == "Dry") %>%
  dplyr::select(Sample_ID, extraction_type, dry_weight)

# Step 5:
#cm3_in_m2 <- 100*100*10 
#bulk_soil_density <- 1.63 
#g_per_m2 <- cm3_in_m2*bulk_soil_density*0.25 

cm3_in_m2 <- 100*100*10 # 15 cm deep, less top 5 cm of OM
bulk_soil_density <- 0.075 # g dwt per cm^3?
# https://sis.agr.gc.ca/cansis/publications/manuals/1984-peat/ca1984_2e_report.pdf
g_per_m2 <- cm3_in_m2*bulk_soil_density # cm3_in_m2 * bulk_soil_density 

temp_dry <- merge(x = mesofauna_ab_temp, y = dry_extraction, 
                  by = "Sample_ID") %>%
  mutate(total_meso_ab_standardized = g_per_m2*total_ab/dry_weight) %>%
  dplyr::select(-total_ab)
temp_dry$collembola <- g_per_m2*temp_dry$collembola/temp_dry$dry_weight
temp_dry$juv_oribatids <- g_per_m2*temp_dry$juv_oribatids/temp_dry$dry_weight
temp_dry$oribatids <- g_per_m2*temp_dry$oribatids/temp_dry$dry_weight
temp_dry$prostigs <- g_per_m2*temp_dry$prostigs/temp_dry$dry_weight
temp_dry$astigmata <- g_per_m2*temp_dry$astigmata/temp_dry$dry_weight
temp_dry$zerconidae <- g_per_m2*temp_dry$zerconidae/temp_dry$dry_weight
temp_dry$mesostigs <- g_per_m2*temp_dry$mesostigs/temp_dry$dry_weight
temp_dry$others <- g_per_m2*temp_dry$others/temp_dry$dry_weight

meso_abundances <- merge(x = sample_metadata, y = temp_dry, 
                         by = "Sample_ID") 

write.csv(meso_abundances, "Outputs/meso_abundances.csv")

temp_wet <- microfauna_ab %>% 
  mutate(total_micro_ab_standardized = g_per_m2*mean_abundance/dry_weight)
# note: do not have node level abundances for microfauna as proportional 
# abundances from previous BRACE work (Kamath et al. 2022) was used to 
# estimate node level abundances

micro_abundances <- merge(x = sample_metadata, y = temp_wet, 
                          by = "Sample_ID") 

write.csv(micro_abundances, "Outputs/micro_abundances.csv")

# Step 6: 
temp0 <- merge(x = temp_wet, y = temp_dry, by = "Sample_ID") %>%
  dplyr::select(Sample_ID, total_micro_ab_standardized, 
                total_meso_ab_standardized) %>%
  mutate(total_ab_standardized = 
           total_micro_ab_standardized + total_meso_ab_standardized)

# Step 7: 
abundance_data <- merge(x = sample_metadata, y = temp0, by = "Sample_ID") %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time)) %>%
  mutate(total_meso_ab_standardized = log(total_meso_ab_standardized))

abundance_data$destructive_time <- as.factor(abundance_data$destructive_time)
abundance_data$moisture_tx <- as.factor(abundance_data$moisture_tx)
abundance_data$temp_tx <- as.factor(abundance_data$temp_tx)
abundance_data$block_effect <- as.factor(abundance_data$block_effect)
abundance_data$tx <- as.factor(abundance_data$tx)

# Step 8:
write.csv(abundance_data, "Outputs/abundances.csv")

# Step 9:
abundance_figs <- abundance_data %>% filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5") # remove outlier

# create plot
meso.ab_plot <- 
  ggplot(abundance_figs, aes(x=interaction(destructive_time, 
                                           temp_tx, moisture_tx), 
                             y=total_meso_ab_standardized,
                             fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create LM
abundance_figs$block_effect <- factor(abundance_figs$block_effect, ordered = F)

meso.ab_lm <- lm(total_meso_ab_standardized ~ (moisture_tx * temp_tx * destructive_time) 
                 + block_effect, data=abundance_figs)

meso.ab_lm2_temp <- lm(total_meso_ab_standardized ~
                         moisture_tx * temp_tx * destructive_time + 
                         moisture_tx * temp_tx +
                         moisture_tx * destructive_time + 
                         temp_tx * destructive_time +
                         moisture_tx +
                         temp_tx +
                         destructive_time +
                         block_effect,
                       data = abundance_figs)

meso.ab_lm2 <- lm(total_meso_ab_standardized ~
                    
                    
                    
                    temp_tx * destructive_time +
                    moisture_tx +
                    temp_tx +
                    destructive_time +
                    block_effect,
                  data = abundance_figs)

summary(meso.ab_lm2)

# create aov
meso.ab_aov <- 
  Anova(meso.ab_lm, type = 'III')

meso.ab_aov2 <- 
  Anova(meso.ab_lm2, type = 'III')

meso.ab_aov2

# summary table
meso.ab_summary <- abundance_data %>%
  dplyr::group_by(interaction(destructive_time, moisture_tx, temp_tx)) %>%
  dplyr::summarize(mean = mean(total_meso_ab_standardized), 
                   se = sd(total_meso_ab_standardized)/
                     sqrt(length(total_meso_ab_standardized)))

# test assumptions
ggqqplot(residuals(meso.ab_lm2))
shapiro_test(residuals(meso.ab_lm2))

plot(meso.ab_lm2, 1)
abundance_figs %>% levene_test(total_meso_ab_standardized ~ 
                                 (moisture_tx * temp_tx * destructive_time))

# call everything
meso.ab_plot; meso.ab_summary; summary(meso.ab_lm); meso.ab_aov