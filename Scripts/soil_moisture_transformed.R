# Trevor Pettit
# Sept 14, 2023
# soil_moisture_transformed.R

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

sample_metadata <- read.csv("Data/sample metadata v2.csv")
pH_data <- read.csv("Data/pH data.csv")
soil_moisture_data <- read.csv("Data/soil moisture data.csv")
CNS_data <- read.csv("Data/CNS data.csv")

sample_metadata$block_effect <- ordered(sample_metadata$block_effect, levels = c("low", "high"))

# tidy data
aggregated_soil_moisture_data <- merge(x = sample_metadata, 
                                       y = soil_moisture_data,
                                       by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5") %>%
  filter(!Sample_ID == "ExRes 17") %>% # NEW
  filter(!Sample_ID == "ExRes 27") # NEW
aggregated_soil_moisture_data <- aggregated_soil_moisture_data %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time)) %>%
  mutate(pct_moisture = log(pct_moisture))
aggregated_soil_moisture_data$destructive_time <- 
  as.factor(aggregated_soil_moisture_data$destructive_time)
aggregated_soil_moisture_data$block_effect <- 
  as.factor(aggregated_soil_moisture_data$block_effect)
aggregated_soil_moisture_data$moisture_tx <- 
  as.factor(aggregated_soil_moisture_data$moisture_tx)
aggregated_soil_moisture_data$temp_tx <- 
  as.factor(aggregated_soil_moisture_data$temp_tx)
aggregated_soil_moisture_data$tx <- 
  as.factor(aggregated_soil_moisture_data$tx)
aggregated_soil_moisture_data$tx <- 
  droplevels(aggregated_soil_moisture_data$tx)
key_soil_moisture <- levels(aggregated_soil_moisture_data$tx)

# create plot
soil_moisture_plot <- 
  ggplot(aggregated_soil_moisture_data, 
         aes(x=interaction(destructive_time, temp_tx, moisture_tx), 
             y=pct_moisture, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
soil_moisture_summary <- aggregated_soil_moisture_data %>%
  group_by(tx) %>%
  summarise(mean = mean(pct_moisture), 
            se = sd(pct_moisture)/sqrt(length(pct_moisture))) %>%
  mutate(key = key_soil_moisture)

# create LM
aggregated_soil_moisture_data$block_effect <- factor(aggregated_soil_moisture_data$block_effect, ordered = F)
soil_moisture_lm <- lm(pct_moisture ~ (moisture_tx * temp_tx * destructive_time) + block_effect, 
                       data = aggregated_soil_moisture_data)

soil_moisture_lm2_temp <- lm(pct_moisture ~
                               moisture_tx * temp_tx * destructive_time + 
                               moisture_tx * temp_tx +
                               moisture_tx * destructive_time + 
                               temp_tx * destructive_time +
                               moisture_tx +
                               temp_tx +
                               destructive_time +
                               block_effect,
                             data = aggregated_soil_moisture_data)

soil_moisture_lm2 <- lm(pct_moisture ~
                          
                          
                          
                          
                          moisture_tx +
                          temp_tx +
                          destructive_time +
                          block_effect,
                        data = aggregated_soil_moisture_data)

summary(soil_moisture_lm2)

# create aov
soil_moisture_aov <- 
  Anova(soil_moisture_lm, type = 'III')

soil_moisture_aov2 <- 
  Anova(soil_moisture_lm2, type = 'III')

soil_moisture_aov2

# test assumptions
ggqqplot(residuals(soil_moisture_lm2))
shapiro_test(residuals(soil_moisture_lm2))

plot(soil_moisture_lm2, 1)
aggregated_soil_moisture_data %>% levene_test(pct_moisture ~ 
                                                (moisture_tx * temp_tx * destructive_time))

# call everything
soil_moisture_plot; soil_moisture_summary; summary(soil_moisture_lm); soil_moisture_aov

temp <- lm(pct_moisture~moisture_tx*temp_tx*destructive_time+block_effect, data = aggregated_soil_moisture_data)
tukey <- aov(temp)
TukeyHSD(tukey)

