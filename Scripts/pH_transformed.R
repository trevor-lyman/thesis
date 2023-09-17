# Trevor Pettit
# Sept 14, 2023
# pH_transformed.R

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
aggregated_pH_data <- merge(x = sample_metadata, y = pH_data,
                            by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5") %>%
  filter(!Sample_ID == "ExRes 17") %>% # NEW
  filter(!Sample_ID == "ExRes 25") # NEW

aggregated_pH_data <- aggregated_pH_data %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time)) %>%
  mutate(pH = log(pH))
aggregated_pH_data$tx <- as.factor(aggregated_pH_data$tx)
aggregated_pH_data$destructive_time <- 
  as.factor(aggregated_pH_data$destructive_time)
aggregated_pH_data$moisture_tx <- as.factor(aggregated_pH_data$moisture_tx)
aggregated_pH_data$temp_tx <- as.factor(aggregated_pH_data$temp_tx)
aggregated_pH_data$block_effect <- as.factor(aggregated_pH_data$block_effect)
aggregated_pH_data$tx <- droplevels(aggregated_pH_data$tx)
key_pH <- levels(aggregated_pH_data$tx)

# create plot
pH_plot <- ggplot(aggregated_pH_data, aes(x=interaction(destructive_time, 
                                                        temp_tx, moisture_tx), 
                                          y=log(pH), 
                                          fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
pH_summary <- aggregated_pH_data %>%
  group_by(tx) %>%
  summarise(mean = mean((pH)), se = sd(pH)/sqrt(length((pH)))) %>%
  mutate(key = key_pH)

# create LM
aggregated_pH_data$block_effect <- factor(aggregated_pH_data$block_effect, ordered = F)
pH_lm <- lm((pH) ~ (moisture_tx * temp_tx * destructive_time) + block_effect, 
            data = aggregated_pH_data)

pH_lm2_temp <- lm((pH) ~ moisture_tx * temp_tx * destructive_time + 
                    moisture_tx * temp_tx +
                    moisture_tx * destructive_time + 
                    temp_tx * destructive_time +
                    moisture_tx +
                    temp_tx +
                    destructive_time +
                    block_effect,
                  data = aggregated_pH_data)

pH_lm2 <- lm((pH) ~ 
               
               
               
               moisture_tx +
               temp_tx +
               destructive_time +
               block_effect,
             data = aggregated_pH_data)

summary(pH_lm2)

# create aov
pH_aov <- Anova(pH_lm, type = 'III')
pH_aov2 <- Anova(pH_lm2, type = 'III')

pH_aov2

# test assumptions
ggqqplot(residuals(pH_lm2))
shapiro_test(residuals(pH_lm2))

plot(pH_lm2, 1)
aggregated_pH_data %>% levene_test((pH) ~ 
                                     (moisture_tx * temp_tx * destructive_time))

# call everything
pH_plot; pH_summary; summary(pH_lm); pH_aov

temp <- lm(pH~moisture_tx*temp_tx*destructive_time+block_effect, data = aggregated_pH_data)
tukey <- aov(temp)
TukeyHSD(tukey)

