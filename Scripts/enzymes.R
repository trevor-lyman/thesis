# Trevor Pettit
# June 24, 2023
# enzymes.R

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
enzyme_data <- na.omit(read.csv("Data/enzyme data.csv"))

# tidy data
aggregated_enzyme_data <- merge(x = sample_metadata, y = enzyme_data,
                                by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5")
aggregated_enzyme_data <- aggregated_enzyme_data %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time))
aggregated_enzyme_data$tx <- as.factor(aggregated_enzyme_data$tx)
aggregated_enzyme_data$destructive_time <- 
  as.factor(aggregated_enzyme_data$destructive_time)
aggregated_enzyme_data$block_effect <- 
  as.factor(aggregated_enzyme_data$block_effect)
aggregated_enzyme_data$temp_tx <- 
  as.factor(aggregated_enzyme_data$temp_tx)
aggregated_enzyme_data$moisture_tx <- 
  as.factor(aggregated_enzyme_data$moisture_tx)
aggregated_enzyme_data$tx <- droplevels(aggregated_enzyme_data$tx)
key_enzymes <- levels(aggregated_enzyme_data$tx)

# create plot
phenox_plot <- 
  ggplot(aggregated_enzyme_data, 
         aes(x=interaction(destructive_time, temp_tx, moisture_tx), 
             y=Phenox, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
phenox_summary <- aggregated_enzyme_data %>%
  group_by(tx) %>%
  summarise(mean = mean(Phenox), 
            se = sd(Phenox)/sqrt(length(Phenox))
  ) %>%
  mutate(key = key_enzymes)

# create aov
phenox_aov <- 
  aov(Phenox ~ (moisture_tx * temp_tx * destructive_time) + block_effect, 
      data = aggregated_enzyme_data)

# call everything
phenox_plot; phenox_summary; summary(phenox_aov)

# create plot
perox_plot <- 
  ggplot(aggregated_enzyme_data, 
         aes(x=interaction(destructive_time, temp_tx, moisture_tx), 
             y=Perox, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
perox_summary <- aggregated_enzyme_data %>%
  group_by(tx) %>%
  summarise(
    mean = mean(Perox), 
    se = sd(Perox)/sqrt(length(Perox))
  ) %>%
  mutate(key = key_enzymes)

# create aov
perox_aov <- 
  aov(Perox ~ (moisture_tx * temp_tx * destructive_time) + block_effect, 
      data = aggregated_enzyme_data)

# call everything
perox_plot; perox_summary; summary(perox_aov)

write.csv(aggregated_enzyme_data, "Outputs/enzyme_outputs.csv")

# check assumptions
model_phenox <- lm(Phenox ~ 
                 (moisture_tx * temp_tx * destructive_time) 
               + block_effect, data=aggregated_enzyme_data)

ggqqplot(residuals(model_phenox))
shapiro_test(residuals(model_phenox))

plot(model_phenox, 1)
aggregated_enzyme_data %>% levene_test(Phenox ~ 
                                     (moisture_tx * temp_tx * destructive_time))

model_perox <- lm(Perox ~ 
                     (moisture_tx * temp_tx * destructive_time) 
                   + block_effect, data=aggregated_enzyme_data)

ggqqplot(residuals(model_perox))
shapiro_test(residuals(model_perox))

plot(model_perox, 1)
aggregated_enzyme_data %>% levene_test(Perox ~ 
                                         (moisture_tx * temp_tx * destructive_time))
