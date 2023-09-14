# Trevor Pettit
# Sept 14, 2023
# soil_data_revised_v2.R

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
pH_data <- read.csv("Data/pH data.csv")
soil_moisture_data <- read.csv("Data/soil moisture data.csv")
CNS_data <- read.csv("Data/CNS data.csv")

sample_metadata$block_effect <- ordered(sample_metadata$block_effect, levels = c("low", "high"))

# tidy data
aggregated_pH_data <- merge(x = sample_metadata, y = pH_data,
                            by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5")
aggregated_pH_data <- aggregated_pH_data %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time))
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
                                          y=pH, 
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
  summarise(mean = mean(pH), se = sd(pH)/sqrt(length(pH))) %>%
  mutate(key = key_pH)

# create LM
aggregated_pH_data$block_effect <- factor(aggregated_pH_data$block_effect, ordered = F)
pH_lm <- lm(pH ~ (moisture_tx * temp_tx * destructive_time) + block_effect, 
            data = aggregated_pH_data)

pH_lm2_temp <- lm(pH ~ moisture_tx * temp_tx * destructive_time + 
                  moisture_tx * temp_tx +
                  moisture_tx * destructive_time + 
                  temp_tx * destructive_time +
                  moisture_tx +
                  temp_tx +
                  destructive_time +
                  block_effect,
                data = aggregated_pH_data)

pH_lm2 <- lm(pH ~ 
                  

                  
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
aggregated_pH_data %>% levene_test(pH ~ 
                                     (moisture_tx * temp_tx * destructive_time))

# call everything
pH_plot; pH_summary; summary(pH_lm); pH_aov

# tidy data
aggregated_soil_moisture_data <- merge(x = sample_metadata, 
                                       y = soil_moisture_data,
                                       by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5")
aggregated_soil_moisture_data <- aggregated_soil_moisture_data %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time))
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

# tidy data
aggregated_CNS_data <- merge(x = sample_metadata, y = CNS_data,
                             by = "Sample_ID") %>%
  filter(!destructive_time == "T0") %>%
  filter(!Sample_ID == "ExRes 5")
aggregated_CNS_data <- aggregated_CNS_data %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time))
aggregated_CNS_data$destructive_time <- 
  as.factor(aggregated_CNS_data$destructive_time)
aggregated_CNS_data$block_effect <- 
  as.factor(aggregated_CNS_data$block_effect)
aggregated_CNS_data$moisture_tx <- 
  as.factor(aggregated_CNS_data$moisture_tx)
aggregated_CNS_data$temp_tx <- 
  as.factor(aggregated_CNS_data$temp_tx)
aggregated_CNS_data$tx <- as.factor(aggregated_CNS_data$tx)
aggregated_CNS_data$tx <- droplevels(aggregated_CNS_data$tx)
key_CN <- levels(aggregated_CNS_data$tx)

# create plot
CN_plot <- 
  ggplot(aggregated_CNS_data, 
         aes(x=interaction(destructive_time, temp_tx, moisture_tx), 
             y=C_N_ratio, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_fill_brewer(palette = "Greys")

# create summary table
CN_summary <- aggregated_CNS_data %>%
  group_by(tx) %>%
  summarise(mean = mean(C_N_ratio), 
            se = sd(C_N_ratio)/sqrt(length(C_N_ratio))) %>%
  mutate(key = key_CN)

# create LM
aggregated_CNS_data$block_effect <- factor(aggregated_CNS_data$block_effect, ordered = F)

CN_lm <- lm(C_N_ratio ~ (moisture_tx * temp_tx * destructive_time) + block_effect, 
            data = aggregated_CNS_data)

CN_lm2_temp <- lm(C_N_ratio ~
                    moisture_tx * temp_tx * destructive_time + 
                    moisture_tx * temp_tx +
                    moisture_tx * destructive_time + 
                    temp_tx * destructive_time +
                    moisture_tx +
                    temp_tx +
                    destructive_time +
                    block_effect,
               data = aggregated_CNS_data)

CN_lm2 <- lm(C_N_ratio ~
              
              
              
              
              moisture_tx +
              temp_tx +
              destructive_time +
              block_effect,
          data = aggregated_CNS_data)

summary(CN_lm2)

# create aov
CN_aov <- 
  Anova(CN_lm, type = 'III')

CN_aov2 <- 
  Anova(CN_lm2, type = 'III')

CN_aov2

# test assumptions
ggqqplot(residuals(CN_lm2))
shapiro_test(residuals(CN_lm2))

plot(CN_lm2, 1)
aggregated_CNS_data %>% levene_test(C_N_ratio ~ 
                                                (moisture_tx * temp_tx * destructive_time))

# call everything
CN_plot; CN_summary; summary(CN_lm); CN_aov

# aggregate everything
soil_data_temp <- pH_data %>%
  select(!pH_weight) %>%
  mutate("pct_moisture" = soil_moisture_data$pct_moisture) 

soil_data <- merge(x = soil_data_temp, y = CNS_data, 
                   by = "Sample_ID", all = T) %>%
  select(Sample_ID, pH, pct_moisture, C_N_ratio)

write.csv(soil_data, "Outputs/soil_outputs.csv")
