# Trevor Pettit
# July 1, 2023
# T0.R

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
resp_data <- (read.csv("Outputs/respiration_output.csv"))
soil_data <- read.csv("Outputs/soil_outputs.csv")
enzyme_data <- na.omit(read.csv("Data/enzyme data.csv"))
ab_data <- read.csv("Outputs/abundances.csv")
biom_data <- read.csv("Outputs/biom_wide.csv")  %>%
  mutate(total_biomass = mesostigs_biomass + zerconidae_biomass + 
           prostig_astig_biomass + juv_oribatids_biomass + oribatids_biomass +
           collembola_biomass + n.predator_biomass + n.bacterivore_biomass +
           n.fungivore_biomass + n.omnivore_biomass) %>%
  mutate(meso_biomass = mesostigs_biomass + zerconidae_biomass + 
           prostig_astig_biomass + juv_oribatids_biomass + oribatids_biomass +
           collembola_biomass ) %>%
  mutate(micro_biomass = n.predator_biomass + n.bacterivore_biomass +
           n.fungivore_biomass + n.omnivore_biomass) %>%
  dplyr::select(c(Sample_ID, total_biomass, meso_biomass, micro_biomass))
flux_data <- read.csv("Outputs/T0_models.csv")

aggregated_resp_data <- resp_data %>% filter(destructive_time == "T0") 
 
aggregated_resp_data$cum_resp <- rowSums(log(aggregated_resp_data[,7:13]),  na.rm = T)
aggregated_resp_data$cum_resp2 <- rowSums(aggregated_resp_data[,7:13], na.rm = T)

aggregated_resp_data$destructive_time <- 
  as.factor(aggregated_resp_data$destructive_time)
aggregated_resp_data$block_effect <- 
  as.factor(aggregated_resp_data$block_effect)
aggregated_resp_data$temp_tx <- 
  as.factor(aggregated_resp_data$temp_tx)
aggregated_resp_data$moisture_tx <- 
  as.factor(aggregated_resp_data$moisture_tx)

resp_plot <- 
  ggplot(aggregated_resp_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
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
resp_summary <- aggregated_resp_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(cum_resp), 
            se = sd(cum_resp)/sqrt(length(cum_resp))
  )

resp_plot; resp_summary

# create summary table
resp_summary2 <- aggregated_resp_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(cum_resp2), 
            se = sd(cum_resp2)/sqrt(length(cum_resp2))
  )

resp_summary2

aggregated_soil_data <- merge(x = sample_metadata, y = soil_data,
                           by = "Sample_ID") %>%
  filter(destructive_time == "T0") 

aggregated_soil_data$destructive_time <- 
  as.factor(aggregated_soil_data$destructive_time)
aggregated_soil_data$block_effect <- 
  as.factor(aggregated_soil_data$block_effect)
aggregated_soil_data$temp_tx <- 
  as.factor(aggregated_soil_data$temp_tx)
aggregated_soil_data$moisture_tx <- 
  as.factor(aggregated_soil_data$moisture_tx)

moisture_plot <- 
  ggplot(aggregated_soil_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
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
moisture_summary <- aggregated_soil_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(pct_moisture), 
            se = sd(pct_moisture)/sqrt(length(pct_moisture))
  )

moisture_plot; moisture_summary

pH_plot <- 
  ggplot(aggregated_soil_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=pH, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
pH_summary <- aggregated_soil_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(pH), 
            se = sd(pH)/sqrt(length(pH))
  )

pH_summary; pH_plot

CN_plot <- 
  ggplot(aggregated_soil_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
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
CN_summary <- aggregated_soil_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(C_N_ratio), 
            se = sd(C_N_ratio)/sqrt(length(C_N_ratio))
  )

CN_summary; CN_plot

aggregated_enzyme_data <- merge(x = sample_metadata, y = enzyme_data,
                                by = "Sample_ID") %>%
  filter(destructive_time == "T0") 

aggregated_enzyme_data$destructive_time <- 
  as.factor(aggregated_enzyme_data$destructive_time)
aggregated_enzyme_data$block_effect <- 
  as.factor(aggregated_enzyme_data$block_effect)
aggregated_enzyme_data$temp_tx <- 
  as.factor(aggregated_enzyme_data$temp_tx)
aggregated_enzyme_data$moisture_tx <- 
  as.factor(aggregated_enzyme_data$moisture_tx)

# create plot
phenox_plot <- 
  ggplot(aggregated_enzyme_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
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
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(Phenox), 
            se = sd(Phenox)/sqrt(length(Phenox))
  )

# call everything
phenox_plot; phenox_summary

perox_plot <- 
  ggplot(aggregated_enzyme_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
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
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(
    mean = mean(Perox), 
    se = sd(Perox)/sqrt(length(Perox))
  )

perox_plot; perox_summary

aggregated_ab_data <- ab_data %>% filter(destructive_time == "T0")

aggregated_ab_data$destructive_time <- 
  as.factor(aggregated_ab_data$destructive_time)
aggregated_ab_data$block_effect <- 
  as.factor(aggregated_ab_data$block_effect)
aggregated_ab_data$temp_tx <- 
  as.factor(aggregated_ab_data$temp_tx)
aggregated_ab_data$moisture_tx <- 
  as.factor(aggregated_ab_data$moisture_tx)

ab_plot <- 
  ggplot(aggregated_ab_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=total_ab_standardized, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
ab_summary <- aggregated_ab_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(total_ab_standardized), 
            se = sd(total_ab_standardized)/sqrt(length(total_ab_standardized))
  )

ab_plot; ab_summary

meso_ab_plot <- 
  ggplot(aggregated_ab_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=total_meso_ab_standardized, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
meso_ab_summary <- aggregated_ab_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(total_meso_ab_standardized), 
            se = sd(total_meso_ab_standardized)/sqrt(length(total_meso_ab_standardized))
  )

meso_ab_plot; meso_ab_summary

micro_ab_plot <- 
  ggplot(aggregated_ab_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=total_micro_ab_standardized, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
micro_ab_summary <- aggregated_ab_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(total_micro_ab_standardized), 
            se = sd(total_micro_ab_standardized)/sqrt(length(total_micro_ab_standardized))
  )

micro_ab_plot; micro_ab_summary

aggregated_biom_data <- merge(x = sample_metadata, y = biom_data, 
                              by = "Sample_ID") %>%
  filter(destructive_time == "T0")

aggregated_biom_data$destructive_time <- 
  as.factor(aggregated_biom_data$destructive_time)
aggregated_biom_data$block_effect <- 
  as.factor(aggregated_biom_data$block_effect)
aggregated_biom_data$temp_tx <- 
  as.factor(aggregated_biom_data$temp_tx)
aggregated_biom_data$moisture_tx <- 
  as.factor(aggregated_biom_data$moisture_tx)

biom_plot <- 
  ggplot(aggregated_biom_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=total_biomass, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
biom_summary <- aggregated_biom_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(total_biomass), 
            se = sd(total_biomass)/sqrt(length(total_biomass))
  )

biom_plot; biom_summary

meso_biom_plot <- 
  ggplot(aggregated_biom_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=meso_biomass, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
meso_biom_summary <- aggregated_biom_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(meso_biomass), 
            se = sd(meso_biomass)/sqrt(length(meso_biomass))
  )

meso_biom_plot; meso_biom_summary

micro_biom_plot <- 
  ggplot(aggregated_biom_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=micro_biomass, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
micro_biom_summary <- aggregated_biom_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(micro_biomass), 
            se = sd(micro_biomass)/sqrt(length(micro_biomass))
  )

micro_biom_plot; micro_biom_summary

flux_data$destructive_time <- 
  as.factor(flux_data$destructive_time)
flux_data$block_effect <- 
  as.factor(flux_data$block_effect)
flux_data$temp_tx <- 
  as.factor(flux_data$temp_tx)
flux_data$moisture_tx <- 
  as.factor(flux_data$moisture_tx)

flux_plot <- 
  ggplot(flux_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=total_consumption, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
flux_summary <- flux_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(total_consumption), 
            se = sd(total_consumption)/sqrt(length(total_consumption))
  )

flux_plot; flux_summary

Cmin_plot <- 
  ggplot(flux_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=total_Cmin, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
Cmin_summary <- flux_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(total_Cmin), 
            se = sd(total_Cmin)/sqrt(length(total_Cmin))
  )

Cmin_plot; Cmin_summary

Nmin_plot <- 
  ggplot(flux_data, 
         aes(x=interaction(temp_tx, moisture_tx), 
             y=total_Nmin, fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys")

# create summary table
Nmin_summary <- flux_data %>%
  group_by(interaction(destructive_time, temp_tx, moisture_tx)) %>%
  summarise(mean = mean(total_Nmin), 
            se = sd(total_Nmin)/sqrt(length(total_Nmin))
  )

Nmin_plot; Nmin_summary

