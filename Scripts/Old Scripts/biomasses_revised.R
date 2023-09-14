# Trevor Pettit
# May 26, 2023
# biomasses.R

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
abundance_data <- read.csv("Outputs/abundances.csv") 
sample_metadata <- read.csv("Data/sample metadata v2.csv")
meso_abundances <- 
  read.csv("Outputs/meso_abundances.csv") 
micro_abundances <- read.csv("Outputs/micro_abundances.csv")
meso_body_mass <- 
  read.csv("Data/Barreto et al complete dataset .csv") 
# Carlos' body mass data, in ug of C
micro_body_mass <- read.csv("Data/Kamath et al dataset.csv") %>% 
  dplyr::select(-prop_ab) # body mass from Dev's thesis, in ug wwt
micro_prop_ab <- read.csv("Data/Kamath et al dataset.csv") %>% 
  dplyr::select(-body_mass) # proportional abundances from Dev's thesis

# Step 3:
# Barreto et al. body mass data is in ug of C
# biomass = abundance per m2 * body mass (ug of C) * g of C/ug of C 

# Calculate biomass
mesofauna_biomasses <- data.frame("Sample_ID" = abundance_data$Sample_ID) %>%
  mutate("collembola_biomass" =
           meso_abundances$collembola * 
           # abundance per m2
           meso_body_mass$body_mass[meso_body_mass$Node_ID=="collembola"] * 
           # body mass (ug of C)
           1e-6 
         # g of C/ug of C 
         ) %>%  
  mutate("juv_oribatids_biomass" =
           meso_abundances$juv_oribatids *
           meso_body_mass$body_mass[meso_body_mass$Node_ID=="juv_oribatids"] *
           1e-6
         ) %>%
  mutate("oribatids_biomass" =
           meso_abundances$oribatids *
           meso_body_mass$body_mass[meso_body_mass$Node_ID=="oribatids"] *
           1e-6
         ) %>%
  mutate("prostigs_biomass" = 
           meso_abundances$prostigs *
           meso_body_mass$body_mass[meso_body_mass$Node_ID=="prostigs"] *
           1e-6
         ) %>% 
  mutate("astigmata_biomass" = 
           meso_abundances$astigmata * 
           meso_body_mass$body_mass[meso_body_mass$Node_ID=="astigmata"] *
           1e-6
         ) %>%
  mutate("zerconidae_biomass" =
           meso_abundances$zerconidae *
           meso_body_mass$body_mass[meso_body_mass$Node_ID=="zerconidae"] *
           1e-6
         ) %>%
  mutate("mesostigs_biomass" =
           meso_abundances$mesostigs *
           meso_body_mass$body_mass[meso_body_mass$Node_ID=="mesostigmata"] *
           1e-6
         ) %>%
  mutate("total_mesofauna_biomass" = collembola_biomass + 
  juv_oribatids_biomass + oribatids_biomass + prostigs_biomass + 
  astigmata_biomass + zerconidae_biomass + mesostigs_biomass)

# Calculate biomass
# Kamath et al. body mass data is in ug wwt
microfauna_biomasses <- 
  data.frame("Sample_ID" = micro_abundances$Sample_ID) %>%
  # standardized biomass = group proportional abundance * abundance per m2 * 
  # body mass (ug wwt) * ug to g * g wwt to g dwt * g dwt to g of C  
  mutate(n.bacterivore_biomass =
           micro_prop_ab$prop_ab[micro_prop_ab$Node_ID=="Bacterivore"] * 
           # group proportional abundance
           micro_abundances$total_micro_ab_standardized * 
           # abundance per m2
           micro_body_mass$body_mass[micro_body_mass$Node_ID=="Bacterivore"] * 
           # body mass (ug wwt)
           1e-6 * 0.25 * 0.5 
         #ug to g * g wwt to g dwt * g dwt to g of C
         ) %>%
  mutate(n.fungivore_biomass =
           micro_prop_ab$prop_ab[micro_prop_ab$Node_ID=="Fungivore"] *
           micro_abundances$total_micro_ab_standardized *
           1e-6 * 0.25 * 0.5
         ) %>%
  mutate(n.herbivore_biomass =
           micro_prop_ab$prop_ab[micro_prop_ab$Node_ID=="Herbivore"] *
           micro_abundances$total_micro_ab_standardized *
           micro_body_mass$body_mass[micro_body_mass$Node_ID=="Herbivore"] *
           1e-6 * 0.25 * 0.5
         ) %>%
  mutate(n.omnivore_biomass =
           micro_prop_ab$prop_ab[micro_prop_ab$Node_ID=="Omnivore"] *
           micro_abundances$total_micro_ab_standardized *
           micro_body_mass$body_mass[micro_body_mass$Node_ID=="Omnivore"] *
           1e-6 * 0.25 * 0.5
         ) %>%
  mutate(n.predator_biomass =
           micro_prop_ab$prop_ab[micro_prop_ab$Node_ID=="Predator"] *
           micro_abundances$total_micro_ab_standardized *
           micro_body_mass$body_mass[micro_body_mass$Node_ID=="Predator"] *
           1e-6 * 0.25 * 0.5 
         ) %>%
  mutate(total_microfauna_biomass = 
           n.bacterivore_biomass + n.fungivore_biomass + n.omnivore_biomass +
           n.predator_biomass)

# Step 4: 
# Set basal resources, microbial community
labile_biomass <- 135707
recal_biomass <- 135707

set.seed(1)
fungi_biomass <- rnorm(45, mean = 81.21, sd = 0.10*81.21)

set.seed(1)
bacteria_biomass <- rnorm(45, mean = 10.94, sd = 0.10*10.94)

set.seed(1)
protists_biomass <- rnorm(45, mean = 2.59, sd = 0.10*2.59)

# Step 5: 
# Aggregate data
temp <- merge(x = mesofauna_biomasses, y = microfauna_biomasses, 
              by = "Sample_ID") %>%
  mutate(labile_biomass = labile_biomass) %>%
  mutate(recal_biomass = recal_biomass) %>%
  mutate(fungi_biomass = fungi_biomass) %>%
  mutate(bacteria_biomass = bacteria_biomass) %>%
  mutate(protists_biomass = protists_biomass)

# Step 6: 
# Tidy data
temp2 <- temp %>%
  # combine prostig and astig nodes
  mutate(prostig_astig_biomass = prostigs_biomass + astigmata_biomass) %>%
  # remove original cols
  dplyr::select(-prostigs_biomass, -astigmata_biomass)

temp2 <- # re-organize order of nodes, in top-down order
  temp2[, c("Sample_ID", "mesostigs_biomass", "zerconidae_biomass",
            "prostig_astig_biomass", "juv_oribatids_biomass", 
            "oribatids_biomass", "collembola_biomass", 
            "n.predator_biomass", "n.bacterivore_biomass",
            "n.fungivore_biomass", "n.omnivore_biomass", 
            "protists_biomass", "bacteria_biomass", "fungi_biomass", 
            "labile_biomass", "recal_biomass")]

# Step 7: 
# for use in model:
model_biomasses <- as.data.frame(t(temp2[, -(seq(from=0, to=1, by=1))])) 
# transpose data
colnames(model_biomasses) <- temp2$Sample_ID
# fix colnames 

write.csv(model_biomasses, "Outputs/model_biomasses_v2.csv")
write.csv(temp2, "Outputs/biom_wide_v2.csv")

# Step 8:
temp3 <- merge(x = sample_metadata, y = temp2, by = "Sample_ID") 

biomass_data2 <- temp3 %>%
  mutate(total_biomass = mesostigs_biomass + zerconidae_biomass + 
           prostig_astig_biomass + juv_oribatids_biomass + oribatids_biomass +
           collembola_biomass + n.predator_biomass + n.bacterivore_biomass +
           n.fungivore_biomass + n.omnivore_biomass) %>%
  mutate(bact.ch_biomass = zerconidae_biomass + prostig_astig_biomass + 
           collembola_biomass + n.predator_biomass + n.bacterivore_biomass +
           n.omnivore_biomass) %>%
  mutate(fung.ch_biomass = mesostigs_biomass + juv_oribatids_biomass +
           oribatids_biomass + collembola_biomass + n.predator_biomass +
           n.fungivore_biomass + n.omnivore_biomass) %>%
  mutate(meso_biomass = mesostigs_biomass + zerconidae_biomass + 
           prostig_astig_biomass + juv_oribatids_biomass + oribatids_biomass +
           collembola_biomass ) %>%
  mutate(micro_biomass = n.predator_biomass + n.bacterivore_biomass +
           n.fungivore_biomass + n.omnivore_biomass) %>%
  mutate(tx = interaction(temp_tx, moisture_tx, destructive_time)) %>%
  filter(!Sample_ID == "ExRes 5")

biomass_data <- biomass_data2 %>%
  filter(!destructive_time == "T0") 

biomass_data$tx <- droplevels(biomass_data$tx)

## a = T1/Ambient/12C, b = T2/Ambient/12C, c = T1/High/12C,
## d = T2/High/12C, e = T1/Ambient/20C, f = T2/Ambient/20C, 
## g = T1/High/20C, h = T2/High/20C
biomass_data$destructive_time <- as.factor(biomass_data$destructive_time)
biomass_data$moisture_tx <- as.factor(biomass_data$moisture_tx)
biomass_data$temp_tx <- as.factor(biomass_data$temp_tx)
biomass_data$block_effect <- as.factor(biomass_data$block_effect)
biomass_data$tx <- as.factor(biomass_data$tx)

write.csv(temp2, "Outputs/webs_v2.csv")


test <- merge(x = sample_metadata, y = temp2, by = "Sample_ID")

need <- test %>%
  group_by(interaction(destructive_time, moisture_tx, temp_tx)) %>%
  summarize(mean_meso = mean(mesostigs_biomass), 
            se_meso = sd(mesostigs_biomass)/sqrt(length(mesostigs_biomass)),
            mean_zerc = mean(zerconidae_biomass), 
            se_zerc = sd(zerconidae_biomass)/sqrt(length(zerconidae_biomass)),
            mean_miteNP = mean(prostig_astig_biomass), 
            se_miteNP = sd(prostig_astig_biomass)/sqrt(length(prostig_astig_biomass)),
            mean_juvorib = mean(juv_oribatids_biomass), 
            se_juvorib = sd(juv_oribatids_biomass)/sqrt(length(juv_oribatids_biomass)),
            mean_orib = mean(oribatids_biomass), 
            se_orib = sd(oribatids_biomass)/sqrt(length(oribatids_biomass)),
            mean_coll = mean(collembola_biomass), 
            se_coll = sd(collembola_biomass)/sqrt(length(collembola_biomass)),
            mean_prednem = mean(n.predator_biomass), 
            se_prednem = sd(n.predator_biomass)/sqrt(length(n.predator_biomass)),
            mean_bactnem = mean(n.bacterivore_biomass), 
            se_bactnem = sd(n.bacterivore_biomass)/sqrt(length(n.bacterivore_biomass)),
            mean_fungnem = mean(n.fungivore_biomass), 
            se_fungnem = sd(n.fungivore_biomass)/sqrt(length(n.fungivore_biomass)),
            mean_omninem = mean(n.omnivore_biomass), 
            se_omninem = sd(n.omnivore_biomass)/sqrt(length(n.omnivore_biomass)),
            mean_protists = mean(protists_biomass), 
            se_protists = sd(protists_biomass)/sqrt(length(protists_biomass)),
            mean_bact = mean(bacteria_biomass), 
            se_bact = sd(bacteria_biomass)/sqrt(length(bacteria_biomass)),
            mean_fung = mean(fungi_biomass), 
            se_fung = sd(fungi_biomass)/sqrt(length(fungi_biomass)),
            mean_labile = mean(labile_biomass), 
            se_labile = sd(labile_biomass)/sqrt(length(labile_biomass)),
            mean_recal = mean(recal_biomass), 
            se_recal = sd(recal_biomass)/sqrt(length(recal_biomass)),
            )

write.csv(need, "Outputs/mean_node_biomasses.csv")

# Step 9:
# create plot
total.biom_plot <- 
  ggplot(biomass_data, aes(x=interaction(destructive_time, 
                                         temp_tx, moisture_tx), 
                           y=total_biomass, 
                           fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys") 

# Step 10:
# create LM
biomass_data$block_effect <- factor(biomass_data$block_effect, ordered = F)
biomass_data$block_effect <- relevel(biomass_data$block_effect, "low")

total.biom_lm <- lm(total_biomass ~ (moisture_tx * temp_tx * destructive_time) 
                    + block_effect, data=biomass_data)

# create aov
total.biom_aov <- 
  Anova(total.biom_lm, type = 'III')

# Step 11:
# summary table
total.biom_summary <- biomass_data %>%
  dplyr::group_by(interaction(destructive_time, moisture_tx, temp_tx)) %>%
  dplyr::summarize(mean = mean(total_biomass), 
            se = sd(total_biomass)/
              sqrt(length(total_biomass)))

# Step 12:
# call everything
total.biom_plot; total.biom_summary; summary(total.biom_lm); total.biom_aov

# create plot
micro.biom_plot <- 
  ggplot(biomass_data, aes(x=interaction(destructive_time, 
                                         temp_tx, moisture_tx), 
                           y=micro_biomass, 
                           fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys") 

# Step 10:
# create LM
micro.biom_lm <- lm(micro_biomass ~ (moisture_tx * temp_tx * destructive_time) 
                    + block_effect, data=biomass_data)

# create aov
micro.biom_aov <- 
  Anova(micro.biom_lm, type = "III")

# Step 11:
# summary table
micro.biom_summary <- biomass_data %>%
  dplyr::group_by(interaction(destructive_time, moisture_tx, temp_tx)) %>%
  dplyr::summarize(mean = mean(micro_biomass), 
            se = sd(micro_biomass)/
              sqrt(length(micro_biomass)))

# Step 12:
# call everything
micro.biom_plot; micro.biom_summary; summary(micro.biom_lm); micro.biom_aov

# create plot
meso.biom_plot <- 
  ggplot(biomass_data, aes(x=interaction(destructive_time, 
                                         temp_tx, moisture_tx), 
                           y=meso_biomass, 
                           fill = as.factor(destructive_time))) + 
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(aes(pch = as.factor(destructive_time)), 
              position=position_jitter(width=.1, height=0), size = 3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer(palette = "Greys") 

# Step 10:
# create LM
meso.biom_lm <- lm(meso_biomass ~ (moisture_tx * temp_tx * destructive_time) 
                   + block_effect, data=biomass_data)

# create aov
meso.biom_aov <- 
  Anova(meso.biom_lm, type = 'III')

# Step 11:
# summary table
meso.biom_summary <- biomass_data %>%
  dplyr::group_by(interaction(destructive_time, moisture_tx, temp_tx)) %>%
  dplyr::summarize(mean = mean(meso_biomass), 
            se = sd(meso_biomass)/
              sqrt(length(meso_biomass)))

# Step 12:
# call everything
meso.biom_plot; meso.biom_summary; summary(meso.biom_lm); meso.biom_aov
