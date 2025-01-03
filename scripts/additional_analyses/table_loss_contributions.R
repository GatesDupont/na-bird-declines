library(tidyverse)
library(readxl)
library(bbsBayes2)
library(ggpmisc)
library(quantreg)
library(rcartocolor)
source("scripts/00_Functions.R")
sci_names_df <- get_sci_names()
dat <- load_bbs_data(level = "state", release = 2023, sample = F, quiet = F)


ests <- list()
species_with_masses <- list()
dynamics <- list()
samps <- list()
for(i in 1:6){
  
  # ---- Load BCR results ----
  
  bcr_num <- c(14, 5, 28, 22, 10, 18)[i]
  
  # Hash table for BCR names and numbers
  bcrs <- data.frame(
    region = c("Atlantic Northern Forest","Northern Pacific Rainforest",
               "Appalachian Mountains", "Eastern Tallgrass Prairie",
               "Northern Rockies", "Shortgrass Prairie"),
    number = c(14, 5, 28, 22, 10, 18)) %>%
    mutate(bcr_code = paste0("BCR", number))
  
  
  # Get bcr name
  bird_conservation_region <- bcrs %>% 
    filter(number == bcr_num) %>% 
    pull(region)
  
  BCR_name_for_file <- gsub(" ", "_", toupper(bird_conservation_region))
  BCR_results_filename <- paste0("data/BCR-results/results_bcr_", BCR_name_for_file, ".rds")
  results <- readRDS(BCR_results_filename)
  
  ests[[i]] <- results$SpeciesEstimates %>%
    mutate(BCR = bird_conservation_region)
  
  samps[[i]] <- results$SpeciesSamples %>%
    mutate(BCR = bird_conservation_region)
  
  
  # ---- Load species body masses ----
  
  avonet <- read_excel("data/AVONET2_eBird.xlsx", sheet = 2) %>%
    select(sci_name = Species2, Mass)
  
  species_with_masses[[i]] <- data.frame(species = ests[[i]]$species) %>%
    mutate(species_clean = str_replace(species, " \\(all forms\\)$", "")) %>%
    mutate(species_clean = if_else(
      species_clean == "House Wren", "Northern House Wren", species_clean)) %>%
    left_join(sci_names_df, by = c(species_clean = "Common Name")) %>%
    as_tibble() %>%
    left_join(avonet, by = join_by(sci_name)) %>%
    left_join(dat$species %>% select(english, aou) %>% distinct(), 
              by = c("species" = "english")) %>%
    select(species, Mass) %>%
    na.omit() %>%
    mutate(BCR = bird_conservation_region)
  
  
  dynamics[[i]] <- results$SpeciesTimeSeries %>%
    mutate(BCR = bird_conservation_region) %>%
    left_join(species_with_masses[[i]], by = join_by(species))
  
}

samps_df <- do.call(rbind, samps)
dynamics_df <- do.call(rbind, dynamics)

dynamics_df %>%
  # Organize the BCR data
  rename(bcr = BCR.x) %>%
  mutate(bcr = factor(bcr, levels = c(
    "Appalachian Mountains", "Atlantic Northern Forest","Eastern Tallgrass Prairie",
    "Northern Rockies", "Northern Pacific Rainforest", "Shortgrass Prairie"))) %>%
  # Calculate the number of species in each region
  group_by(bcr) %>%
  mutate(nspp = n_distinct(species)) %>%
  ungroup() %>%
  # Calculate the abundance lost for each population
  group_by(species, bcr) %>%
  arrange(year) %>%
  summarise(start = first(density_fit),
            end = last(density_fit),
            diff = end-start,
            nspp = unique(nspp)) %>%
  ungroup() %>%
  # Keep only the declining species
  filter(diff < 0) %>% 
  # Rank species by losses and take cumulative sum
  group_by(bcr) %>%
  arrange(diff) %>%
  mutate(loss_rank = row_number()) %>%
  mutate(cs = cumsum(diff),
         csp = cs/sum(diff)) %>%
  ungroup() %>%
  # Rank species by former abundance
  group_by(bcr) %>%
  arrange(desc(start)) %>%
  mutate(init_rank = row_number()) %>%
  ungroup() %>%
  # Keep only the species that contribute to 80% of total losses
  filter(csp <= 0.8) %>%
  group_by(bcr) %>%
  summarise(n = n(),
            nspp = unique(nspp),
            p = n/nspp,
            max_rank = max(init_rank)/nspp) %>%
  ungroup() %>%
  select(bcr, n, p, max_rank)

