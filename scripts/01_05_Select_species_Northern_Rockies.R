library(tidyverse)
library(vroom)
library(bbsBayes2)
library(readxl)
source("scripts/00_Functions.R")
# fetch_bbs_data()
bbs_dat <- load_bbs_data()
sci_names_df <- get_sci_names()


# Hash table for BCR names and numbers
bcrs <- data.frame(
  region = c("Atlantic Northern Forest","Northern Pacific Rainforest",
             "Appalachian Mountains", "Eastern Tallgrass Prairie",
             "Northern Rockies"),
  number = c(14, 5, 28, 22, 10)) %>%
  mutate(bcr_code = paste0("BCR", number))

# Reproducibility
set.seed(1)


# ---- Select BCR region ----

# Select region
bird_conservation_region <- "Northern Rockies"

# Get number ID
bcr_number_id <- bcrs %>% 
  filter(region == bird_conservation_region) %>% 
  pull(number)


# ---- Avonet for bird habitats ----

# Read the Excel file into R
birds_and_habs <- read_excel("data/AVONET2_eBird.xlsx", sheet = 2) %>%
  select(sci_name = Species2, Habitat) %>%
  left_join(sci_names_df, by = "sci_name") %>%
  select(species = `Common Name`, everything()) %>%
  mutate(Habitat = case_when(
    Habitat == "Marine" ~ "Coastal",
    TRUE ~ Habitat))


# ---- Get species with route frequency > 0.05 ----

# Load species AOU codes and english names
sp_names_aou <- bbs_dat %>%
  pluck("species") %>%
  select(aou, english)

# Organize and select species
sp_freqs <- bbs_dat %>%
  pluck("birds") %>%
  # Select BCR of interest
  filter(bcr == bcr_number_id) %>%
  select(route_data_id, aou, year) %>%
  distinct() %>%
  # Calculate number of route-years with that species
  mutate(n_routes = n_distinct(route_data_id)) %>%
  group_by(aou) %>%
  summarise(n = n(), n_routes = unique(n_routes)) %>%
  ungroup() %>%
  # Get species english names
  left_join(sp_names_aou, by = "aou") %>%
  distinct() %>%
  # Remove the species not identified to species level
  filter(!grepl("unid. ", english)) %>%
  # Fix specific species 
  filter(!grepl("Northern Oriole", english)) %>%
  filter(!grepl("Solitary Vireo", english)) %>%
  filter(!grepl("Traill's Flycatcher", english)) %>%
  filter(!grepl("Sapsuckers", english)) %>%
  distinct() %>%
  # Calculate species frequencies
  mutate(prop = n/n_routes) %>%
  arrange(desc(prop)) %>%
  filter(prop >= 0.05) %>%
  mutate(bcr = bcr_number_id)

# These guys have to be in the case_when and in the filter english
sp_freqs %>% filter(grepl("\\(all forms\\)", english))

# Connect birds to habs
selected_birds_w_habs <- sp_freqs %>%
  mutate(species = english) %>%
  # Fixing alternative names (manually)
  mutate(species = case_when(
    species == "(Myrtle Warbler) Yellow-rumped Warbler" ~ "Yellow-rumped Warbler",
    species == "Yellow-rumped Warbler (all forms)" ~ "Yellow-rumped Warbler",
    species == "Northern Flicker (all forms)" ~ "Northern Flicker",
    species == "(Yellow-shafted Flicker) Northern Flicker" ~ "Northern Flicker",
    species == "(Slate-colored Junco) Dark-eyed Junco" ~ "Dark-eyed Junco",
    species == "Dark-eyed Junco (all forms)" ~ "Dark-eyed Junco",
    species == "Red-tailed Hawk (all forms)" ~ "Red-tailed Hawk",
    species == "Yellow-rumped Warbler (all forms)" ~ "Yellow-rumped Warbler",
    TRUE ~ species)) %>%
  filter(!(english %in% c(
    "(Myrtle Warbler) Yellow-rumped Warbler",
    "(Yellow-shafted Flicker) Northern Flicker",
    "(Slate-colored Junco) Dark-eyed Junco",
    "(Oregon Junco) Dark-eyed Junco",
    "Western Flycatcher (Cordilleran/Pacific-slope)",
    "Pacific-slope Flycatcher",
    "(Red-shafted Flicker) Northern Flicker",
    "(Audubon's Warbler) Yellow-rumped Warbler",
    "Blue Grouse (Dusky/Sooty)",
    "Dark-eyed Junco",
    "Red-tailed Hawk"))) %>%
  filter(!is.na(english)) %>%
  # Remove coastal species
  left_join(birds_and_habs, by = "species")

# Check that no birds have NA habs (this should be empty)
selected_birds_w_habs %>% filter(is.na(Habitat))

# Remove coastal birds
selected_birds_from_good_habs <- selected_birds_w_habs %>%
  filter(Habitat != "Coastal")


# ---- Save species selections ----

# Grab species
selected_species <- selected_birds_from_good_habs %>% select(species = english, bcr)

# Save species as csv
BCR_file_suffix <- gsub(" ", "_", toupper(bird_conservation_region))
BCR_file <- paste("data/BCR-species/species_bcr_", BCR_file_suffix, ".csv", sep = "")
vroom_write(selected_species, BCR_file)

