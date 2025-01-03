library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(cowplot)
library(rcartocolor)
library(lubridate)
library(readxl)
library(vroom)
library(bbsBayes2)
library(mgcv)
library(gratia)
library(bbmle)
source("scripts/00_Functions.R")
sci_names_df <- get_sci_names()
dat <- load_bbs_data(level = "state", release = 2023, sample = F, quiet = F)
# Define the function for density calculation
count_to_density <- function(N, C_D, C_P) {
  N / (50 * pi * 400^2) * ((400 / C_D)^2) * C_P * (1.5e7 / 1)
}


p1 <- list()
p2 <- list()
p3 <- list()
chg_dat <- list()
for(i in 1:6){
  
  bcr_num <- c(14, 5, 28, 22, 10, 18)[i]
  
  # ---- Get BCR name ----
  
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
  
  cat(i, " ", bird_conservation_region, "\n")
  
  
  # ---- Species ----
  
  detectability_factors <- read_xlsx("data/detectability-factors/detectability_factors.xlsx")
  
  # File name
  BCR_file_suffix <- gsub(" ", "_", toupper(bird_conservation_region))
  BCR_file <- paste("data/BCR-species/species_bcr_", BCR_file_suffix, ".csv", sep = "")
  
  spp_initial <- vroom(BCR_file, show_col_types=F) %>% select(species)
  
  sp_to_rm <- spp_initial %>% 
    mutate(species = gsub(" \\(all forms\\)", "", species)) %>%
    left_join(detectability_factors, by = "species") %>%
    filter(is.na(C_D)) %>%
    pull(species) %>%
    c(., "Eurasian Collared-Dove")
  
  spp <- spp_initial %>% 
    filter(!(species %in% sp_to_rm)) %>% 
    pull(species)
  
  
  # ---- Load data ----
  
  aou_nums <- sapply(spp, FUN = function(x){check_species(
    species = x, 
    species_list = dat$species, 
    combine_species_forms = T, 
    quiet = T)}) %>%
    as.numeric()
  
  
  # ---- Load species body masses ----
  
  avonet <- read_excel("data/AVONET2_eBird.xlsx", sheet = 2) %>%
    select(sci_name = Species2, Mass)
  
  species_with_masses <- data.frame(species = spp) %>%
    mutate(species_clean = str_replace(species, " \\(all forms\\)$", "")) %>%
    mutate(species_clean = if_else(
      species_clean == "House Wren", "Northern House Wren", species_clean)) %>%
    left_join(sci_names_df, by = c(species_clean = "Common Name")) %>%
    as_tibble() %>%
    left_join(avonet, by = join_by(sci_name)) %>%
    left_join(dat$species %>% select(english, aou) %>% distinct(), 
              by = c("species" = "english")) %>%
    select(aou, Mass) %>%
    na.omit()
  
  
  # ---- Organize data ----
  
  routes <- dat$routes %>%
    filter(bcr == bcr_num)
  
  birds <- dat$birds %>%
    # Filter to the BCR
    filter(bcr == bcr_num) %>%
    # Keep only the species in our study
    filter(aou %in% aou_nums) %>%
    # Join the route data which includes observation data
    select(country_num, state_num, route, year, aou, species_total) %>%
    left_join(routes, by = join_by(country_num, state_num, route, year)) %>%
    # Format the time
    filter(!is.na(start_time)) %>%
    filter(!is.na(end_time)) %>%
    mutate(start_time = sprintf("%04d", start_time),
           start_time = str_replace(start_time, "^(..)(..)$", "\\1:\\2"),
           start_time = hm(start_time)) %>%
    mutate(end_time = sprintf("%04d", end_time),
           end_time = str_replace(end_time, "^(..)(..)$", "\\1:\\2"),
           end_time = hm(end_time)) %>%
    mutate(duration = as.numeric(end_time - start_time)/(60*60),
           duration = scale(duration)[,1]) %>%
    # Format the date to julian dat
    mutate(date = make_date(year, month, day),
           yday = yday(date)) %>%
    mutate(yday = scale(yday)[,1]) %>%
    # Identify the observer ID
    rename(observer = obs_n) %>%
    # Create the route ID
    mutate(route = str_c(country_num, "-", state_num, "-", route)) %>%
    select(route, observer, year, aou, species_total, yday, duration) %>%
    # Join the mass and detectability data
    left_join(species_with_masses, by = "aou") %>%
    left_join(detectability_factors, by = "aou") %>%
    # Convert count to denstiy
    mutate(density = count_to_density(species_total, C_D, C_P)) %>%
    # Get biomass from mass x density
    mutate(biomass = density * Mass) %>%
    # Remove missing observations
    na.omit() %>%
    # Make observer first year effect
    group_by(observer) %>%
    mutate(first_year_effect = if_else(year == min(year), 1, 0)) %>%
    ungroup() %>%
    # Aggregate to get route level density and biomass
    group_by(route, observer, year) %>%
    summarise(bird_density = sum(density), 
              bird_biomass = sum(biomass),
              duration = unique(duration),
              yday = unique(yday),
              first_year_effect = unique(first_year_effect)) %>%
    ungroup() %>%
    # Format the random effects as factors
    mutate(observer = as.factor(observer),
           route = as.factor(route)) %>%
    # Start the years at 0
    mutate(min_year = min(year),
           year = year - min_year)
  
  
  # ---- Fit the model ----
  
  m_density <- bam(
    bird_density ~ s(year, k = 6) + 
      s(yday, k = 5) + 
      s(duration) + 
      first_year_effect +
      s(observer, bs = "re") + 
      s(route, bs = "re"),
    data = birds,
    family = gaussian(link = "log"),
    control = list(trace = F),
    discrete = T,
    gamma = 1.4)
  
  m_biomass <- bam(
    bird_biomass ~ s(year, k = 6) + 
      s(yday, k = 5) + 
      s(duration) + 
      first_year_effect +
      s(observer, bs = "re") + 
      s(route, bs = "re"),
    data = birds,
    family = gaussian(link = "log"),
    control = list(trace = F),
    discrete = T,
    gamma = 1.4)
  
  
  # ---- Make predictions data frame ----
  
  predictions_df <- birds %>%
    select(-bird_density, -bird_biomass) %>%
    mutate(
      yday = 0,
      duration = 0,
      first_year_effect = 0,
      observer = unique(observer)[1],
      route = unique(route)[1]) %>%
    distinct() %>%
    arrange(year)
  
  nuisance_covariates = c(
    "s(yday)", 
    "s(duration)",
    "first_year_effect",
    "s(observer)",
    "s(route)")
  
  
  # ---- Plot density predictions ----
  
  predictions_density <- predictions_df %>%
    fitted_values(m_density, data = ., exclude = nuisance_covariates) %>%
    mutate(init_y = if_else(year == min(year), .fitted, NA)) %>%
    mutate(init_y = max(init_y, na.rm = T)) %>%
    mutate(diff = (.fitted - init_y)/init_y) %>%
    mutate(diff_lwr = (.lower_ci - init_y)/init_y) %>%
    mutate(diff_upr = (.upper_ci - init_y)/init_y) %>%
    select(year, diff, diff_lwr, diff_upr) %>%
    mutate(variable = "Density")
  
  
  # ---- Plot biomass predictions ----
  
  predictions_biomass <- predictions_df %>%
    fitted_values(m_biomass, data = ., exclude = nuisance_covariates) %>%
    mutate(init_y = if_else(year == min(year), .fitted, NA)) %>%
    mutate(init_y = max(init_y, na.rm = T)) %>%
    mutate(diff = (.fitted - init_y)/init_y) %>%
    mutate(diff_lwr = (.lower_ci - init_y)/init_y) %>%
    mutate(diff_upr = (.upper_ci - init_y)/init_y) %>%
    select(year, diff, diff_lwr, diff_upr) %>%
    mutate(variable = "Biomass")
  
  chg_dat[[i]] <- rbind(predictions_density, 
                    predictions_biomass) %>%
    mutate(bcr = bird_conservation_region)
  
  
  # ---- Plot together ----
  
  p2[[i]] <- rbind(predictions_density, predictions_biomass) %>%
    ggplot(aes(x = year + unique(birds$min_year), y = diff, group = variable, color = variable, fill = variable)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4, color = "white") +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4) +
    geom_ribbon(aes(ymin = diff_lwr, ymax = diff_upr), alpha = 0.4, color = NA) +
    geom_line(linewidth = 1.25) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(breaks = scales::pretty_breaks(), 
                       labels = scales::percent_format(),
                       limits = c(-0.45, 0.1)) +
    scale_color_manual(values = c("orange2","coral2")) +
    scale_fill_manual(values = c("orange2","coral2")) +
    theme_light(16) +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          aspect.ratio = 1) 
  
  # ---- Plot RAD ----
  
  BCR_name_for_file <- gsub(" ", "_", toupper(bird_conservation_region))
  BCR_results_filename <- paste0("data/BCR-results/results_bcr_", BCR_name_for_file, ".rds")
  results <- readRDS(BCR_results_filename)
  species_preds_df <- results$SpeciesTimeSeries %>%
    group_by(year) %>%
    arrange(desc(density_fit)) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    select(species, year, density_fit, rank)
  
  p1[[i]] <- ggplot(species_preds_df, aes(x = rank, y = density_fit, group = year, color = year)) +
    geom_line(linewidth = 1) +
    scale_x_log10(breaks = c(1,10,100)) +
    scale_color_carto_c(palette = "ag_GrnYl") +
    labs(x = NULL, y = NULL) +
    theme_light(16) +
    theme(panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          aspect.ratio = 1) 
  
  cat("done \n\n")
  
}

# library(cowplot)
# 
# ptop <- plot_grid(p1[[2]] + theme(axis.text.x = element_text(size = 11)),
#                   p1[[5]] + theme(axis.text.x = element_text(size = 11)), 
#                   p1[[6]] + theme(axis.text.x = element_text(size = 11)),
#                   p1[[4]] + theme(axis.text.x = element_text(size = 11)),
#                   p1[[3]] + theme(axis.text.x = element_text(size = 11)), 
#                   p1[[1]] + theme(axis.text.x = element_text(size = 11)), 
#                   nrow = 1, align = "hv")
# 
# pbottom <- plot_grid(p2[[2]] + theme(axis.text.x = element_text(size = 9)),
#                      p2[[5]] + theme(axis.text.x = element_text(size = 9)), 
#                      p2[[6]] + theme(axis.text.x = element_text(size = 9)),
#                      p2[[4]] + theme(axis.text.x = element_text(size = 9)),
#                      p2[[3]] + theme(axis.text.x = element_text(size = 9)), 
#                      p2[[1]] + theme(axis.text.x = element_text(size = 9)), 
#                      nrow = 1, align = "hv")
# 
# plot_grid(ptop, pbottom, nrow = 2, align = "hv")
# 
# 
# ptop <- plot_grid(p1[[2]] +  labs(x = "Rank") + theme(axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 11)),
#                   p1[[5]] +  labs(x = "Rank") + theme(axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 11)), 
#                   p1[[6]] +  labs(x = "Rank") + theme(axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 11)),
#                   nrow = 1, align = "hv")
# 
# pbottom <- plot_grid(p2[[2]] +  labs(x = "Year") + theme(axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 9)),
#                      p2[[5]] +  labs(x = "Year") + theme(axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 9)), 
#                      p2[[6]] +  labs(x = "Year") + theme(axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 9)),
#                      nrow = 1, align = "hv")
# 
# plot_grid(ptop, pbottom, nrow = 2, align = "hv")
# 
# 
# 
# (p1[[3]] + p1[[1]] + p1[[4]]) / (p1[[5]] + p1[[2]] + p1[[6]])
# 
# 
# (p2[[2]] + p2[[2]] + p2[[1]])/(p2[[2]] + p2[[2]] + p2[[3]])/(p2[[4]] + p2[[5]] + p2[[6]])
# p2[[2]] + geom_line(linewidth = 2) + theme(axis.text = element_text(size = 18))
# 
# 
# bcrs$number
# 
# bcrs$number[2]; bcrs$number[2]; bcrs$number[1]
# bcrs$number[2]; bcrs$number[2]; bcrs$number[3]
# bcrs$number[4]; bcrs$number[5]; bcrs$number[6]

chg_dat %>% 
  do.call(rbind,.) %>%
  mutate(bcr = factor(bcr, levels =  c(
    "Appalachian Mountains", "Atlantic Northern Forest", "Eastern Tallgrass Prairie",
    "Northern Rockies", "Northern Pacific Rainforest", "Shortgrass Prairie"))) %>%
  group_by(bcr) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  filter(variable == "Biomass") %>%
  select(-year) %>%
  mutate(across(is.numeric, function(x) 100 * round(x,2))) %>%
  arrange(bcr)
  






