library(tidyverse)
library(readxl)
library(bbsBayes2)
library(ggpmisc)
library(rcartocolor)
source("scripts/00_Functions.R")
sci_names_df <- get_sci_names()
dat <- load_bbs_data(level = "state", release = 2023, sample = F, quiet = F)


ests <- list()
species_with_masses <- list()
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
  
}


# Get the species
the_species <- read.csv("data/allspecies.csv") %>% 
  select(-1) %>%
  as_tibble() %>%
  left_join(sci_names_df, by = join_by("species_clean" == "Common Name")) %>%
  distinct() %>%
  mutate(species_bcrForm = if_else(
    species_bcrForm == "Northern House Wren", "House Wren", 
    species_bcrForm))

# Get masses and trophic levels
avonet <- read_excel("data/AVONET2_eBird.xlsx", sheet = 2) %>%
  select(sci_name = Species2, Mass = Mass, trophic_level = Trophic.Level) %>%
  semi_join(the_species, by = join_by(sci_name)) %>%
  left_join(the_species, by = join_by(sci_name)) %>%
  select(species = species_bcrForm, Mass, trophic_level) %>%
  distinct()

df <- do.call(rbind, ests) %>%
  left_join(avonet, by = "species") %>%
  mutate(Nstar = 1920.1 * (Mass ^ -0.81)) %>%
  mutate(distance_from_K = Nstar -  exp(N0_mean)) %>%
  select(species, r_mean, N0_mean, Nstar, distance_from_K, BCR) %>%
  mutate(prop_dist_from_K = exp(N0_mean)/Nstar)

ggplot(df, aes(x = Nstar, y = exp(N0_mean), color = r_mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_point(size = 3, pch = 21, fill = NA, color = 1) +
  stat_ma_line(method = "RMA", range.x = "interval", range.y = "interval") +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(y = "Initial abundance", x = "Allomertic prediction of abundance",
       color = "Change rate") +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  theme_light(14) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5, size = 9),
        legend.key.width = unit(2, "cm"),
        legend.key.height= unit(0.25, "cm"),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())

df %>% 
  mutate(distance_from_K = if_else(distance_from_K > 999, 999, distance_from_K)) %>%
  ggplot(., aes(x = distance_from_K, y = r_mean, color = BCR)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  # geom_text(aes(label = species)) +
  geom_point(size = 3, alpha = 0.5) +
  # geom_point(size = 3, color = "white") +
  # geom_point(size = 3, pch = 21, fill = NA, color = 1, stroke = 0.5) +
  coord_cartesian(xlim = c(-1000,1000), ylim = c(-0.125, 0.125)) +
  scale_x_continuous(labels = c("-1000", "-500", "0", "500", "≥1000"),
                     breaks = c(-1000, -500, 0, 500, 1000)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Distance from carrying capacity (via Damuth 1987)", 
       y = "Average annual rate of population change") +
  theme_light(14) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5, size = 9),
        legend.key.width = unit(2, "cm"),
        legend.key.height= unit(0.25, "cm"),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())

df %>%
  filter(r_mean < 0.11) %>%
  ggplot(aes(x = prop_dist_from_K, y = r_mean)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1) +
  
  geom_point(pch = 21, size = 3, stroke = 0.2, fill = "skyblue") +
  
  # geom_point(aes(color = BCR), size = 3) +
  # geom_point(size = 3, pch = 21, fill = NA, color = "blue1", stroke = 0.1) +
  
  stat_ma_line(method = "RMA", range.y = "interval", range.x = "interval",
               fill = "red3", color = "red3", alpha = 0.3,
               linewidth = 1.25) +
  
  # geom_smooth(method = "lm") +
  
  scale_x_log10(labels = scales::label_comma(), breaks = c(0.0001, 0.001,0.01,0.1,1,10)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_carto_d(palette = 6) +
  labs(x = "Initial abundance as a proportion of K", 
       y = "Average annual rate of population change") +
  theme_light(14) +
  theme(aspect.ratio = 1,
        axis.title = element_text(family = "Arial", face = "bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5, size = 12),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height= unit(0.2, "cm"),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())

df %>%
  filter(r_mean < 0.1) %>%
  ggplot(aes(x = prop_dist_from_K, y = r_mean)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1) +
  geom_point(aes(color = N0_mean), size = 3) +
  # geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp"),
  #             method.args = list(gamma = 1.4), color = "red3", fill = "red") +
  stat_ma_line(method = "RMA", range.x = "interval", range.y = "interval", color = "red3", fill = "red") +
  scale_color_viridis_c(breaks = c(log(0.01),log(0.1),log(1),log(10),log(100),log(1000)),
                        labels = c(0.01, 0.1, 1, 10, 100, 1000)) +
  scale_x_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10,100)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Former abundance as a proportion of\nscaling-derived carrying capacity", 
       y = "Population trend",
       color = "Former abundance") +
  theme_light(14) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        legend.title.position = "bottom",
        legend.title = element_text(hjust = 0.5, size = 12, family = "Roboto"),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height= unit(0.2, "cm"),
        axis.title = element_text(family = "Roboto"),
        axis.text = element_text(family = "Roboto"),
        legend.text = element_text(family = "Roboto"),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())


# ---- Plot the distirbution of declines near K ----

# Get the trends above K
r_abv_K <- df %>%
  filter(r_mean < 0.1) %>%
  filter(prop_dist_from_K >= 0.95) %>%
  pull(r_mean)

# Calculate mean and sd
(K_erode.mean <- mean(r_abv_K))
(K_erode.sd <- sd(r_abv_K))

# Plot the data
hist(r_abv_K)

# Fit a normal density curve
x <- seq(min(r_abv_K), max(r_abv_K), length = 100)
y <- dnorm(x, mean = K_erode.mean, sd = sd(r_abv_K))

# Add the normal density curve
abline(v = 0, col = "red", lwd = 4)
lines(x, y * 0.35, col = "blue", lwd = 3)


