# Tidyverse and ggpmisc
library(tidyverse)
suppressWarnings(library(ggpmisc)) 
# Data Import and Handling
library(vroom) 
library(readxl) 
# Statistical Modeling
library(mgcv)
library(gratia)
library(lmodel2) 
library(car)
# Progress Monitoring
library(progress)
# BBS data
library(bbsBayes2)
# Custom functions & data
source("scripts/00_Functions.R")
source("scripts/fit_model_with_fallback.R")
source("scripts/species_abundance_estimates.R")
sci_names_df <- get_sci_names()


# --- Repro ----

set.seed(1)


# ---- Setup ----

# Load all BBS data
dat <- load_bbs_data(level = "state", release = 2023, sample = F, quiet = F)

# Initial points
min_year <- 1966
detectability_factors <- read_xlsx("data/detectability-factors/detectability_factors.xlsx")

# Hash table for BCR names and numbers
bcrs <- data.frame(
  region = c("Atlantic Northern Forest","Northern Pacific Rainforest",
             "Appalachian Mountains", "Eastern Tallgrass Prairie",
             "Northern Rockies", "Shortgrass Prairie"),
  number = c(14, 5, 28, 22, 10, 18)) %>%
  mutate(bcr_code = paste0("BCR", number))


# ---- Select BCR region ----

# bird_conservation_region <- "Atlantic Northern Forest"
# bird_conservation_region <- "Northern Pacific Rainforest"
# bird_conservation_region <- "Appalachian Mountains"
bird_conservation_region <- "Eastern Tallgrass Prairie"
# bird_conservation_region <- "Northern Rockies"
# bird_conservation_region <- "Shortgrass Prairie"


# Get number ID
bcr_code <- bcrs %>% 
  filter(region == bird_conservation_region) %>% 
  pull(bcr_code)

# Identify only the routes in the BCR
routes_bcr <- dat$routes %>%
  filter(bcr == bcrs %>% 
           filter(region == bird_conservation_region) %>%
           pull(number))

# # Get time zones
# tzs <- routes_bcr %>%
#   select(latitude, longitude) %>%
#   distinct() %>%
#   mutate(tz = tz_lookup_coords(lat = latitude, lon = longitude,  method = "accurate"))


# ---- Get species ----

# File name
BCR_file_suffix <- gsub(" ", "_", toupper(bird_conservation_region))
BCR_file <- paste("data/BCR-species/species_bcr_", BCR_file_suffix, ".csv", sep = "")

spp_initial <- vroom(BCR_file) %>% select(species)

sp_to_rm <- spp_initial %>% 
  mutate(species = gsub(" \\(all forms\\)", "", species)) %>%
  left_join(detectability_factors, by = "species") %>%
  filter(is.na(C_D)) %>%
  pull(species) %>%
  c(., "Eurasian Collared-Dove")

spp <- spp_initial %>% 
  filter(!(species %in% sp_to_rm)) %>% 
  pull(species)

n_spp <- length(spp)


# ---- AVONET for masses ----

avonet <- read_excel("data/AVONET2_eBird.xlsx", sheet = 2) %>%
  select(sci_name = Species2, Mass)

species_with_masses <- data.frame(species = spp) %>%
  mutate(species_clean = str_replace(species, " \\(all forms\\)$", "")) %>%
  mutate(species_clean = if_else(
    species_clean == "House Wren", "Northern House Wren", species_clean)) %>%
  left_join(sci_names_df, by = c(species_clean = "Common Name")) %>%
  as_tibble() %>%
  left_join(avonet)


# ---- Loop through species ----

pb <- progress_bar$new(
  format = "[:bar] :percent elapsed: :elapsed eta: :eta",
  total = n_spp, clear = FALSE)

options(warn = 1)
warnings_by_iteration <- list()

species_samples <- list()
species_estimates <- list()
species_parameter_variances <- list()
species_preds <- list()
species_preds_sims <- list()
species_FirstLastDensities <- list()
for(i in 1:n_spp){
  
  # ---- Fetch the BBS data ----
  
  # Get the aou number for the species
  sp_aou_num <- check_species(
    species = spp[i], 
    species_list = dat$species, 
    combine_species_forms = T, 
    quiet = T)
  
  # Keep species i records only in the BCR
  sp_i_inBCR <- dat$birds %>%
    filter(aou == sp_aou_num) %>%
    semi_join(routes_bcr, by = c("country_num", "state_num", "route")) %>%
    left_join(routes_bcr, by = c("country_num", "state_num", "route", "year")) %>%
    mutate(observer = obs_n) %>%
    select(country_num, state_num, route, 
           species_total,
           route, observer,
           year, month, day, 
           start_time, end_time)
  
  # Get route-years in the BCR that don't have species i
  zeros_routes <- routes_bcr %>%
    anti_join(sp_i_inBCR, by = c("country_num", "state_num", "route", "year")) %>%
    mutate(observer = obs_n) %>%
    mutate(species_total = 0) %>%
    select(country_num, state_num, route, 
           species_total,
           route, observer,
           year, month, day, 
           start_time, end_time)
  
  # Combined data for all of BCR
  df <- sp_i_inBCR %>%
    # Combine zero and non-zero surveys
    rbind(., zeros_routes) %>%
    # Convert route and observer to factors for modeling
    mutate(route = str_c(country_num, state_num, route, sep = "-")) %>%
    mutate(route = as.factor(route),
           observer = as.factor(observer)) %>%
    # Make a first year effect
    group_by(observer) %>%
    mutate(first_year = min(year)) %>%
    ungroup() %>%
    mutate(first_year = if_else(year == first_year, 1, 0)) %>%
    # Convert times to proper format to calculate duration
    mutate(start_time = sprintf("%04d", start_time),
           start_time = strptime(start_time, format = "%H%M"),
           start_time = hour(start_time) * 60 + minute(start_time)) %>%
    mutate(end_time = sprintf("%04d", end_time),
           end_time = strptime(end_time, format = "%H%M"),
           end_time = hour(end_time) * 60 + minute(end_time)) %>%
    # Date to yday
    mutate(date = make_date(year, month, day),
           yday = yday(date)) %>%
    mutate(yday = scale(yday)[,1]) %>%
    # Arranging how I like the columns ordered
    select(observer, first_year, route,
           year, yday,
           start_time, end_time,
           species_total) %>%
    # Making a random effect for year
    mutate(year_factor = as.factor(as.character(year))) %>%
    mutate(duration = scale(end_time - start_time)[,1]) %>%
    mutate(year_full = year) %>%
    mutate(year = year - min_year) %>%
    na.omit() 
  
  
  # ---- Modelling ----
  
  m <- fit_model_with_fallback(df)
  
  
  # ---- Calculate N0 and r ----
  
  # Get detectability-adjusted abundance
  det_facts_sp <- detectability_factors %>%
    filter(aou == sp_aou_num)
  
  # Assign det adj facts to variables
  C_P <- det_facts_sp$C_P
  C_D <- det_facts_sp$C_D
  
  # Number of samples
  n_samples <- 500
  
  # Get growth rate
  m_est_r <- coef(m)["year"]; names(m_est_r) <- "r"
  m_est_r_se <- summary(m, re.test = F)$se[2]; names(m_est_r_se) <- "se(r)"
  m_est_r_samples_link <- rnorm(n = n_samples, mean = m_est_r, sd = m_est_r_se)
  m_est_r_samples <- exp(m_est_r_samples_link)-1
  
  # Get count intercept and adjust for detection
  m_est_y0 <- coef(m)["(Intercept)"]; names(m_est_y0) <- "y0"
  m_est_y0_se <- summary(m, re.test = F)$se[1]; names(m_est_y0_se) <- "se(y0)"
  N <- exp(rnorm(n = n_samples, mean = m_est_y0, sd = m_est_y0_se))
  m_est_N0_samples <- log(N/(50 * pi * 400^2) * ((400/C_D)^2) * C_P * (1.5e7/1))
  
  result <- data.frame(
    r = m_est_r_samples,
    N0 = m_est_N0_samples,
    species = spp[i]) %>%
    mutate(sample = row_number())
  
  # Save samples
  species_samples[[i]] <- result
  
  # Save point estimates
  species_estimates[[i]] <- result %>%
    group_by(species) %>%
    summarise(r_mean = mean(r),
              r_lwr = quantile(r, 0.025),
              r_upr = quantile(r, 0.975),
              N0_mean = mean(N0),
              N0_lwr = quantile(N0, 0.025),
              N0_upr = quantile(N0, 0.975)) %>%
    ungroup()
  
  
  # Variances
  variance_initial_abundance <- var(m_est_N0_samples)
  variance_growth_rate <- var(m_est_r_samples)
  
  species_parameter_variances[[i]] <- data.frame(
    var_initial_abundance = variance_initial_abundance,
    var_growth_rate = variance_growth_rate,
    species = spp[i])
  
  
  # ---- Prepare to predict fitted dynamics ----
  
  count_to_density <- function(N) N/(50 * pi * 400^2) * ((400/C_D)^2) * C_P * (1.5e7/1)
  
  biomass_i <- species_with_masses %>%
    filter(species == spp[i]) %>%
    pull(Mass)
  
  pred_df <- df %>%
    select(year, year_factor, 
           route, observer, 
           first_year, duration, yday) %>%
    mutate(duration = 0, yday = 0, 
           route = unique(route)[1],
           observer = unique(observer)[1],
           first_year = 0) %>%
    distinct() %>%
    arrange(year)
  
  
  # ---- Predict annual dynamics - SIMS ----
  
  species_preds_sims[[i]] <- pred_df %>%
    fitted_values(m, ., exclude = c("s(observer)", "s(route)", "s(duration)", "s(yday)", "first_year"),
                  scale = "link") %>%
    mutate(year = as.numeric(as.character(year_factor))) %>%
    mutate(species = spp[i]) %>%
    select(species, year, .fitted, .lower_ci, .upper_ci, .se) %>%
    rowwise() %>%
    mutate(sim_logC = list(rnorm(n = 1000, mean = .fitted, sd = .se))) %>%
    unnest(sim_logC) %>%
    mutate(sim_C = exp(sim_logC),
           sim_D = count_to_density(sim_C),
           sim_B = biomass_i * sim_D) %>%
    select(species, year, sim_C, sim_D, sim_B) %>%
    group_by(species, year) %>%
    mutate(sim = row_number()) %>%
    ungroup()
  
  
  # ---- Predict annual dynamics - FITS ----
  
  species_preds[[i]] <- pred_df %>%
    fitted_values(m, ., exclude = c("s(observer)", "s(route)", "s(duration)", "s(yday)", "first_year"), 
                  scale = "response") %>%
    mutate(year = as.numeric(as.character(year_factor))) %>%
    mutate(species = spp[i]) %>%
    
    select(species, year, 
           count_fit = .fitted, 
           count_lwr = .lower_ci, 
           count_upr = .upper_ci) %>%
    
    mutate(density_fit = count_to_density(count_fit),
           density_lwr = count_to_density(count_lwr), 
           density_upr = count_to_density(count_upr)) %>%
    
    mutate(biomassDensity_fit = biomass_i * density_fit,
           biomassDensity_lwr = biomass_i * density_lwr,
           biomassDensity_upr = biomass_i * density_upr)
  
  # ggplot(preds, aes(x = year, y = count_fit)) +
  #   geom_pointrange(aes(ymin = count_lwr, ymax = count_upr)) +
  #   coord_cartesian(ylim = c(0,NA))
  # 
  # ggplot(preds, aes(x = year, y = density_fit)) +
  #   geom_pointrange(aes(ymin = density_lwr, ymax = density_upr)) +
  #   coord_cartesian(ylim = c(0,NA))
  # 
  # ggplot(preds, aes(x = year, y = biomassDensity_fit)) +
  #   geom_pointrange(aes(ymin = biomassDensity_lwr, ymax = biomassDensity_upr)) +
  #   coord_cartesian(ylim = c(0,NA))
  
  
  # ---- Calculate initial and final density ----
  
  species_FirstLastDensities[[i]] <- calculate_first_last_density_perArea(
    m = m, species = spp[i],
    C_D = C_D, 
    C_P = C_P,
    area_original = area_original, 
    area_wanted = 1, 
    min_year = min_year)
  
  # Capture and store warnings generated in this iteration
  iter_warnings <- warnings()
  if (length(iter_warnings) > 0) {
    # Store warnings in the list with the iteration number as the identifier
    warnings_by_iteration[[as.character(i)]] <- iter_warnings
  }
  
  # Reset warning state for the next iteration
  options(warn = 1)
  
  pb$tick()
  
}

# Save results
# sp_to_rm_from_lm2 <- do.call(rbind, species_estimates) %>% filter(N0_mean < log(0.01)) %>% pull(species)
sp_to_rm_from_lm2 <- NA
species_samples_df <- do.call(rbind, species_samples) %>% filter(!(species %in% sp_to_rm_from_lm2))
species_estimates_df <- do.call(rbind, species_estimates) %>% filter(!(species %in% sp_to_rm_from_lm2))
species_variances_df <- do.call(rbind, species_parameter_variances) %>% filter(!(species %in% sp_to_rm_from_lm2))
species_preds_df <- do.call(rbind, species_preds) %>% mutate(BCR = bcr_code) %>% filter(!(species %in% sp_to_rm_from_lm2))
species_preds_sims_df <- do.call(rbind, species_preds_sims) %>% mutate(BCR = bcr_code) %>% filter(!(species %in% sp_to_rm_from_lm2))
species_FirstLast_df <- do.call(rbind, species_FirstLastDensities) %>% filter(!(species %in% sp_to_rm_from_lm2))


# ---- Get resampled type 2 regression ----

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent elapsed: :elapsed eta: :eta",
  total = 500, clear = FALSE)

# Get samples from each regression
regr_samples <- list()
for(i in 1:500){
  
  # Select data
  df_i <- species_samples_df %>%
    group_by(species) %>%
    slice(i) %>%
    ungroup()
  
  # Fit model
  mod <- lmodel2(
    formula = r ~ N0,
    data = df_i,
    range.x = "interval", 
    range.y = "interval",
    nperm = 99)
  
  # Get estimated means
  means <- mod$regression.results %>%
    filter(Method == "RMA")
  mean_intercept <- means$Intercept
  mean_slope <- means$Slope
  
  # Get estimated intervals
  intervals <- mod$confidence.intervals %>%
    filter(Method == "RMA")
  upper_intercept <- intervals["97.5%-Intercept"]
  upper_slope <- intervals["97.5%-Slope"]
  
  # Get standard errors
  se_intercept <- as.numeric((upper_intercept - mean_intercept) / 1.96)
  se_slope <- as.numeric((upper_slope - mean_slope) / 1.96)
  
  # Type2Regr Intercept samples
  regr_intercept_samples <- rnorm(
    n = 100, 
    mean = mean_intercept, 
    sd = se_intercept)
  
  # Type2Regr Slope samples
  regr_slope_samples <- rnorm(
    n = 100, 
    mean = mean_slope, 
    sd = se_slope)
  
  # Type2Regr output
  results <- data.frame(
    intercept = regr_intercept_samples,
    slope = regr_slope_samples) %>%
    mutate(sample = i)
  
  # Save results
  regr_samples[[i]] <- results
  
  pb$tick()
  
}

regr_samples_df <- do.call(rbind, regr_samples)

hist(regr_samples_df$intercept, breaks = 50)
hist(regr_samples_df$slope, breaks = 50)


# ---- Calculating regressions ----

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent elapsed: :elapsed eta: :eta",
  total = nrow(regr_samples_df), clear = FALSE)

range_lwr <- min(species_estimates_df$N0_mean)
range_upr <- max(species_estimates_df$N0_mean)
x <- seq(range_lwr, range_upr, length.out = 100)

y_results <- list()
for(i in 1:nrow(regr_samples_df)){
  
  beta0 <- regr_samples_df$intercept[i]
  beta1 <- regr_samples_df$slope[i]
  y <- beta0 + beta1 * x
  
  y_df <- data.frame(y, x) %>%
    mutate(sample = i)
  
  y_results[[i]] <- y_df
  
  pb$tick()
}

y_results_df <- do.call(rbind, y_results)

final_relationship <- y_results_df %>%
  group_by(x) %>%
  summarise(mean = mean(y),
            lwr = quantile(y, 0.025),
            upr = quantile(y, 0.975)) %>%
  ungroup()

ggplot(final_relationship, aes(x = x, y = mean)) +
  
  geom_hline(yintercept = 0) +
  
  geom_errorbar(data = species_estimates_df,
                aes(x = N0_mean, y = r_mean, 
                    ymin = r_lwr, ymax = r_upr), 
                color = "gray70") +
  
  geom_errorbarh(data = species_estimates_df, 
                 aes(x = N0_mean, y = r_mean,  
                     xmin = N0_lwr, xmax = N0_upr), 
                 color = "gray70") +
  
  geom_point(data = species_estimates_df,
             aes(x = N0_mean, y = r_mean), 
             color = "gray70", size = 1) +
  
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "purple", alpha = 0.5) +
  geom_line(col = "purple", linewidth = 1.5) +
  
  scale_x_continuous(breaks = c(log(0.01),log(0.1),log(1),log(10),log(100)),
                     labels = c(0.01, 0.1, 1, 10, 100)) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  
  labs(x = "Former abundance per 1500 hectares ", 
       y = "Annual growth rate",
       subtitle = bird_conservation_region) +
  
  theme_light(16) +
  theme(aspect.ratio = 1,
        text = element_text(family = "Lato"),
        panel.grid.minor = element_blank())


# ---- Results ----

# Combine your data objects into a single list
combined_list <- list(
  FinalRelationship = final_relationship,
  SpeciesSamples = species_samples_df, 
  SpeciesEstimates = species_estimates_df,
  SpeciesVariances = species_variances_df,
  SpeciesTimeSeries = species_preds_df,
  SpeciesTimeSeriesSims = species_preds_sims_df,
  FirstLastDensities = species_FirstLast_df
)

# Save the combined list as an RDS file in your directory
BCR_output_file <- paste("data/BCR-results/results_bcr_", BCR_file_suffix, ".rds", sep = "")
saveRDS(combined_list, file = BCR_output_file)

