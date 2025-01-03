# BBS data
library(tidyverse)
library(bbsBayes2)
library(readxl) 
library(mgcv)
library(lubridate)
library(progress)
library(vroom)
library(ggpmisc)
# Define the function for density calculation
count_to_density <- function(N, C_D, C_P) {
  N / (50 * pi * 400^2) * ((400 / C_D)^2) * C_P * (1.5e7 / 1)
}


# ---- SELECT BCR ----

bird_conservation_region <- "Shortgrass Prairie"

# Hash table for BCR names and numbers
bcrs <- data.frame(
  region = c("Atlantic Northern Forest","Northern Pacific Rainforest",
             "Appalachian Mountains", "Eastern Tallgrass Prairie",
             "Northern Rockies", "Shortgrass Prairie"),
  number = c(14, 5, 28, 22, 10, 18)) %>%
  mutate(bcr_code = paste0("BCR", number))

# Now retrieving number
bcr_number <- bcrs %>% 
  filter(region == bird_conservation_region) %>%
  pull(number)


# ---- Model stuff ----

# Function to attempt fitting the model with negative binomial, fallback to Poisson on warning
fit_model_with_fallback <- function(data) {
  tryCatch({
    # Try fitting with negative binomial
    model <- bam(species_total ~ year*latitude + 
                   s(yday, k = 4) +
                   s(duration, k = 4) +
                   s(year_factor, bs = "re") +
                   s(route, bs = "re") + 
                   s(observer, bs = "re"),
                 family = nb(),
                 data = df_m,
                 gamma = 1.4,
                 discrete = T)
    return(model)
  }, warning = function(w) {
    # If a warning is caught, fit with Poisson
    model <- bam(species_total ~ year*latitude + 
                   s(yday, k = 4) +
                   s(duration, k = 4) +
                   s(year_factor, bs = "re") +
                   s(route, bs = "re") + 
                   s(observer, bs = "re"),
                 family = poisson(),
                 data = df_m,
                 gamma = 1.4,
                 discrete = T)
    return(model)
  })
}


# ---- BBS data ----

# Load all BBS data
dat <- load_bbs_data(level = "state", release = 2023, sample = F, quiet = F)

# Detectability factors
detectability_factors <- read_xlsx("data/detectability-factors/detectability_factors.xlsx")

# Select species - species in study
species <- toupper(bird_conservation_region) %>%
  str_replace_all(" ", "_") %>%
  paste0("data/BCR-species/species_bcr_", ., ".csv") %>%
  vroom() %>%
  select(species) %>%
  left_join(dat$species, by = c("species" = "english")) %>%
  select(species, aou) %>%
  distinct()

# AOU codes
spp_aou <- species$aou

# Number of species
nspp <- length(spp_aou)

# Routes data
routes_dat <- dat$routes %>%
  filter(bcr == bcr_number) %>%
  mutate(route = str_c(country_num, state_num, route, sep = "-")) %>%
  select(route, year, observer = obs_n, month, day, latitude, start_time, end_time) %>%
  mutate(start_time = strptime(sprintf("%04d", start_time), format = "%H%M"),
         end_time = strptime(sprintf("%04d", end_time), format = "%H%M")) %>%
  mutate(duration = as.numeric(end_time - start_time)) %>%
  mutate(date = make_date(year, month, day),
         yday = yday(date)) %>%
  select(year, route, observer, yday, duration, latitude) %>%
  group_by(observer) %>%
  mutate(first_year = if_else(year == min(year), 1, 0)) %>%
  ungroup()

# Bird data
bird_dat <- dat$birds %>%
  filter(bcr == bcr_number) %>%
  filter(aou %in% spp_aou) %>%
  select(country_num, state_num, route, year, aou, species_total) %>%
  mutate(route = str_c(country_num, state_num, route, sep = "-")) %>%
  select(route, year, aou, species_total) %>%
  complete(nesting(route, year), aou, fill = list(species_total = 0)) %>%
  select(species_total, aou, year, route) %>%
  left_join(routes_dat, by = join_by(route, year)) %>%
  mutate(yday = scale(yday)[,1],
         latitude = scale(latitude)[,1],
         duration = scale(duration)[,1],
         year_factor = as.factor(year),
         observer = as.factor(as.character(observer)),
         year = year - min(year),
         route = as.factor(route)) %>%
  na.omit()


# ---- Loop through species ----

pb <- progress_bar$new(
  format = "[:bar] :percent elapsed: :elapsed eta: :eta",
  total = nspp, clear = FALSE)

ests <- list()
for(i in 1:nspp){
  
  df_m <- bird_dat %>%
    filter(aou == spp_aou[i])
  
  m <- fit_model_with_fallback(df_m)
  
  summary_m <- summary(m, re.test = F)
  
  ests[[i]] <- data.frame(
    aou = spp_aou[i],
    term = c("latitude", "year:latitude","AARPC","InitN"),
    estimate = c(
      coef(m)["latitude"],
      coef(m)["year:latitude"],
      exp(coef(m)["year"])-1,
      exp(coef(m)["(Intercept)"])),
    pv = c(
      summary_m$p.pv["latitude"],
      summary_m$p.pv["year:latitude"],
      NA,
      NA)) %>%
    `rownames<-`(NULL) %>%
    mutate(pv = if_else(pv < 0.05, "sig.", "insig.")) %>%
    as_tibble() %>%
    left_join(detectability_factors %>% select(2:4),
              by = join_by(aou)) %>%
    mutate(estimate = if_else(
      term == "InitN", 
      count_to_density(estimate, C_D, C_P), 
      estimate)) %>%
    select(-c(C_D,C_P))
  
  pb$tick()
  
}

ests_df <- do.call(rbind, ests)

ests_df_clean <- ests_df %>%
  select(aou, term, estimate, pv) %>%
  distinct() %>%
  filter(!(term %in% c("AARPC", "InitN"))) %>% 
  mutate(estimate = if_else(pv == "insig.", 0, estimate)) %>%
  pivot_wider(id_cols = aou, names_from = term, values_from = estimate) %>%
  select(aou, latitude, yearXlatitude = 3) %>%
  mutate(bcr = bird_conservation_region)
  
path <- paste0("data/RangeShift-results/RangeShift_",str_replace_all(bird_conservation_region, " ", "_"),".csv")
write.csv(ests_df_clean, path)



