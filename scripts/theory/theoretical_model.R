library(tidyverse)
library(lmodel2)
library(deSolve)
library(ggthemes)
library(ggpmisc)
library(readxl)
library(sn)
source("scripts/00_Functions.R")
# sci_names_df <- get_sci_names()
# bbsdat <- bbsBayes2::load_bbs_data(release = "2023")
sp_hash <- bbsdat$species %>%
  select(aou, species = english) %>%
  distinct()

set.seed(0)


# ---- MTE functions ----

r <- function(r_0, M){
  
  r_0 * (M^(-0.25))
  
}

K <- function(K_0, M){
  
  K_0 * (M^(-0.82))
  
} 


# ---- Mass data ----

# Get the species
the_species <- read.csv("data/allspecies.csv") %>% 
  select(-1) %>%
  as_tibble() %>%
  left_join(sci_names_df, by = join_by("species_clean" == "Common Name")) %>%
  distinct()

# Get mass distribution
spp_masses <- read_excel("data/AVONET2_eBird.xlsx", sheet = 2) %>%
  select(sci_name = Species2, mass = Mass) %>%
  semi_join(the_species, by = join_by(sci_name)) %>%
  left_join(the_species, by = join_by(sci_name)) %>%
  select(species = species_bcrForm, mass) %>%
  distinct() %>%
  pull(mass)

# Log masses
log_spp_masses <- log(spp_masses)

# I use this a bit below to simulate the body mass distribution
fit <- selm(log_spp_masses ~ 1, family = "SN")


# ---- Model -----

# Define the model for logistic growth of each species
spp_mod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dN <- numeric(n_species)
    
    for (i in 1:n_species) {
      
      N <- state[i]
      r <- parameters[["r"]][i]
      
      K_initial <- parameters[["K"]][i]
      K_decay_rate <- parameters[["K_decay_rate"]][i]
      # K <- K_initial * exp(K_decay_rate * t)
      K <- K_initial * ((1 + K_decay_rate) ^ t)
      
      dN[i] <- r * N * ((K - N) / K)
    }
    
    return(list(dN))
  })
}


# ---- Population parameters ----

# Bird species
n_species <- 10000

# Mass in kilograms
Mass <- rsn(
  n = n_species, 
  xi = coef(fit, "DP")[1], 
  omega = coef(fit, "DP")[2], 
  alpha = coef(fit, "DP")[3]) %>%
  exp()

# Intrinsic growth rate
r_values <- r(r_0 = 0.1, M = Mass)

# Carrying capacity
K_values <- K(K_0 = 1920.1, M = Mass)

# Initial abundance
N0_values <- runif(n_species, min = 0.0001, max = 1.1) * K_values

# Decay rate for carrying capacity (4% per year)
K_decay_rate <- rnorm(n_species, mean = -0.014, sd = 0.015)

# Hash table
spp_params_hash <- data.frame(
  Species_numeric = 1:n_species,
  initNfromK = N0_values/K_values
)


# ----- Organize parts of running model ----

# Create named vectors for initial state and parameters
initial_state <- c(N0_values)
parameters <- list(
  r = r_values, 
  K = K_values,
  K_decay_rate = K_decay_rate)

# Time points to solve the ODE
observed_times <- seq(0, 56, by = 1)


# ---- Run model ----

# Run the model and organize results
model_output <- ode(y = initial_state, times = observed_times, func = spp_mod, parms = parameters) %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "Species", values_to = "Population") %>%
  mutate(Species = as.factor(Species))

if(F) {
  
  ggplot(model_output, aes(x = time, y = Population, color = Species)) +
    # facet_wrap(~Species, scales = "free") +
    # geom_line() +
    scale_y_log10() +
    theme(legend.position = "none")
  
}


# ---- Analyze results for common vs rare ----

# Function to get initial abundance and growth rate for each species
get_ests <- function(y,x){
  
  m <- glm(y ~ x, family = gaussian(link = "log"))
  
  N0 <- exp(coef(m)[1])
  r <- exp(coef(m)[2])-1
  
  result <- c(N0, r)
  
  return(result)
  
}

# Calculate trends to get initial abundance and growth rate for each species
recovered_trends <- model_output %>%
  group_by(Species) %>%
  arrange(time) %>%
  summarise(N0_true = first(Population),
            N0_hat = get_ests(Population, time)[1],
            r_hat = get_ests(Population, time)[2]) %>%
  ungroup() %>%
  mutate(log_N0_hat = log(N0_hat)) %>%
  mutate(Species_numeric = as.numeric(as.character(Species))) %>%
  arrange(Species_numeric) %>%
  left_join(spp_params_hash, by = "Species_numeric")


# ---- Plot the results ----

# Plot species N0_hat and r_hat
ggplot(recovered_trends, aes(x = N0_hat, y = r_hat)) +
  geom_hline(yintercept = 0, col = 1) +
  geom_point(aes(color = initNfromK), size = 3, alpha = 1) +
  stat_ma_line(method = "RMA", range.x = "interval", range.y = "interval",
               fill = 1, color = 1, se = F, linewidth = 2) +
  labs(x = "Former abundance", 
       y = "Population trend") +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(NA, 0.07),
                     breaks = c(0.07, 0.05,0.03,0.01,-0.01,-0.03,-0.05,-0.07)) +
  scale_color_viridis_c("Initial distance from carrying capacity", option = "H") +
  theme_light(15) +
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.5, 
                                ticks.colour = NA)) +
  theme(aspect.ratio = 1,
        text = element_text(family = "Roboto"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(angle = 0, vjust = 1),
        legend.key.height = unit(0.3, "cm"),  # Thinner bar
        legend.key.width = unit(2, "cm"))

