library(tidyverse)
library(lmodel2)
library(deSolve)
library(ggthemes)
library(ggpmisc)

set.seed(0)



# ---- Model -----

# Define the logistic growth model for a species
spp_mod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Logistic growth equation
    dN <- r * N * (1 - (N / K))
    return(list(dN))
  })
}


# ---- Parameters ----

# Parameters and initial state

state_rare <- c(N = 0.05)
state_common <- c(N = 0.98)
params_Kstatic <- list(r = 0.065, K = 0.98)
params_Kdrop <- list(r = 0.065, K = 0.7)

times <- seq(0, 50, by = 1)

# Solve the ODE
output_rare_Kstatic   <- ode(y = state_rare, times = times, func = spp_mod, parms = params_Kstatic)
output_rare_Kdrop     <- ode(y = state_rare, times = times, func = spp_mod, parms = params_Kdrop)
output_common_Kstatic <- ode(y = state_common, times = times, func = spp_mod, parms = params_Kstatic)
output_common_Kdrop   <- ode(y = state_common, times = times, func = spp_mod, parms = params_Kdrop)


# ---- Compile results ----

results <- rbind(
  output_rare_Kstatic %>% as_tibble() %>% mutate(species = "Rare", K = "static"),
  output_rare_Kdrop %>% as_tibble() %>% mutate(species = "Rare", K = "drop"),
  output_common_Kstatic %>% as_tibble() %>% mutate(species = "Common", K = "static"),
  output_common_Kdrop %>% as_tibble() %>% mutate(species = "Common", K = "drop")) %>%
  mutate(time = as.numeric(time), N = as.numeric(N)) %>%
  mutate(species = factor(species, levels = c("Rare", "Common")),
         K = factor(K, levels = c("static", "drop")))

ggplot(results, aes(x = time, y = N, group = K, color = interaction(species, K))) +
  facet_wrap(~species) +
  geom_hline(yintercept = 0, color = NA) +
  geom_line(linewidth = 2, lineend = "round") +
  scale_color_manual(values = c(
    "Rare.static" = "skyblue",
    "Rare.drop" = "blue3",
    "Common.static" = "pink",
    "Common.drop" = "red3")) +
  geom_hline(yintercept = 1.0, linewidth = 1, color = "gray80", linetype = "longdash") +
  geom_hline(yintercept = 0.7, linewidth = 1, color = "black", linetype = "longdash") +
  geom_segment(
    aes(x = 45, y = 0.96, xend = 45, yend = 0.74),
    lineend = "round", linejoin = "mitre",
    arrow = arrow(length = unit(0.1, "inches")), 
    color = "black", 
    linewidth = 0.6
  ) +
  annotate(
    geom = "text",
    x = 43, y = 0.865, hjust = 1,
    label = "K", parse = T,
    size = 6, color = "black") +
  labs(x = "Year", y = "Abundance") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_light(18) +
  theme(
    axis.text = element_text(color = 1),
    axis.title = element_text(color = 1),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, color = 1, linewidth = 0.5),
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = 1),
    panel.grid = element_blank())



