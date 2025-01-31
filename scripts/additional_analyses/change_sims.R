library(tidyverse)
library(mgcv)
# set.seed(1)


# ---- Survey setup ----

n_sites <- 10
n_years <- 50
n_species <- 80
C <- array(NA, dim = c(n_years, n_sites, n_species))


# ---- Simulate counts ----

beta0 <- c()
beta1 <- c()
for(spp in 1:n_species){
  
  beta0[spp] <- rnorm(n = 1, mean = -1, sd = 2.5)
  beta1[spp] <- rnorm(n = 1, mean = 0, sd = 0.02)
  pln_noise <- rexp(n = 1, rate = 10)
  
  for(s in 1:n_sites){
    
    eps_site <- rnorm(n = 1, mean = 0, sd = 1)
    
    for(t in 1:n_years){
      
      log_mu <- (beta0[spp] + eps_site) + beta1[spp] * (t-1)
      # C[t,s,spp] <- rpois(n = 1, lambda = exp(log_mu))
      log_lambda <- rnorm(n = 1, mean = log_mu, sd = pln_noise)
      lambda <- exp(log_lambda)
      C[t,s,spp] <- rpois(n = 1, lambda = lambda)
    }
  }
}


# ---- Convert C matrix to data frame list for modeling ----

df_list <- list()
for(spp in 1:n_species){
  df_list[[spp]] <- C[,,spp] %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(year = row_number()-1) %>%
    pivot_longer(cols = -year, names_to = "site", values_to = "count") %>%
    mutate(species = spp,
           site = as.factor(site))
}


# ---- Recover beta0 and beta1 estimates ----

hat_beta0 <- c()
hat_beta1 <- c()
for(spp in 1:n_species){

  m <- bam(
    formula = count ~ year + s(site, bs = "re"),
    # formula = count ~ year, 
    family = "poisson", 
    data = df_list[[spp]], 
    discrete = T,
    gamma = 1.4)
  
  hat_beta0[spp] <- coef(m)[1]
  hat_beta1[spp] <- coef(m)[2]
  
}


# ---- Plot results ----

data.frame(
  N = hat_beta0,
  delta = exp(hat_beta1)-1) %>%
  filter(delta >= -0.05, delta <= 0.05,
         N >= log(0.0001), N <= log(1000)) %>%
  ggplot(aes(x = N, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  geom_point() +
  scale_x_continuous(
    breaks = c(log(0.0001),log(0.001), log(0.01),log(0.1),log(1),log(10), log(100), log(1000)),
    labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000")) +
  scale_y_continuous(limits = c(-0.05, 0.05), labels = scales::percent_format()) +
  theme(aspect.ratio = 1)


# ---- Fit a community-level model ----

fit_corr <- glmmTMB(
  count ~ year + 
    (1 + year | species) + 
    (1 | species:site),
  family = nbinom2,
  data = df
)

fit_nocorr <- glmmTMB(
  count ~ year + 
    (1 | species) + 
    (0 + year | species) + 
    (1 | species:site),
  family = nbinom2,
  data = df
)

AIC(fit_corr, fit_nocorr)
anova(fit_nocorr, fit_corr)








