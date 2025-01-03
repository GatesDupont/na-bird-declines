# Load necessary libraries
library(ape)
library(auk)
library(parallel)
library(MuMIn)
library(tidyverse)
library(MRFtools)
library(readxl)
library(mgcv)
source("scripts/00_Functions.R")
sci_names_df <- get_sci_names()
bbsdat <- bbsBayes2::load_bbs_data(release = "2023")
sp_hash <- bbsdat$species %>%
  select(aou, species = english) %>%
  distinct()


# ---- Mass data ----

# Get the species
the_species <- read.csv("data/allspecies.csv") %>% 
  select(-1) %>%
  as_tibble() %>%
  left_join(sci_names_df, by = join_by("species_clean" == "Common Name")) %>%
  distinct()

# Get mass distribution
Mass <- read_excel("data/AVONET2_eBird.xlsx", sheet = 2) %>%
  select(sci_name = Species2, Mass, 
         trophic_level = Trophic.Level,
         trophic_niche = Trophic.Niche) %>%
  semi_join(the_species, by = join_by(sci_name)) %>%
  left_join(the_species, by = join_by(sci_name)) %>%
  mutate(log_mass = log(Mass)) %>%
  select(species = species_bcrForm, log_mass) %>%
  mutate(species = if_else(species == "Northern House Wren", "House Wren", species))


# ---- Loop -----

df_bcr <- list()
for(i in 1:6){
  
  # BCR names
  bird_conservation_region <- c("Atlantic Northern Forest","Northern Pacific Rainforest",
                                "Appalachian Mountains", "Eastern Tallgrass Prairie",
                                "Northern Rockies", "Shortgrass Prairie")[i]
  
  # ---- Load in model results ----
  
  # Read the combined list as an RDS file in your directory
  BCR_file_suffix <- gsub(" ", "_", toupper(bird_conservation_region))
  BCR_output_file <- paste("data/BCR-results/results_bcr_", BCR_file_suffix, ".rds", sep = "")
  combined_results <- readRDS(file = BCR_output_file)
  
  # ---- Broad Groups categories -----
  file_name <- "data/BCR-species/complete/species_bcr_APPALACHIAN_MOUNTAINS_categories.xlsx"
  categories_dataset <- paste("data/BCR-species/complete/species_bcr_", BCR_file_suffix, "_categories.xlsx", sep = "")
  broad_groups <- read_excel(categories_dataset)
  
  if("Sagebrush_obligates" %in% colnames(broad_groups)){
    broad_groups <- broad_groups %>% select(-Sagebrush_obligates)
  }
  
  # ---- Load range shift stuff ----
  
  path <- paste0("data/RangeShift-results/RangeShift_",str_replace_all(bird_conservation_region, " ", "_"),".csv")
  df_RangeShift <- read.csv(path) %>%
    select(-X)
  
  df_RangeShift_wSppp <- combined_results$SpeciesEstimates %>%
    left_join(sp_hash, by = join_by(species)) %>%
    select(species, aou) %>%
    left_join(df_RangeShift, by = join_by(aou)) %>%
    select(-c(aou, bcr))
  
  
  # ---- Assign habitats to species ----
  
  # Pull out important data
  df_bcr[[i]] <- combined_results$SpeciesEstimates %>%
    # Append variance of estimates and make weights
    left_join(combined_results$SpeciesVariances, by = "species") %>%
    rename(Variance = var_growth_rate) %>%
    mutate(weight = 1 / Variance) %>% # Calculate weights
    mutate(weight = round(weight/min(weight), 0)) %>%
    # Keep only required data
    select(species, growth_rate = r_mean, Initial_abundance = N0_mean, weight) %>%
    # Append my species groups
    mutate(species_clean = sub(" \\(all forms\\)", "", species)) %>%
    left_join(broad_groups, by = "species") %>%
    mutate(Migratory_species = if_else(Residents == 1, 0, 1)) %>%
    left_join(df_RangeShift_wSppp, by = join_by(species)) %>%
    na.omit()
  
}

# Load PIF data for climate and ag data
PIF_dat <- read_xlsx("/Users/gatesdupont/Desktop/databases/ACAD Global 2024.05.23.xlsx") %>%
  select(species_clean = `Common Name`, 
         agriculture = Agriculture, 
         urban = Urban,
         clim_vuln = `CV-b`)
  # mutate(agriculture = if_else(str_detect(agriculture, "b"), 1, NA))

# Combine the datasets                                                                                                                 
df <- do.call(rbind, df_bcr) %>%
  mutate(bcr = as.factor(bcr),
         species = as.factor(species)) %>%
  left_join(Mass, by = join_by(species)) %>%
  left_join(PIF_dat, by = c("species_clean")) %>%
  mutate(agriculture = if_else(is.na(agriculture), 0, 1)) %>%
  mutate(urban = if_else(is.na(urban), 0, 1))

# This is from the 03_drivers_combinedBCRs script
my_species <- df %>%
  select(species) %>%
  distinct() %>%
  mutate(species_clean = gsub(" \\(all forms\\)", "", species)) %>%
  mutate(sci_name = ebird_species(.$species_clean)) %>%
  mutate(sci_name_dash = gsub(" ", "_", sci_name))

# Read the phylogenetic tree - I think this is Clements 4 taxonomy
phylo_tree <- read.tree("/Users/gatesdupont/Desktop/databases/Stage2_Hackett_MCC_no_neg.tre")

# Data frame of tree tip names (seems to be from Clements 4)
phylo_tree_spp <- data.frame(
  sci_name_dash = phylo_tree$tip.label,
  exists = 1)

# Converting 2023 eBird taxonomy to Clements 4 taxonomy, manually
species_names_hash <- my_species %>%
  mutate(sci_name_dash = gsub("Setophaga_americana", "Parula_americana", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Setophaga_citrina", "Wilsonia_citrina", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Setophaga", "Dendroica", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Dendroica_ruticilla", "Setophaga_ruticilla", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Poecile", "Parus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Haemorhous", "Carpodacus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Antrostomus", "Caprimulgus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Cardellina", "Wilsonia", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Spinus", "Carduelis", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Dryobates", "Picoides", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Leiothlypis", "Vermivora", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Parkesia", "Seiurus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Melozone_fusca", "Pipilo_fuscus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Melozone", "Pipilo", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Rhynchophanes", "Calcarius", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Corthylio", "Regulus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Centronyx", "Ammodramus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Geothlypis_philadelphia", "Oporornis_philadelphia", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Vermivora_cyanoptera", "Vermivora_pinus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Troglodytes_pacificus", "Troglodytes_troglodytes", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Troglodytes_hiemalis", "Troglodytes_troglodytes", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Spatula_cyanoptera", "Anas_cyanoptera", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Circus_hudsonius", "Circus_cyaneus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Artemisiospiza_nevadensis", "Amphispiza_belli", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Peucaea_cassinii", "Aimophila_cassinii", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Ardea_alba", "Casmerodius_albus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Cistothorus_stellaris", "Cistothorus_platensis", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Geothlypis_formosa", "Oporornis_formosus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Geothlypis_tolmiei", "Oporornis_tolmiei", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Ixoreus_naevius", "Zoothera_naevia", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Selasphorus_calliope", "Stellula_calliope", sci_name_dash)) %>%
  left_join(phylo_tree_spp) %>%
  select(species, tree_sci_name = sci_name_dash)

# Final data
df2 <- df %>%
  left_join(species_names_hash) %>%
  mutate(tree_sci_name = as.factor(tree_sci_name)) %>%
  mutate(Initial_abundance = scale(Initial_abundance)[,1])

# Pruned tree
pruned <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, df2$tree_sci_name))

# Penalty matrix from phylogeny
pmat <- mrf_penalty(pruned, type = "individual")

# Make sure the scientific species names are in the same order as the tip labels
spp_names_ordered_ref <- df2 %>%
  mutate(tree_sci_name = factor(tree_sci_name, levels = pruned$tip.label)) %>%
  select(tree_sci_name, species) %>%
  distinct() %>%
  arrange(tree_sci_name) %>%
  pull(species)

# Make sure the scientific species names are in the same order as the tip labels
df2 <- df2 %>%
  mutate(tree_sci_name = factor(tree_sci_name, levels = pruned$tip.label)) %>%
  arrange(tree_sci_name) %>%
  mutate(species = factor(species, levels = spp_names_ordered_ref)) %>%
  mutate(Initial_abundance = scale(Initial_abundance)[,1],
         log_mass = scale(log_mass)[,1],
         clim_vuln = scale(clim_vuln)[,1],
         latitude = scale(latitude)[,1],
         yearXlatitude = scale(yearXlatitude)[,1])


# ---- Global model ----

# # Global model
# m.global <- bam(growth_rate ~ Introduced_species + Aerial_insectivores +
#                   Secondary_growth_specialists + Urban_adapted_species +
#                   Forest_obligates + Grassland_specialists + Wetland_obligates + 
#                   Obligate_migrants + Facultative_migrants +
#                   Tree_cavity_nesters + log_mass +
#                   latitude * yearXlatitude +
#                   s(bcr, bs = "re") +
#                   s(tree_sci_name, bs = "mrf", xt = list(penalty = pmat)),
#                 discrete = T, select = F, gamma = 1.4,
#                 weights = df2$weight, data = df2, na.action = "na.fail")

# Global model
m.global <- bam(growth_rate ~ log_mass +
                  Aerial_insectivores + Obligate_migrants +
                  Urban_adapted_species + Forest_obligates +
                  Grassland_specialists + agriculture + 
                  latitude * yearXlatitude +
                  s(bcr, bs = "re") +
                  s(tree_sci_name, bs = "mrf", xt = list(penalty = pmat)),
                discrete = T, select = F, gamma = 1.4,
                weights = df2$weight, data = df2, na.action = "na.fail")

summary(m.global)


# ----- Use global model for model selection ----

# Set up the cluster and fit models
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))
clusterEvalQ(clust, library(mgcv))
options(na.action = "na.fail")
clusterExport(clust, c("pmat", "df2"))
model_set <- dredge(m.global, cluster = clust, trace = 2, 
                    fixed = ~ s(bcr, bs = "re") + s(tree_sci_name, bs = "mrf", xt = list(penalty = pmat)))
stopCluster(clust)


# ---- Mod avg ----

mavg <- model.avg(model_set, subset = delta <= 2)
# mavg <- model.avg(model_set, subset = cumsum(weight) <= 0.95)

model_set %>%
  as_tibble() %>%
  mutate(rank = row_number()) %>%
  select(delta, weight) %>%
  mutate(cs = cumsum(weight)) %>%
  mutate_all(.funs = function(x) round(x, 3))
  
cov_npp <- df %>%
  select(Introduced_species, Aerial_insectivores, urban,
         Secondary_growth_specialists, Urban_adapted_species, 
         Forest_obligates, Grassland_specialists, Wetland_obligates, 
         Obligate_migrants, Facultative_migrants, Residents,
         Tree_cavity_nesters, log_mass, clim_vuln, agriculture, 
         latitude, yearXlatitude) %>%
  mutate(`latitude:yearXlatitude` = latitude * yearXlatitude) %>%
  pivot_longer(cols = everything()) %>%
  mutate(value = case_when(
    value != 0 ~ 1,
    .default = value)) %>%
  group_by(name) %>%
  summarise(total = sum(value)) %>%
  ungroup() %>%
  select(coef = name, nspp = total) %>%
  rbind(data.frame(coef = "(Intercept)", nspp = 675))


# ----  Prepare and plot results ----

numeric_covs <- c("log mass", "clim vuln", "(Intercept)",
                  "latitude", "yearXlatitude", 
                  "latitude:yearXlatitude")

# Prepare results
confints_df <- data.frame(
  coef = names(coef(mavg)),
  estimate = coef(mavg),
  lwr = confint(mavg, level = 0.85)[,1],
  upr = confint(mavg, level = 0.85)[,2]) %>%
  filter(!str_detect(coef, "s\\(tree_sci_name\\)")) %>%
  filter(!str_detect(coef, "s\\(bcr\\)")) %>%
  filter(!str_detect(coef, "s\\(species\\)")) %>%
  # filter(coef != "(Intercept)") %>%
  `row.names<-`(NULL) %>%
  as_tibble() %>%
  left_join(cov_npp, by = join_by(coef)) %>%
  mutate(coef = str_replace_all(coef, "_", " ")) %>%
  mutate(shape = case_when(
    coef %in% numeric_covs ~ 21,
    .default = 22)) %>%
  mutate(coef_name = case_when(
    coef == "log mass" ~ "log(Mass)",
    coef == "urban" ~ "Urban associates",
    coef == "(Intercept)" ~ "Intercept",
    coef == "agriculture" ~ "Agriculture associates",
    coef == "clim vuln" ~ "Climate vulnerability",
    coef == "yearXlatitude" ~ "Shifted north",
    coef == "latitude" ~ "Northern ranged",
    coef == "latitude:yearXlatitude" ~ "Shifted out of region",
    coef == "Initial abundance" ~ "log(Initial abundance)",
    coef == "Urban adapted species" ~ "Urban-adapted species",
    .default = coef)) %>%
  mutate(color = case_when(
    ((lwr < 0) & (upr < 0)) ~ "#c9c1fa",
    ((lwr > 0) & (upr > 0)) ~ "#4734bc",
    .default = "#785ef0")) %>%
  mutate(nspp = as.character(nspp)) %>%
  mutate(nspp = if_else(
    estimate == min(estimate),
    paste("n =", nspp),
    nspp))

# Plot
ggplot(confints_df, aes(x = estimate, y = reorder(coef_name, -estimate),
                        color = color, fill = color, shape = shape)) +
  geom_hline(yintercept = c(seq(1,20, by = 2)), color = "white") +
  geom_hline(yintercept = c(seq(1,20, by = 2)), color = "gray70", linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, linewidth = 1) +
  geom_point(size = 5, color = 1, stroke = 0.6) +
  geom_label(aes(x = 0.018, label = nspp), hjust = 1, size = 4.5,
             color = "white", fill = "white") +
  geom_text(aes(x = 0.018, label = nspp),  size = 4.5,
             hjust = 1, color = 8) +
  scale_shape_identity(guide = "legend", breaks = c(21,22),
                       name = NULL,
                       labels = c("Numeric", "Binary")) +
  scale_color_identity(guide = "legend", name = NULL,
                       breaks = c("#c9c1fa", "#785ef0","#4734bc"),
                       labels = c("Negative", "Insig.", "Positive")) +
  scale_fill_identity() +
  guides(color = guide_legend(
    order = 1,
    override.aes = list(linewidth = 3))) +
  labs(x = "Model-averaged coefficient estimate", y = NULL) +
  theme_light(16) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    axis.text = element_text(family = "Roboto", color = 1),
    axis.title = element_text(family = "Roboto", color = 1),
    panel.grid.major.x = element_blank(),
    axis.ticks.y = element_blank(),
    aspect.ratio = 1,
    panel.grid.minor = element_blank()
  )


# Interpreting range shifts covs

df_bcr[[1]] %>%
  select(species, latitude, yearXlatitude) %>%
  mutate(intx = (scale(latitude) * scale(yearXlatitude))[,1]) %>%
  arrange((intx))

# positive yearXlatitude is shifting north
# negative yearXlatitude is shifting south

# positive latitude is more abundant in the north
# negative latitude is more abundant in the south

# positive interaction is species shifting out of the bcr
# negative interaction is species shifting into the bcr
