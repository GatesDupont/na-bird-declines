# Changes in NA bird communities

This repository contains the code and analyses for the project. The study investigates long-term bird abundance trends across North American birds, focusing on changes to community structure.

## Contents

- `data/`: Processed data files used in the analysis.
- `scripts/`: R scripts for data preprocessing and modeling.

- **scripts/**: Contains all the R scripts used for data processing and analysis.
  - `00_Functions.R`: Utility functions for analysis.
  - `01_01_Select_species_Atlantic_Northern_Forest_.R` to `01_06_Select_species_Shortgrass_Prairie.R`: Scripts for selecting species in different BCRs.
  - `02_Model_density.R`: Script for density modeling.
  - `03_Drivers_phylo.R`: Script analyzing phylogenetic drivers.
  - `fit_model_with_fallback.R`: Script for fitting a negative binomial model with a Poisson model as a fallback when the dispersion parameter approaches infinity.
  - `theory/`: Contains scripts for theoretical analysis.
    - `theoretical_models.R`: Script for theoretical analysis.
    - `concept_figure.R`: Script for theoretical analysis.
  - `additional_analyses/`: Contains scripts for additional analyses.
