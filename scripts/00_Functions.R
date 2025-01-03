library(tidyverse)
library(jsonlite)
library(vroom)


# ---- Functions ----

# Get taxnomy CSV file from eBird API
get_taxonomy <- function(){
  
  url <- "https://api.ebird.org/v2/ref/taxonomy/ebird"
  
  data <- data.table::fread(url) %>%
    select(`Common Name` = COMMON_NAME, 
           species_code = SPECIES_CODE)
  
}

# Get taxnomy CSV file from eBird API
get_sci_names <- function(){
  
  # url <- "https://api.ebird.org/v2/ref/taxonomy/ebird"
  url = "https://api.ebird.org/v2/ref/taxonomy/ebird?&version=2023"
  
  data <- data.table::fread(url) %>%
    select(`Common Name` = COMMON_NAME, 
           sci_name = SCIENTIFIC_NAME)
  
}

check_species <- function(species, species_list, combine_species_forms,
                          quiet = FALSE) {
  nm <- deparse(substitute(species))
  if(missing(species) || is.null(species)) {
    stop("No `", nm, "` specified", call. = FALSE)
  }
  
  species <- tolower(species)
  s_list_combo <- dplyr::filter(species_list, .data$unid_combined == TRUE)
  s_list_no_combo <- dplyr::filter(species_list, .data$unid_combined == FALSE)
  
  s_combo <- dplyr::filter(
    s_list_combo,
    .env$species == tolower(.data$english) |
      .env$species == tolower(.data$french) |
      .env$species == tolower(paste(.data$genus, .data$species)) |
      .env$species == .data$aou)
  
  s_no_combo <- dplyr::filter(
    s_list_no_combo,
    .env$species == tolower(.data$english) |
      .env$species == tolower(.data$french) |
      .env$species == tolower(paste(.data$genus, .data$species)) |
      .env$species == .data$aou)
  
  
  if(nrow(s_combo) == 0 & nrow(s_no_combo) == 0) {
    stop("Invalid species. Ensure `", nm, "` matches a ",
         "English or French species name or AOU code.\n",
         "See `search_species()` for a generic search.",
         call. = FALSE)
  } else if(nrow(s_no_combo) == 0 &
            !combine_species_forms & nrow(s_combo) == 1){
    stop("`combine_species_forms = FALSE` but '",
         s_combo$english, "' is a combined form...\nif you want this species, ",
         "set `combine_species_forms = TRUE`", call. = FALSE)
  } else if(nrow(s_combo) == 0 & combine_species_forms & nrow(s_no_combo) == 1){
    stop("`combine_species_forms = TRUE` but '",
         s_no_combo$english, "' is an unidentified,\nnon-combined form...",
         "if you want this species, set ",
         "`combine_species_forms = FALSE`", call. = FALSE)
  } else if((combine_species_forms & nrow(s_combo) > 1) ||
            (!combine_species_forms & nrow(s_no_combo) > 1)) {
    stop("Multiple species matched. ",
         "See `search_species()`for a generic search.", call. = FALSE)
  }
  if(combine_species_forms) s <- s_combo else s <- s_no_combo
  
  if(!quiet) message("Filtering to species ", s$english, " (", s$aou, ")")
  s$aou
}
