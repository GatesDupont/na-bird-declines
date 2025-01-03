# Function to attempt fitting the model with negative binomial, fallback to Poisson on warning
fit_model_with_fallback <- function(data) {
  tryCatch({
    # Try fitting with negative binomial
    model <- bam(species_total ~ year + s(year_factor, bs = "re") +
                   s(observer, bs = "re") + s(route, bs = "re") + 
                   s(duration, k = 4) + s(yday, k = 4) + first_year,
                 data = data,
                 family = nb(),
                 control = gam.control(trace = FALSE),
                 gamma = 1.4,
                 discrete = TRUE)
    # message("\n Model fitted with Negative Binomial distribution.")
    return(model)
  }, warning = function(w) {
    # If a warning is caught, fit with Poisson
    # message("\n Warning caught; fitting with Poisson distribution instead.")
    model <- bam(species_total ~ year + s(year_factor, bs = "re") +
                   s(observer, bs = "re") + s(route, bs = "re") + 
                   s(duration, k = 4) + s(yday, k = 4) + first_year,
                 data = data,
                 family = poisson(),
                 control = gam.control(trace = FALSE),
                 gamma = 1.4,
                 discrete = TRUE)
    return(model)
  })
}