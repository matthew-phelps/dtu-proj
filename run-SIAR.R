# SIAR model
# library(pacman)
# pacman::p_load(magrittr, tidyverse, purr, plotly)
# pacman::p_loaded() # Check loaded packages
rm(list = ls())

source("SIAR-functions.R")
# Each day is composed of two steps:
# 1. Infect and update within each city due to internal transmission
# 2. Move people between cities and update each city's population


# INITIALIZE STARTING VALUES -----------------------------------------------
# Create some dummy data for testing. Once the model functions are working we
# will separate the functions, and call the functions with the data in a
# separate script.

n_days <- 200 # time frame (in days) of model
n_cities <- 100
days_infectious <- 5 # days spent infectious on average
days_asymptomatic <- 2 # days spend in asymptomatic stage
start_infectious <- 3
beta_I <-0.2
beta_A <-0.07
fraction_who_travel <- 0.05
asymptomatic_ratio <- 3 # corresponds to 75% asymptomatic

set.seed(19)
cities_population_vec <- round(rlnorm(100, meanlog = 6.2, sdlog = 1.9))

# MAKE CITIES --------------------------------------------------------

city <- makeCity(n_people, n_time = n_days,
                 seed_infectious = start_infectious,
                 n_cities = n_cities,
                 n_ppl = cities_population_vec)

city <- seedOneCity(x = city,
                    start_infectious = start_infectious)



# RUN MODEL ----------------------------------------------------------------
# set.seed(13)
# loopOverDays(x = city, n_steps = n_days, n_cities = n_cities,
#              frac_travel = fraction_who_travel,
#              days_asymptomatic)




itr <- 12
out_raw <- map(1:itr, function(x){
  loopOverDays(x = city, n_steps = n_days, n_cities = n_cities,
               frac_travel = fraction_who_travel)
} )

out <- prettyOutput(out_raw, itr, n_days, n_cities)
# out[which.max(out$I), ]

# ggploting is sloooow with large number of lines. May take ~30 secs.
plot1 <- plotCitiesOverTime(out, alpha = .5)
ggplotly(plot1)

#
cum_I <- cumInfected(out_raw)
hist(cum_I)

city_locations <- data.frame(city_id = 1:n_cities,
           x_val = rnorm(n_cities, 724354, sd = 1000),
           y_val = rnorm(n_cities, 6175803, sd = 1000)
)


plot(city_locations$x_val, city_locations$y_val)
