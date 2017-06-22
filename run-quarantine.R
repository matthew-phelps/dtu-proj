# Within-city infection
library(pacman)
pacman::p_load(magrittr, tidyverse, purr)
pacman::p_loaded() # Check loaded packages
library(ggplot2)
rm(list = ls())
source("quarantine-functions.R")
# Each day is composed of two steps:
# 1. Infect and update within each city due to internal transmission
# 2. Move people between cities and update each city's population

# This script will run both steps inside "loopOverDays()" 



# INITIALIZE STARTING VALUES -----------------------------------------------
# Create some dummy data for testing. Once the model functions are working we
# will separate the functions, and call the functions with the data in a
# separate script.

n_days <- 100 # time frame (in days) of model
n_cities <- 10
days_infectious <- 5 # days spent infectious on average
start_infectious <- 45
beta <-0.2
fraction_who_travel <- 0.05

set.seed(13)
cities_population_vec <- sample(1000:30000, n_cities, replace = TRUE)



# HELPER FUNCTIONS --------------------------------------------------------

city <- makeCity(n_people, n_time = n_days,
                 seed_infectious = start_infectious,
                 n_cities = n_cities,
                 n_ppl = cities_population_vec)

# movement <- oneStepMovingMatrix(n.cities = n_cities,
#                                 max.people = 5)
# infectionByMovement(movement, city)
# RUN MODEL ----------------------------------------------------------------
set.seed(13)
# out <- loopOverDays(x = city, n_steps = n_days, n_cities = n_cities,
#              frac_travel = fraction_who_travel, start.quarantine = 75, end.quarantine = 10)
# 
# out
# sum(out$in.quara)

itr <- 2
out_raw <- map(1:itr, function(x){
  loopOverDays(x = city, n_steps = n_days, n_cities = n_cities,
               frac_travel = fraction_who_travel, start.quarantine = 75, end.quarantine = 10)
  } 
)

out <- prettyOutput(out_raw, itr, n_days, n_cities)
head(out)
plot(out$I[out$city_id==1], type ="l")



interaction()
# ggploting is sloooow with large number of lines. May take ~30 secs.
plotCitiesOverTime(out, alpha = .1)
plots <- plotMeanInfection(out)
multiplot(plots)



cum_I <- cumInfected(out_raw)
hist(cum_I)
