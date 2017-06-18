# Within-city infection
library(pacman)
pacman::p_load(magrittr, tidyverse, purr)
pacman::p_loaded() # Check loaded packages
rm(list = ls())
source("SIR-functions.R")
# Each day is composed of two steps:
# 1. Infect and update within each city due to internal transmission
# 2. Move people between cities and update each city's population

# This script is for step 1.



# INITIALIZE STARTING VALUES -----------------------------------------------
# Create some dummy data for testing. Once the model functions are working we
# will separate the functions, and call the functions with the data in a
# separate script.

n_days <- 100 # time frame (in days) of model
n_people <- 1000
n_cities <- 1000
beta <-0.2
days_infectious <- 5 # days spent infectious on average

start_infectious <- 5

# HELPER FUNCTIONS --------------------------------------------------------

city <- makeCity(n_people, n_time = n_days,
                 seed_infectious = start_infectious,
                 n_cities = n_cities)

# RUN MODEL ----------------------------------------------------------------
itr <- 12
out_raw <- map(1:itr, function(x){
  loopOverDays(x = city, n_steps = n_days, n_cities = n_cities)
} )

out <- prettyOutput(out_raw, itr, n_days, n_cities)
out[which.max(out$I), ]

# ggploting is sloooow with large number of lines. May take ~30 secs.
plotCitiesOverTime(out)

#
cum_I <- cumInfected(out_raw)
hist(cum_I)