# SIAR model
# library(pacman)
# pacman::p_load(magrittr, tidyverse, plotly, RColorBrewer,
#                magick)

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
n_cities <- 90
days_infectious <- 5 # days spent infectious on average
days_asymptomatic <- 2 # days spend in asymptomatic stage
start_infectious <- 3
beta_I <-0.2
beta_A <-0.07
fraction_who_travel <- 0.01
asymptomatic_ratio <- 3 # corresponds to 75% asymptomatic

set.seed(18)
cities_population_vec <- round(rlnorm(n_cities, meanlog = 6.2, sdlog = 2.3))
cities_population_vec <- cities_population_vec + 1000
set.seed(10)
city_locations <- data.frame(city_id = 1:n_cities,
                             x_val = rnorm(n_cities, 724354, sd = 10000),
                             y_val = rnorm(n_cities, 6175803, sd = 10000)
)


# MAKE CITIES --------------------------------------------------------

city <- makeCity(n_people, n_time = n_days,
                 seed_infectious = start_infectious,
                 n_cities = n_cities,
                 n_ppl = cities_population_vec)

city <- seedOneCity(x = city,
                    start_infectious = start_infectious)




# RUN WITHOUT TRAVEL ------------------------------------------------------
# # 
# itr <- 200
# out_raw <- map(1:itr, function(x){
#   loopOverDays(x = city, n_steps = n_days, n_cities = n_cities,
#                frac_travel = fraction_who_travel,
#                beta_I = beta_I,
#                beta_A = beta_A,
#                decay = 4,
#                travel = FALSE)
# } )
# 
# out <- byDayAndItr(out_raw, itr, n_days, n_cities=)
# 
# out2 <- summarizeCityI(out, n_cities = 1)
# 
# # Variance on within-city model
# var_within_city <- data.frame(var = varMonitor(out2$cumAR, itr))
# 
# varPlot(var_within_city,"Variance of cumulative AR") %>%
#   ggsave(filename = "out/var_within.jpg", plot = .,  device = "jpg")
# 

# CONVERGENCE CHECKING ------------------------------------------------------
# 
itr <- 400
out_raw <- iterateSim(x = city, n_steps = n_days, n_cities = n_cities,
                      frac_travel = fraction_who_travel,
                      beta_I = beta_I,
                      beta_A = beta_A,
                      decay = 4,
                      travel = TRUE,
                      itr = itr)

out <- byDayAndItr(out_raw, itr, n_days, n_cities=n_cities)
cum_I <- summarizeCityI(out, n_cities = n_cities)
prop_infected <- propInfectPerItr(cum_I)



# Check variance
var_prop_infected <- varMonitor(prop_infected$prop_infected, itr)
varPlot(var_prop_infected, "Variance on proportion of cities infected ") %>%
  ggsave(filename = "out/var_propotion_SAIR.jpg", plot = .,  device = "jpg")

# Need >250 runs for convergence in the smaller towns



# SENSITIVITY: TRAVEL -----------------------------------------------------
travel_vec <- seq(from = 0.005, to = 0.05, length.out = 10)
itr <- 300
set.seed(13)
travel_sensitivity_raw <- map(seq_along(travel_vec), function(i){
  # browser()
  iterateSim(x = city, n_steps = n_days, n_cities = n_cities,
             frac_travel = travel_vec[i],
             beta_I = beta_I,
             beta_A = beta_A,
             decay = 4,
             travel = TRUE,
             itr = itr)
})

travel_sens <- propAllTravel(travel_sensitivity_raw,
                             n_cities = n_cities,
                             itr = itr) %>%
  extractMeanAndCI(travel_vec)

c <- plotSensitivity(travel_sens, "Fraction who travel")
%>%
  ggsave(filename = "out/frac_travel.jpg", plot = .)



# SENSITIVITY BETA_A --------------------------------------------------------
beta_A_vec <- seq(from = 0, to = 0.1, length.out = 10)
itr <- 250

beta_sensitivity_raw <- map(seq_along(beta_A_vec), function(i){
  # browser()
  iterateSim(x = city, n_steps = n_days, n_cities = n_cities,
             frac_travel = fraction_who_travel,
             beta_I = beta_I,
             beta_A = beta_A_vec[i],
             decay = 4,
             travel = TRUE,
             itr = itr)
})


beta_sens <- propAllTravel(beta_sensitivity_raw,
                           n_cities = n_cities,
                           itr = itr) %>%
  extractMeanAndCI(beta_A_vec)

b <- plotSensitivity(beta_sens, "Beta A value")
%>%
  ggsave(filename = "out/beta_sens.jpg", plot = .)



# SENSITIVITY TO DECAY ----------------------------------------------------
decay_vec <- seq(from = 1, to = 4, by = 1)
itr <- 250

decay_sensitivity_raw <- map(seq_along(decay_vec), function(i){
  # browser()
  iterateSim(x = city, n_steps = n_days, n_cities = n_cities,
             frac_travel = fraction_who_travel,
             beta_I = beta_I,
             beta_A = beta_A,
             decay = decay_vec[i],
             travel = TRUE,
             itr = itr)
})


decay_sens <- propAllTravel(decay_sensitivity_raw,
                            n_cities = n_cities,
                            itr = itr) %>%
  extractMeanAndCI(decay_vec)

a <- plotSensitivity(decay_sens, "Decay value")

  ggsave(filename = "out/decay_sens.jpg", plot = .)



library(cowplot)

cowplot::plot_grid(a,b,c, labels = ("AUTO")) %>%
  ggsave(filename = "out/sen.jpg", plot = .)



# MAPS OF SINGLE OUTBREAKS BY DAY -----------------------------------------
itr = 1
set.seed(14)
n_days <- 120
city <- makeCity(n_people, n_time = n_days,
                 seed_infectious = start_infectious,
                 n_cities = n_cities,
                 n_ppl = cities_population_vec) %>%
  seedOneCity(start_infectious = start_infectious)

fixed_parameters_out <- iterateSim(x = city,
                                   n_steps = n_days,
                                   n_cities = n_cities,
                                   frac_travel = fraction_who_travel,
                                   beta_I = beta_I,
                                   beta_A = beta_A,
                                   decay = 3,
                                   travel = TRUE,
                                   itr = itr)

splitByDay <- function(x){
  # Function to split dataset into a list - each element = one day
  split(x[[1]], f = x[[1]]$day)
}

by_day <- splitByDay(fixed_parameters_out)

map_per_day <- lapply(by_day, mapPrepEachDay, utm = city_locations)


mapCities(map_per_day[[65]])

# Step 1: create png for each day
# Map of cities
seq(from = 1, to = n_days, by=1) %>% 
  map_df(mapCitiesByDay, data_map = map_per_day)


# # Step 2: List those Plots, Read them in, and then make animation
# out <- list.files(path = "gif", pattern = "*.png", full.names = T) %>% 
#   file.info()
#   map(image_read) # reads each path file
#   
# 
# image_join(out) %>% # joins image
#   image_animate(fps=2) %>% # animates, can opt for number of loops
#   image_write("infection_spread.gif") # write to current dir
# 


# MAP2 --------------------------------------------------------------------


itr = 1
set.seed(13)
fixed_parameters_out <- iterateSim(x = city,
                                   n_steps = n_days,
                                   n_cities = n_cities,
                                   frac_travel = fraction_who_travel,
                                   beta_I = beta_I,
                                   beta_A = beta_A,
                                   decay = 3,
                                   travel = TRUE,
                                   itr = itr)

map_tmp <- byDayAndItr(fixed_parameters_out, itr = itr, n_days = n_days,
                       n_cities = n_cities) %>%
  summarizeCityI() %>%
  summarizeOverAllItr()

map_out <- merge(meanAR, city_locations, by = "city_id")
map_out$tot_N <- cities_population_vec



# PROFILING CODE ----------------------------------------------------------
devtools::install_github("hadley/lineprof")
library("lineprof")
l <- lineprof(map(1:itr, function(x){
  loopOverDays(x = city, n_steps = n_days, n_cities = n_cities,
               frac_travel = fraction_who_travel)
} ))
shine(l)

l <- lineprof(cum_I <- summarizeCityI(out))
shine(l)

ci <- merge(city, city_locations, by="city_id")
ci <- ci[!duplicated(ci$city_id), ]
city_map <- ggplot(ci) +
  geom_point(aes(x = x_val, y= y_val,
                 size = tot_N),
             alpha = 0.6) +
  theme_classic()
ggsave(filename = "out/city_map.jpg", plot = city_map)



# plot over time ----------------------------------------------------------
itr = 10
out_raw <- iterateSim(x = city,
           n_steps = n_days,
           n_cities = n_cities,
           frac_travel = fraction_who_travel,
           beta_I = beta_I,
           beta_A = beta_A,
           decay = 3,
           travel = TRUE,
           itr = itr)

out <- prettyOutput(out_raw, itr = itr, n_days = n_days, n_cities = n_cities)

plotCitiesOverTime(out, alpha = .5)
