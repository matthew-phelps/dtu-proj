# This file contains all the functions to run the model for within-city
# infection. It is called by the "run-SIR-within.R" script. There are 3 sections of functions - SIMULATION, DATA WRANGLING, and MISCELLANEOUS


# SIMULATION --------------------------------------------------------------

# These functions are the core of the within-city transmission simulation

numberOfNewI <- function(beta, I, S, N){
  # The number of infections on time t + 1 is determined by a poisson
  # distribution with a mean(lambda) defined below
  lambda <- (beta*I) * (S / N)
  out <- rpois(n_cities, lambda)
  
  # If we randomly select more I than there are S -> set I to equal S
  inx <- out > S
  if(any(inx)){
    browser()
    out[inx] <- S[inx]
  }
  
  return(out)
}

numberOfNewR <- function(I, time_in_I){
  # Recovery is poisson distributed with a mean(lambda) = (number of infected) *
  # (rate of recovery). Rate of recovery = 1 / (time spent in infectious state).
  out_R <- rpois(n_cities, I * 1/time_in_I)
  
  # If we randomly select more R than there are I -> set R to equal I
  inx <- out_R > I
  if(any(inx)){
    # browser()
    out_R[inx] <- I[inx]
  }
  
  return(out_R)
}



simOneTimeStep<- function(x,
                          time_in_I = days_infectious){
  # browser()
  # In this function we infect people from within-city transmission.
  
  # If there are no Infected in all cities, nothing can happen so we
  # simply increment the day by 1 and end the loop.
  if(all(x$num_I==0)){
    x$day <- x$day + 1
    return(x)
  }
  
  # Total population current time-step for each city.
  N <- rowSums(x[, c("num_S", "num_I", "num_R")]) 
  x$new_I <- numberOfNewI(beta, x$num_I, x$num_S, N = N) # New Infected
  newR <- numberOfNewR(x$num_I, time_in_I=time_in_I)  # New Recovered
  
  # If any new Infected OR Recovered, update city object.
  if(any(x$newI > 0) | any(newR > 0)){
    x$num_R <- x$num_R + newR
    x$num_I <- x$num_I + x$new_I - newR
    x$num_S <- N - rowSums(x[, c("num_I", "num_R")])
  }
  
  # Increment day by 1
  x$day <- x$day + 1
  return(x)
}


loopOverDays <- function(x = city,
                         n_steps = n_days,
                         n_cities = n_cities){
  # Function to loop over all time-steps.
  # Two steps happen during each iteration of this loop:
  # 1. Internal transmission and updates
  # 2. Move people between cities
  
  # Function is vectorized with respected to multiple cities - i.e. internal
  # transmission in all cities will be updated at the same time in this
  # function.
  
  for(i in 1:n_steps){
    browser()
    # Subset row of current time-step.
    sub_x <- x[x$day==i, ]
    
    # Step 1 - updates from internal transmission
     tmp  <- simOneTimeStep(x = sub_x)
     
     tmp$n_move <- rpois(n_cities, 0.05*tmp$tot_N)
     
    x[x$day == (i + 1), ]
    oneStepMovingMatrix(n.cities = n_cities, max.people = max.people)
    infectionByMovement(sub_x)
  }
  return(x)
}


oneStepMovingMatrix <- function(n.cities = n_cities, max.people = max.people){
  #Function will create a data.frame that holds infomation about:
  #The sending city
  #The receiving city
  #The number of people moved
  movement <- data.frame(sending = sample(1:n.cities, n.cities, replace = FALSE),
                         receiving = sample(1:n.cities, n.cities, replace = TRUE),
                         n.people = sample(1:max.people, n.cities, replace = TRUE))
  #indexing the cases where the a city is sending people to themself
  inx <- which(movement$sending == movement$receiving)
  if(length(inx)>0){
    # This might be a problem if there are multiple cities seeding itself. What
    # happens if both the first and the last city seed themselves?
    
    #if the city is not the last city, increment the receiving city with one
    if(movement$receiving[inx] < n.cities) movement$receiving[inx] + 1
    #if the city is not the last city, decrese the receiving city with one
    if(movement$receiving[inx] >= 1) movement$receiving[inx] - 1
  }
  return(movement)
}


probAtMovement <- function(prev = prevalence, n.moved = n.moved){
  prob <- 1 - (1 - prev)^n.moved
  return(prob)
}


infectionByMovement <- function(mM = movementMatrix, cities = cities, prev = 0.1){
  browser()
  #index of infectet cities
  infect.id <- cities$city_id[which(cities$num_S>0)]
  for(i in 1:length(infect.id)){
    #Index on reciving town
    inx.move <- mM$sending==infect.id[i]
    #initialzing the temp.storage for the movement
    rec.temp.store <- mM[inx.move,]
    
    rec.temp.store$tot_animal <- ave(rec.temp.store$n.people,
                                     rec.temp.store$receiving, 
                                     FUN = function(x) sum(x))
    #getting rid of duplicates
    rec.temp.store[!(duplicated(rec.temp.store$receiving)), ]
    #calculating the prob. for an infection due to movement
    rec.temp.store$porb <- probAtMovement(prev, rec.temp.store$tot_animal)
    #see if there is a spread of infection
    rec.temp.store$newInfected <- rbinom(length(rec.temp.store$porb), 1, rec.temp.store$porb)==1
  }
  return()
}











# DATA WRANGLING ----------------------------------------------------------

# These functions deal with the data output of the simulations, including
# reshaping and summarizing.

prettyOutput <- function(x, itr, n_days, n_cities){
  # Take output of simulations and make tidy data for ggplot.
  
  out <- lapply(x, "[[", 2) %>%
    unlist() %>%
    data.frame(S = .)
  
  out$I <- lapply(x, "[[", 3) %>%
    unlist()
  out$new_I <- lapply(x, "[[", "new_I") %>%
    unlist()
  
  out$city_id <- lapply(x, "[[", "city_id") %>%
    unlist()
  
  out$itr <- rep(1:itr, each = n_days*n_cities)
  out$day <- rep(1:n_days)
  
  return(out)
}


cumInfected <- function(x){
  # Get the cummulative number of infected per city per iteration.
  lapply(x, function(x){
    # browser()
    x$cum_sum <- ave(x$new_I, x$city_id, FUN = sum)
    x$cum_sum[!duplicated(x$city_id)]
  }) %>%
    unlist()
}



# MISCELLANEOUS FUNCTIONS --------------------------------------------------

# Other helper functions

makeCity <- function(n_people, n_time, seed_infectious,
                     n_cities, n_ppl){
  # browser()
  out<- lapply(1:n_cities, function(i) {
    out <- data.frame(day = 1:n_time,
                      num_S = n_ppl[i] - seed_infectious,
                      num_I = seed_infectious,
                      new_I = 0,
                      num_R = 0,
                      tot_N = n_ppl[i],
                      city_id = i)
    out[2:n_time, c("num_S", "num_I", "new_I", "num_R")] <- NA
    out
  })
  out <- do.call(rbind.data.frame, out)
}

makeCityContactMatrix <- function(n_cities){
  # This function creates a weight matrix describing the contact between all
  # city-pairs. Matrix is symmetric since we are assuming a symmetric flow of
  # people b/w cities. The columns represent sending weights, and the rows
  # represent receiving weights.
  
  made_up_data <- rpois(n_cities*n_cities, lambda =  10 )
  city_contact_matrix <- matrix(data = made_up_data, n_cities, n_cities)
  
  # Set diagnols to 0 since city will not send to itself
  diag(city_contact_matrix) <- 0
  
  # Normalize matrix column-wise so each column adds to 1.0
  city_contact_matrix <- sweep(city_contact_matrix, 2, colSums(city_contact_matrix), FUN="/")
  city_contact_matrix
}

plotCitiesOverTime <- function(out) {
  ggplot() +
    geom_line(data = out,
              aes(x = day, y = I, group = interaction(city_id, itr)),
              alpha = 0.01,
              col = "darkred") +
    theme_minimal()
}