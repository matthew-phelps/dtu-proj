### Written by Emil Tosti 
## emtosti@gmail.com
## version .01
# 20.06.17

######------- use the day column

###### Quarantine of City function ######
#The following code will

##
# look at quarantine list
# if on list
#   Determine if the quarantine is lifted
###
# else
# determine the number of infected within a city
# if the city this is above a threshold the city will be put under quarantine
# place on quarantine list
# determine the "bleeding" form the city
#   Calculate the #people moving out
#     Determine if they are sick or not
# determine the "bleeding" into the city
#   Calculate the #people moving in
#     Determine if they are sick or not
# 
#
###
# make other cities more "strick" towards people from quarantine cities
#   Determine the "bleeding" into a city from an quarantine city
#     Calculate if people getting in are infected based on prev. "sick or not"
###


makeQuarantineList <-function(n.cities = n.cities){
  quarantine <- data.frame(city.id = 1:n.cities, 
                           quara = rep(x = FALSE, n.cities))
}

determineInfectFreq <- function(x, inx.quara){
  #number of infected in a given city
  temp.store.infec  <- sum(x$num_I[x$city_id==inx.quara],x$new_I[x$id==inx.quara])
  #all people in a city
  temp.store.n.city <- sum(temp.store.infec, x$tot_N)
  #the freq of infection
  infec.freq <- temp.store.infec/temp.store.n.city
  return(infec.freq)
}

bleedingOfCity <- function(city, quarantine, inx.quara = inx.quara, security.level = security.level){
  n.people <- city$tot_N
  ### use an appropriate selection
  
}


putCityUnderQuarantine <- function(x = city,start.threshold = 50, end.threshold = 20){
  #Determine which cities are under quarantine
  #browser()
  
  inx.quara <- which(x$in.quara == 1)
  if(length(inx.quara)>0){
    for(i in 1:length(inx.quara)){
      ##calling the freq of infection within the city
      #infec.freq <- determineInfectFreq(x, inx.quara[i])
      #evaluation of the quarantine should be lifted
      infect.in.city <- sum(x$num_I[x$city_id == inx.quara[i]],x$new_I[x$city_id == inx.quara[i]])
      if(infect.in.city < end.threshold){
        x$in.quara[x$city_id == inx.quara[i]] <- FALSE
      }
    }
  }
  #evaluate freq of infection in all, non quarantined, cites
  inx.non.quara <- which(!x$in.quara)
  if(length(inx.non.quara)>0){
  #look at which cities have a high enough level of infection to be put under quarantine
    for(i in 1:length(inx.non.quara)){
      #browser()
      # <- determineInfectFreq(x, inx.non.quara[i])
      infect.in.city <- sum(x$num_I[x$city_id == inx.non.quara[i]],x$new_I[x$city_id == inx.non.quara[i]])
      if(infect.in.city >= start.threshold){
        x$in.quara[x$city_id == inx.non.quara[i]] <- TRUE
      }
    }
  }
  return(x)
}


# This file contains all the functions to run the model for the SIR model of
# infection. It is called by the "run-SIR-within.R" script. There are 4 sections
# of functions - SIMULATION, DATA WRANGLING, MOVING PEOPLE and MISCELLANEOUS


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
  # browser()
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
                         n_cities = n_cities,
                         frac_travel){
  # Function to loop over all time-steps.
  # Two steps happen during each iteration of this loop:
  # 1. Internal transmission and updates
  # 2. Evalute cities for quarantine
  # 3. Move people between cities
  
  # Function is vectorized with respected to multiple cities - i.e. internal
  # transmission in all cities will be updated at the same time in this
  # function.
  
  for(i in 1:n_steps){
    # Step 1 - updates from internal transmission for day i
    step1_output  <- simOneTimeStep(x = x[x$day==i, ])
    
    #check if any cities needs to go under quarantine
    #browser()
    step1_output <- putCityUnderQuarantine(step1_output)
    
    contact_matrix <- makeCityContactMatrix(n_cities, cities = step1_output)
    
    # Step 2 - Move Infected
    if(any(step1_output$num_I >0)){
      # If there are any infected -> move them between cities
      step2_output <- moveInfected(step1_output = step1_output,
                                   n_cities = n_cities,
                                   contact_matrix,
                                   frac_travel = frac_travel)
      # Increment by 1 time-step
    
      x[x$day == (i + 1), ] <- step2_output
    } else {
      x[x$day == (i + 1), ] <- step1_output
    }

  }
  return(x)
}


# MOVING PEOPLE -----------------------------------------------------------

countMovedInfected <- function(num_infec, n_cities, frac_travel) {
  out <- rpois(n_cities, frac_travel*num_infec) 
  inx <- out > num_infec
  if(any(inx)){
    # browser()
    out[inx] <- num_infec[inx]
  }
  
  return(out)
}

moveInfected <- function(step1_output, n_cities, contact_matrix,
                         frac_travel = frac_travel){
  # browser()
  # Enumerate number of infected who move FROM each city
  tmp <- step1_output
  tmp$n_move_out <- countMovedInfected(num_infec = tmp$num_I, n_cities,
                                       frac_travel = frac_travel)
  
  # Remove infected out of SENDING cities
  tmp$num_I <- tmp$num_I - tmp$n_move_out
  
  # Loop through cities that send at least 1 infected and define where
  # infected move TO
  inx <- tmp[, "n_move_out"] == 0
  move_counter <- list()
  move_counter <- lapply(1:nrow(tmp), function(k){
    return(sample(1:n_cities, tmp$n_move_out[k],
                  replace = TRUE,
                  prob = contact_matrix[, k]))
  }) %>%
    unlist()
  
  # Summarize how many times a city receives an Infected person
  new_I_each_city <- data.frame(table(move_counter))
  #browser()
  # Add Infected to RECEIVING cities 
  inx <- tmp$city_id %in% new_I_each_city$move_counter
  tmp$num_I[inx] <- new_I_each_city$Freq + tmp$num_I[inx]
  # stopifnot(sum(tmp$n_move_out) == sum(tmp2$Freq, na.rm = TRUE))
  
  return(tmp)
}




# DATA WRANGLING ----------------------------------------------------------

# These functions deal with the data output of the simulations, including
# reshaping and summarizing.

prettyOutput <- function(x, itr, n_days, n_cities){
  # Take output of simulations and make tidy data for ggplot.
  
  out <- lapply(x, "[[", "num_S") %>%
    unlist() %>%
    data.frame(S = .)
  
  out$I <- lapply(x, "[[", "num_I") %>%
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
    out <- data.frame(city_id = i,
                      day = 1:n_time,
                      num_S = n_ppl[i] - seed_infectious,
                      num_I = seed_infectious,
                      new_I = 0,
                      num_R = 0,
                      tot_N = n_ppl[i],
                      
                      n_move_out = NA,
                      in.quara = rep(FALSE, n_time))
    out[2:n_time, c("num_S", "num_I", "new_I", "num_R")] <- NA
    out
  })
  out <- do.call(rbind.data.frame, out)
}

makeCityContactMatrix <- function(n_cities, cities = x){
  # This function creates a weight matrix describing the contact between all
  # city-pairs. Matrix is symmetric since we are assuming a symmetric flow of
  # people b/w cities. The columns represent sending weights, and the rows
  # represent receiving weights.
  
  made_up_data <- rpois(n_cities*n_cities, lambda =  10 )
  city_contact_matrix <- matrix(data = made_up_data, n_cities, n_cities)
  
  # Set diagnols to 0 since city will not send to itself
  diag(city_contact_matrix) <- 0
  
  # Normalize matrix column-wise so each column adds to 1.0
  city_contact_matrix <- sweep(city_contact_matrix,
                               2,
                               colSums(city_contact_matrix),
                               FUN="/")
  #check for quarantined cities
  #browser()
  inx.q <- which(cities$in.quara == 1)
  if(length(inx.q)>0 ){
    city_contact_matrix[inx.q, ] = city_contact_matrix[inx.q, ] * 0.1
  }
  return(city_contact_matrix)
}

plotCitiesOverTime <- function(out, alpha = alpha) {
  ggplot() +
    geom_line(data = out,
              aes(x = day, y = I, group = interaction(city_id, itr)),
              alpha = alpha,
              col = "darkred") +
    theme_minimal()
}







