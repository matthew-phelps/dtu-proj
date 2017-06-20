# This file contains all the functions to run the model for SIAR model. It is
# called by the "run-SIR-within.R" script. There are 3 sections of functions -
# SIMULATION, DATA WRANGLING, and MISCELLANEOUS


# SIMULATION --------------------------------------------------------------
# These functions are the core of the within-city transmission simulation

loopOverDays <- function(x = city,
                         n_steps = n_days,
                         n_cities = n_cities,
                         frac_travel,
                         ...){
  # Function to loop over all time-steps.
  # Two steps happen during each iteration of this loop:
  # 1. Internal transmission and updates
  # 2. Move people between cities
  
  # Function is vectorized with respected to multiple cities - i.e. internal
  # transmission in all cities will be updated at the same time in this
  # function.
  contact_matrix <- makeCityContactMatrix(n_cities)
  # browser()
  for(i in 1:n_steps){
    
    
    # Step 1 - updates from internal transmission for day i
    step1_output  <- simOneTimeStep(x = x[x$day==i, ],
                                    time_in_A = days_asymptomatic)
    # browser()
    
    
    # Step 2 - Move Infected
    if(any(step1_output$num_I >0)){
      # If there are any infected -> move them between cities
      step2_output <- moveAsymptomatic(step1_output = step1_output,
                                       n_cities = n_cities,
                                       contact_matrix,
                                       frac_travel = frac_travel)
      # Increment by 1 timse-step
      x[x$day == (i + 1), ] <- step2_output
    } else {
      x[x$day == (i + 1), ] <- step1_output
    }
    
  }
  
  return(x)
}


simOneTimeStep<- function(x,
                          time_in_I = days_infectious,
                          time_in_A = days_asymptomatic){
  browser()
  # In this function we infect people from within-city transmission.
  
  # If there are no Infected in all cities, nothing can happen so we
  # simply increment the day by 1 and end the loop.
  inx <- all(x$num_I==0 & x$num_A==0)
  if(inx){
    x$day <- x$day + 1
    return(x)
  }
  # browser()
  # Total population current time-step for each city.
  N <- rowSums(x[, c("num_S", "num_I", "num_A", "num_R")]) 
  x$new_I <- numberOfNewI(beta_I = beta_I,
                          beta_A = beta_A,
                          I = x$num_I,
                          A = x$num_A,
                          x$num_S,
                          N = N) # New Infected
  newA <- numberOfNewA(x = x)
  newR_from_I <- numberOfNewR(x$num_I, time_in_I=time_in_I)  # New Recovered
  newR_from_A <- numberOfNewRFromA(x$num_A, time_in_A=time_in_A)
  
  # If any new Infected OR Recovered, update city object.
  inx <- any(x$newI > 0) | any(newR_from_I > 0) | any(newR_from_A > 0)
  if(inx){
    x$num_R <- x$num_R + newR_from_I + newR_from_A
    x$num_I <- x$num_I + x$new_I - newR_from_I
    x$num_A <- x$num_A + newA - newR_from_A
    x$num_S <- N - rowSums(x[, c("num_I", "num_A", "num_R")])
  }
  
  # Increment day by 1
  x$day <- x$day + 1
  return(x)
}

numberOfNewI <- function(beta_I, beta_A,
                         I, A, S, N){
  # The number of infections on time t + 1 is determined by a poisson
  # distribution with a mean(lambda) defined below
  
  # I infected by I
  lambda_I <- (beta_I*I) * (S / N)
  out_I <- rpois(n_cities, lambda_I)
  
  # I infected by A
  lambda_A <- (beta_A*A) * (S / N)
  out_A <- rpois(n_cities, lambda_A)
  
  # browser()
  out <- out_I + out_A
  # if(any(is.na(lambda_I))) browser()
  
  # If we select more I than there are S -> set I to equal S
  inx <- out > S
  if(any(inx)){
    # browser()
    out[inx] <- S[inx]
  }
  
  return(out)
}

numberOfNewA <- function(x){

  A <- 3 * x$new_I
  max_A <- x$num_S - x$new_I
  inx <- A > max_A
  if(any(inx)){
    # browser()
    A[inx] <- max_A[inx]
  }
  A
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

numberOfNewRFromA <- function(A, time_in_A){
  # Recovery is poisson distributed with a mean(lambda) = (number of infected) *
  # (rate of recovery). Rate of recovery = 1 / (time spent in infectious state).
  out_R <- rpois(n_cities, A * 1/time_in_A)
  
  # If we randomly select more R than there are I -> set R to equal I
  inx <- out_R > A
  if(any(inx)){
    # browser()
    out_R[inx] <- A[inx]
  }
  
  return(out_R)
}




# MOVING PEOPLE -----------------------------------------------------------

countMovedAsymptomatic <- function(num_asymptomatic,
                                   n_cities,
                                   frac_travel){
  out <- rpois(n_cities, frac_travel*num_asymptomatic) 
  inx <- out > num_asymptomatic
  if(any(inx)){
    # browser()
    out[inx] <- num_asymptomatic[inx]
  }
  
  return(out)
}

moveAsymptomatic <- function(step1_output,
                             n_cities,
                             contact_matrix,
                             frac_travel = frac_travel){
  # browser()
  # Enumerate number of infected who move FROM each city
  tmp <- step1_output
  tmp$n_move_out <- countMovedAsymptomatic(num_asymptomatic = tmp$num_A,
                                           n_cities,
                                           frac_travel = frac_travel)
  
  # Remove infected out of SENDING cities
  tmp$num_A <- tmp$num_A - tmp$n_move_out
  
  # Loop through cities that send at least 1 infected and define where
  # infected move TO
  inx <- tmp[, "n_move_out"] == 0
  move_counter <- list()
  move_counter <- lapply(1:nrow(tmp), function(k){
    # browser()

    return(sample(1:n_cities, tmp$n_move_out[k],
                  replace = TRUE,
                  prob = contact_matrix[, k]))
  }) %>%
    unlist()
  
  # Summarize how many times a city receives an Infected person
  new_I_each_city <- data.frame(table(move_counter))
  
  # Add Infected to RECEIVING cities 
  # browser()
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
  
  out$A <- lapply(x, "[[", "num_A") %>%
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
                      num_S = n_ppl[i],
                      num_I = 0,
                      num_A = 0,
                      num_R = 0,
                      tot_N = n_ppl[i],
                      new_I = 0,
                      n_move_out = NA)
    out[2:n_time, c("num_S", "num_I", "num_A", "num_R", "new_I")] <- NA
    out
  })
  out <- do.call(rbind.data.frame, out)
}


seedOneCity <- function(x = city,
                        start_infectious = start_infectious){
  # This function seeds the largest city with "start_infectious" number of
  # infections + 3x that many asymptomatic
  
  seed_city <- x[which.max(x$tot_N), "city_id"]
  inx <- x$city_id == seed_city & x$day == 1
  x[inx, c("num_I", "num_A")] <- c(start_infectious, (3 * start_infectious))
  x
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
  city_contact_matrix <- sweep(city_contact_matrix,
                               2,
                               colSums(city_contact_matrix),
                               FUN="/")
  return(city_contact_matrix)
}

plotCitiesOverTime <- function(out, alpha = alpha) {
  ggplot() +
    geom_line(data = out,
              aes(x = day, y = I,
                  group = interaction(city_id, itr),
                  col = as.factor(city_id)),
              alpha = alpha) +
    theme_minimal() +
    theme(legend.position = "none")
}

