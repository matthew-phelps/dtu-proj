# This file contains all the functions to run the model for SIAR model. It is
# called by the "run-SIR-within.R" script. There are 3 sections of functions -
# SIMULATION, DATA WRANGLING, and MISCELLANEOUS


# SIMULATION --------------------------------------------------------------
# These functions are the core of the within-city transmission simulation


iterateSim <- function(x, n_steps,
                       n_cities,
                       frac_travel,
                       beta_I,
                       beta_A,
                       decay,
                       travel,
                       itr){
  
  # This function simulates the system "itr" number of times
  map(1:itr, function(x){
    # browser()
    loopOverDays(x = city, n_steps = n_steps, n_cities = n_cities,
                 frac_travel = frac_travel,
                 beta_I = beta_I,
                 beta_A = beta_A,
                 decay = decay,
                 travel = travel)
  } )
  
}


loopOverDays <- function(x = city,
                         n_steps = n_steps,
                         n_cities = n_cities,
                         beta_I,
                         beta_A,
                         frac_travel, fake = FALSE,
                         decay = decay,
                         travel){
  # Function to loop over all time-steps.
  # Two steps happen during each iteration of this loop:
  # 1. Internal transmission and updates
  # 2. Move people between cities
  
  # Function is vectorized with respected to multiple cities - i.e. internal
  # transmission in all cities will be updated at the same time in this
  # function.
  # browser()
  if(travel == FALSE){
    # browser()
    id <- x$city_id[x$num_I > 0 & !is.na(x$num_I)]
    x <- x[which(x$city_id==id), ]
    n_cities <- 1
  }
  contact_matrix <- makeCityGravityMatrix(n_cities = n_cities,
                                          dist_mat = dist_mat,
                                          decay = decay)
  if(fake) contact_matrix <- makeCityContactMatrix(n_cities)
  # browser()
  for(i in 1:n_steps){
    
    
    # Step 1 - updates from internal transmission for day i
    if(any(x[x$day==i, "num_I"] >0)){
      # If there are any infected -> run simulation
      step1_output  <- simOneTimeStep(x = x[x$day==i, ],
                                      time_in_A = days_asymptomatic,
                                      n_cities = n_cities,
                                      beta_I = beta_I,
                                      beta_A = beta_A)
      
      # Make sure sensible numbers
      stopifnot(!any(step1_output$num_I > step1_output$tot_N))
      
      # Step 2 - Move Infected
      if(travel==TRUE){
        step2_output <- moveAsymptomatic(step1_output = step1_output,
                                         n_cities = n_cities,
                                         contact_matrix,
                                         frac_travel = frac_travel)
        # Increment by 1 timse-step
        # browser()
        x[x$day == (i + 1), ] <- step2_output
        
      } else {
        # If we don't simulate travel, only do step 1
        x[x$day == (i + 1), ] <- step1_output
        
      }
    }else{
      # browser()
      # If no infected left in system, simply increment by 1 days
      inx <- which(names(x) %in% ("day"))
      x[x$day == (i + 1), -inx] <- x[x$day == i, -inx]
      
    }
    
  }
  
  return(x)
}


simOneTimeStep<- function(x,
                          time_in_I = days_infectious,
                          time_in_A = days_asymptomatic,
                          n_cities,
                          beta_I,
                          beta_A){
  # browser()
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
                          N = N,
                          n_cities = n_cities) # New Infected
  newA <- numberOfNewA(x = x)
  newR_from_I <- numberOfNewR(x$num_I,
                              time_in_I=time_in_I,
                              n_cities = n_cities)  # New Recovered
  # browser()
  newR_from_A <- numberOfNewRFromA(x$num_A,
                                   time_in_A=time_in_A,
                                   n_cities = n_cities)
  
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
                         I, A, S, N, n_cities){
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

numberOfNewR <- function(I, time_in_I, n_cities){
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

numberOfNewRFromA <- function(A, time_in_A, n_cities){
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
  # browser()
  inx <- tmp[, "n_move_out"] > 0
  if(any(inx)) {
    # browser()
    sub_tmp <- tmp[inx, ]
    move_counter <- list()
    
    move_counter <- lapply(1:nrow(sub_tmp), function(k){
      
      # Catch errors in contact matrix
      if(any(is.na(contact_matrix[, sub_tmp$city_id[k]]))) browser()
      
      return(sample(1:n_cities, sub_tmp$n_move_out[k],
                    replace = TRUE,
                    prob = contact_matrix[, sub_tmp$city_id[k]]))
    }) %>%
      unlist()
    
    # Summarize how many times a city receives an Infected person
    new_I_each_city <- data.frame(table(move_counter))
    
    # Add Infected to RECEIVING cities 
    # browser()
    inx <- tmp$city_id %in% new_I_each_city$move_counter
    tmp$num_I[inx] <- new_I_each_city$Freq + tmp$num_I[inx]
  }
  # stopifnot(sum(tmp$n_move_out) == sum(tmp2$Freq, na.rm = TRUE))
  
  return(tmp)
}


calcDist <- function(x1, x2, y1, y2) sqrt((x1 - x2) ^ 2 + (y1 - y2)^2)

gravityCalc <- function(pop1, pop2, d, decay){
  # This function calculates the relative contact between two cities
  # based upon the population of the two cities and the distance separating
  # them.
  # browser()
  (pop1 * pop2) / d^decay
}



makeCityGravityMatrix <- function(n_cities,
                                  pop_vec = cities_population_vec,
                                  dist_mat,
                                  decay = decay){
  # This function creates a matrix showing the normalized contact between all
  # city-pairs.

  
  # Calculate the distance matrix describing each city-pair
  dist_mat <- matrix(data = NA, nrow = n_cities, ncol = n_cities)
  dist_mat <- lapply(1:n_cities, function(i){
    calcDist(x1 = city_locations$x_val[i],
             x2 = city_locations$x_val,
             y1 = city_locations$y_val[i],
             y2 = city_locations$y_val)
  })%>%
    do.call(rbind.data.frame, .) %>%
    `colnames<-`(paste("v", 1:n_cities, sep = ""))
  
  dist_mat <- dist_mat/1000 # Convert to KM
  # Apply gravityCal() to each city pair
  out <- lapply(1:n_cities, function(i){
    gravityCalc(pop_vec[i], pop_vec,
                d = dist_mat[, i], decay = decay)
    
  }) 
  
  # Reshape and rename columns
  contact_matrix <- do.call(rbind.data.frame, out) %>%
    `diag<-` (0) %>%
    `colnames<-`(paste("v", 1:n_cities, sep = ""))
  
  # Normalize so that each column sums to 1
  contact_matrix <- sweep(contact_matrix,
                          2,
                          colSums(contact_matrix),
                          FUN="/")
  # browser()
  return(contact_matrix)
}
# DATA WRANGLING ----------------------------------------------------------

# These functions deal with the data output of the simulations, including
# reshaping and summarizing.

byDayAndItr <- function(x, itr, n_days, n_cities){
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
  
  out$tot_N <- lapply(x, "[[", "tot_N") %>%
    unlist()
  
  out$city_id <- lapply(x, "[[", "city_id") %>%
    unlist()
  # browser()
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



# SUMMARY FUNCTIONS -------------------------------------------------------

summarizeCityI <- function(x){
  # browser()
  out <- split(x, list(x$city_id, x$itr))
  z <- lapply(out, function(x) sum(x$new_I)) %>%
    unlist()%>%
    data.frame(new_internal_I = .)
  z$tot_N = x[!duplicated(x$city_id), "tot_N"]
  z$city_id = 1:n_cities
  z$itr = rep(1:itr, each = n_cities)
  rownames(x) <- (NULL)
  z$max_day_I <- lapply(out, function(x) max(x$I)) %>%
    unlist()
  z$max_day_A <- lapply(out, function(x) max(x$A)) %>%
    unlist()
  z$cumAR <- round(z$new_internal_I / z$tot_N, digits = 2)
  return(z)
}

summarizeOverAllItr <- function(x){
  # browser()
  out <- split(x, f = x$city)
  z <-lapply(out, function(x)mean(x$cumAR*100)) %>%
    unlist() %>%
    data.frame(mean_AR = .,
               city_id = 1:n_cities)
  z$mean_newI <- lapply(out, function(x)mean(x$new_internal_I)) %>%
    unlist()
  z$mean_newA <- lapply(out, function(x)mean(x$max_day_A)) %>%
    unlist()
  return(z)
}

propInfectPerItr <- function(x){
  # browser()
  out <- split(x, f = x$itr)
  lapply(out, function(z){
    sum(z[, "new_internal_I"] >0) / nrow(z)
  }) %>%
    unlist() %>%
    data.frame(prop_infected = .)
  
}

propAllTravel <- function(x){
  # This function finds the mean proportion of cities infected for each "travel
  # fraction" level tested over. Used t.test to find mean instead of quantile
  # since we want the CI around the mean parameter
  lapply(x, function(i){
    # browser()
    out <-  byDayAndItr(i, itr, n_days, n_cities) %>%
      summarizeCityI() %>%
      propInfectPerItr()

    # If all values are the same, t.test throws error, so work around
    flag <- all(out$prop_infected == out$prop_infected[1])
    if(flag){
      x <- mean(out$prop_infected)
      x[2:3] <- x[1]
      attributes(x) <- NULL
      names(x) <- c("2.5%", "97.5%", "mean")
    } else {
      x <- t.test(out$prop_infected)$conf.int
      x[3] <- mean(out$prop_infected)
      attributes(x) <- NULL
      names(x) <- c("2.5%", "97.5%", "mean")
    }
    return(x)
  })
}


extractMeanAndCI <- function(x, z_vec){
  out <- lapply(x, function(x) x["2.5%"]) %>%
    unlist()%>%
    data.frame(low = .)
  out$hi <- lapply(x, function(x) x["97.5%"]) %>%
    unlist()
  out$mean <- lapply(x, function(x) x["mean"]) %>%
    unlist() 
  out$value <- z_vec
  return(out)
}



mapPrepEachDay <- function(x, utm){
  # Prepare a data set for ggplotting for each day. Function needs to be applied
  # to each day
  # browser()
  x$mean_AR <- x$num_I / x$tot_N * 100
  x <- merge(x, utm, by = "city_id")
  return(x)
}


varMonitor <- function(x, itr){
  var_mon <- numeric(itr)
  # browser()
  for(i in 1:itr){
    var_mon[i] <- var(x[1:i])
  }
  return(var_mon)
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



# PLOTTING ----------------------------------------------------------------

plotSensitivity <- function(x, value){
  ggplot() +
    geom_line(data = x,
              aes(x = value, y = mean)) +
    geom_ribbon(data = x,
                alpha = 0.2,
                fill = "red",
                aes(x = value,
                    ymin = low, ymax = hi)) +
    theme_minimal() +
    ylab("Mean proportion of \n cities infected")+
    xlab(paste(value)) +
    ggtitle("Proportion of cities infected") +
    ylim(c(0, 1))
}


mapCities <- function(x){
  ggplot() +
    geom_point(data = x,
               aes(x = x_val, y = y_val,
                   size = tot_N,
                   # col = mean_newA,
                   col = mean_AR)) + 
    scale_colour_distiller(name = "Mean AR rate \nper 100 people",
                           palette = "Spectral",
                           limits=c(0,20)) +
    geom_point(data = x[x$mean_newI == 0, ],
               aes(x = x_val, y = y_val,
                   size = tot_N),
               col = "pink") +
    theme_minimal()
}



mapCitiesByDay <- function(x, data_map){
  # browser()
  ggplot() +
    geom_point(data = data_map[[x]],
               aes(x = x_val, y = y_val,
                   size = tot_N,
                   # col = mean_newA,
                   col = mean_AR)) + 
    scale_colour_distiller(name = "Mean AR rate \nper 100 people",
                           palette = "Spectral",
                           limits=c(0,10)) +
    geom_point(data = data_map[[x]][(data_map[[x]]$num_I == 0 & data_map[[x]]$num_R==0), ],
               aes(x = x_val, y = y_val,
                   size = tot_N),
               col = "pink") +
    theme_minimal() +
    ggtitle(paste0("Day:", x, sep=" "))
  # browser()
  print(paste0("saving plot ", x))
  ggsave(filename = paste0("gif/", x, ".png", sep=""),
         width = 8,height=8,dpi = 96)
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

