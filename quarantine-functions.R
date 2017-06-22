### Written by Emil Tosti 
## emtosti@gmail.com
## version .01
# 20.06.17

library(reshape2)
library(ggplot2)
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

determineInfectlevel <- function(x, inx.quara, h.level = 50, l.level = 10){
  #browser()
  m.level = round(mean(c(h.level,l.level)))
  #number of infected in a given city
  n.infec  <- sum(x$num_I[x$city_id==inx.quara],x$new_I[x$city_id==inx.quara], na.rm = TRUE)
  if(n.infec >= h.level){
    infect.level = 2
  }
  else if(n.infec >= m.level){
    infect.level = 1
  }
  else{
    infect.level = 0
  }
  return(infect.level)
}


putCityUnderQuarantine <- function(x = city, start.threshold = 50, end.threshold = 20){
  #Determine which cities are under quarantine
  #browser()
  #level of infection in all, non quarantined, cites
  inx.quara <- which(x$in.quara == 1)
  if(length(inx.quara)>0){
    for(i in 1:length(inx.quara)){
      ##calling the freq of infection within the city
      #infec.freq <- determineInfectFreq(x, inx.quara[i])
      #evaluation of the quarantine should be lifted
      infect.in.city <- sum(x$num_I[x$city_id == inx.quara[i]],x$new_I[x$city_id == inx.quara[i]])
      if(infect.in.city < end.threshold){
       # browser()
        x$in.quara[x$city_id == inx.quara[i]] <- FALSE
        x$sec.level[x$city_id == inx.quara[i]] = determineInfectlevel(x, inx.quara[i], 
                                                                       h.level = start.threshold, l.level = end.threshold)
      }
    }
  }
  #level of infection in all, non quarantined, cites
  inx.non.quara <- which(!x$in.quara)
  if(length(inx.non.quara)>0){
  #look at which cities have a high enough level of infection to be put under quarantine
    for(i in 1:length(inx.non.quara)){
      #browser()
      # <- determineInfectFreq(x, inx.non.quara[i])
      infect.in.city <- sum(x$num_I[x$city_id == inx.non.quara[i]],x$new_I[x$city_id == inx.non.quara[i]])
      if(infect.in.city >= start.threshold){
        x$in.quara[x$city_id == inx.non.quara[i]] <- TRUE
        #browser()
        x$sec.level[x$city_id == inx.non.quara[i]] = determineInfectlevel(x, inx.non.quara[i], 
                                                                           h.level = start.threshold, l.level = end.threshold)
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
    #browser()
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
                         frac_travel,
                         start.quarantine, 
                         end.quarantine){
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
    step1_output <- putCityUnderQuarantine(step1_output,
                                           start.threshold = start.quarantine, 
                                           end.threshold  = end.quarantine
                                           )
    
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
  out$in.quara <- lapply(x, "[[", "in.quara") %>%
    unlist()
  out$sec.level <- lapply(x, "[[", "sec.level") %>%
    unlist()
  out$itr <- rep(1:itr, each = n_days*n_cities)
  out$day <- rep(1:n_days)
  
  return(out)
}


cumInfected <- function(x){
  # Get the cummulative number of infected per city per iteration.
  lapply(x, function(x){
    #browser()
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
                      in.quara = rep(FALSE, n_time),
                      sec.level = rep(0, n_time)
                      )
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
  inx.q <- which(cities$in.quara == 1)
  if(length(inx.q)>0 ){
    #browser()
    for( i in 1:length(inx.q)){
      if(cities$sec.level[cities$city_id==inx.q[i]] == 2){
      city_contact_matrix[inx.q[i], ] = city_contact_matrix[inx.q[i], ] * 0.1
      }
      else if(cities$sec.level[cities$city_id==inx.q] == 1){
        city_contact_matrix[inx.q[i], ] = city_contact_matrix[inx.q[i], ] * 0.35
      }
      else{
        city_contact_matrix[inx.q[i], ] = city_contact_matrix[inx.q[i], ] * 0.45
      }
    }
  }
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


plotMeanInfection <- function(out){
  plots <- list()
  for(i in 1:max(out$city_id)){  
    plots[[i]] <- ggplot(data = out[out$city_id==i,], aes(x = day, y = I)) +
    stat_summary(aes(y = I, group=1), fun.y=mean, colour=i, geom="line", group=1) +
    theme_minimal()
  }
  return(plots)
}

# in.quara.city1 <- data.frame(
#           itr1 = rep(0, n_days),
#           itr2 = rep(0, n_days),
#           itr3 = rep(0, n_days),
#           itr4 = rep(0, n_days),
#           itr5 = rep(0, n_days),
#           itr6 = rep(0, n_days),
#           itr7 = rep(0, n_days),
#           itr8 = rep(0, n_days),
#           itr9 = rep(0, n_days),
#           itr10 = rep(0, n_days),
#           itr11 = rep(0, n_days),
#           itr12 = rep(0, n_days)
#            )
# 
# 
# 
# mdf <- melt(out[out$city_id==1,], id.vars=c("day", "itr"), measure.var = c("in.quara"))
# mdf[mdf$day ==1,]
# 
# 
# cdf <- cast(mdf, day )
# head(cdf)
# 
# 
# cdf$in.quara <- as.integer(cdf$in.quara)
# 
# ggplot(data = cdf$in.quara[])
# 
# city.1 <- data.frame(day = 1:length(unique(cdf$day)),
#                      in.quara = rep(0, length(unique(cdf$day))))
# for(i in 1:length(unique(cdf$day))){
#   temp = vector(length =12, mode = "integer")
#   temp2 <- cdf[cdf$day==i,]
#   print(temp2)
#   for(n in 1:length(unique(cdf$itr))){
#     length(cdf$in.quara[cdf$itr==n])
#     # temp[n] <- 
#   }
#     #city.1$in.quara[i] <- Mode(temp)
# }

# reshapeOutput <- function(out){
#   mdf <- melt(out, id.vars=c("day","city_id","itr"), measure.vars = c("S", "I", "new_I"))
# }
# 
# 
# 
# cdf <- cast(mdf, formular = day + day + itr ~ variable)
# ggplot(data = cdf, aes(x = day, y = I)) + 
#   stat_summary(aes(y = I, col=cdf$city_id), fun.y=mean, geom="line") 
#   
###------- From reposetory ------ ######
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

