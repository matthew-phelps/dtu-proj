S <- 0
E <- 1
I <- 2

itr <- 2
beta <- 0.9
E_prob <- 1
I_time <- 5
n_animal <- 100
n_days <- 20000
prob_inf <- function(beta, I, N){
  1 - exp(-beta * I/N)
}


makeHerd <- function(n_animal= n_animal){
  data.frame(id = 1:n_animal,
                   status = S,
                   counter = 0)
}

herd <- makeHerd(n_animal)

simOverTime<- function(x = herd, days = n_days, N = n_animal, counter_time = I_time,
                       seed_I = 1, E_prob = 0.6){
  # browser()
 
  # Initialize 1st infection in herd
  x[1:seed_I, "status"] <- I
  x[1:seed_I, "counter"] <- counter_time
  
  # Set up tracking for each compartment
  S_t <- numeric(days)
  S_t[1] <- sum(x[, "status"]==S)
  
  E_t <- numeric(days)
  E_t[1] <- sum(x[, "status"]==E)

  I_t <- numeric(days)
  I_t[1] <- sum(x[, "status"]==I)
  
  new_I <- numeric(days)
  new_I[1] <- sum(x[, "status"]==I)
  
  # Loop over days
  for(i in 1:(days)){
    # browser()
    
    # Number of current infecitons
    inx_I <- x[, "status"]==I
    inx_S <- x[, "status"]==S
    inx_E <- x[, "status"]==E
    
    # New E
    prob <- pi(beta, sum(I), N)
    newE <- inx_S & rbinom(n_animal, 1, prob = prob)==1
    
    # Move I
    inx_moveI <- x[inx_I, "counter"]==1
    
    # Move E
    inx_moveE <- inx_E & rbinom(n_animal, 1, prob = E_prob)
    
    
    # Update counters
    x[inx_I, "counter"] <-  x[inx_I, "counter"] - 1
    
    # Move animals
    x[inx_moveI, c("status")] <- c(S)
    x[newE, c("status")] <- c(E)
    x[inx_moveE, c("status")] <- (I)
    x[inx_moveE, c("counter")] <- counter_time
    
    # Store counts of each status
    S_t[i] <- sum(x[, "status"]==S)
    E_t[i] <- sum(x[, "status"]==E)
    I_t[i] <- sum(x[, "status"]==I)
    new_I[i] <- sum(newE)
  }
  
  return(list(S_t = S_t,
              E_t = E_t,
              I_t = I_t,
              new_I = new_I))
}


out <- simOverTime()
plot(out$S_t, type = "l", ylim = c(0, n_animal))
lines(out$I_t, col  ="darkred")
lines(out$E_t, col = "darkgreen")


simOverItr <- function(itr){
  browser()
  out <- list()
  for (j in 1:itr){
    out[[j]] <- simOverTime()
  }
  return(out)
}
simOverItr(itr)



