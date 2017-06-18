
# 2 -----------------------------------------------------------------------
rm(list=ls())
n_herds <- 50
S <- FALSE
I <- TRUE
makeHerd <- function(){
  data.frame(id = 1:n_herds,
             status = S)
}

moveAnimal <- function(){
  # browser()
  out <- data.frame(id_send = sample(1:n_herds, 500, replace = TRUE),
                    id_rec = sample(1:n_herds, 500, replace = TRUE),
                    n_animals = sample(1:5, 500, replace = TRUE))
  inx <- out$id_send == out$id_rec
  out$id_rec[inx] <- out$id_rec[inx] + 1
  out
}

m <- moveAnimal()
herd <- makeHerd()
prob <- function(prev, N){
  1 - (1 - prev)^N
}
herd[c(4, 10), "status"] <- I

simInfect <- function(x = herd, move = m){
  
  inx_inf <- x[, "status"] == I
  id_infect <- x[inx_inf, "id"]
  new_inf <- numeric(sum(inx_inf))
  for(i in 1:length(id_infect)){
    id <- id_infect[i]
    inx_move <- move[, "id_send"]==id
    
    # Prev in infecting herd
    prev <- .10
    
    rec <- move[inx_move, ]
    # Total animals by receiving ID
    rec$tot_animal <- ave(rec$n_animals, rec$id_rec, FUN = function(x) sum(x))
    
    # Remove duplicate rec_ids
    rec <- rec[!(duplicated(rec$id_rec)), ]
    rec$prob <- prob(prev = .1, N = rec$tot_animal)
    inx_infc <- rbinom(nrow(rec), 1, prob = rec$prob)==1
    
    inx_inf_herd <- (x$id %in% rec$id_rec[inx_infc])
    # browser()
    x[inx_inf_herd, "status"] <- I
    new_inf[i] <- sum(inx_infc)
  }
  return(sum(new_inf))
}

out <- numeric(100)
for(j in 1:100){
  out[j] <- (simInfect())
}