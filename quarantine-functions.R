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

determineInfectFreq <- function(city, inx.quara){
  #number of infected in a given city
  temp.store.infec  <- city$num_I[city$id==inx.quara] + city$num_I[city$id==inx.quara]
  #all people in a city
  temp.store.n.city <- temp.store.infec + city$tot_N
  #the freq of infection
  infec.freq <- temp.store.infec/temp.store.n.city
  return(infec.freq)
}

bleedingOfCity <- function(city, quarantine, inx.quara = inx.quara, security.level = security.level){
  n.people <- sum(city$num_S[city$city_id==inx.quara],
                  city$num_I[city$city_id==inx.quara],
                  city$new_I[city$city_id==inx.quara],
                  city$num_R[city$city_id==inx.quara]
                  )
  ### use an appropriate selection
  
}


putCityUnderQuarantine <- function(city = city, quarantine = quarantine, start.threshold = 0.55, end.threshold = 0.01){
  #Determine which cities are under quarantine
  inx.quara <- which(!quarantine$quara)
  for(i in 1:length(inx.quara)){
    #calling the freq of infection within the city
    infec.freq <-determineInfectFreq(city, inx.quara[i])
    #evaluation of the quarantine should be lifted
    if(infec.freq < end.threshold){
      quarantine$quara[quarantine$quara == inx.quara[i]] <- FALSE
    }
  }
  #evaluate freq of infection in all, non quarantined, cites
  inx.non.quara <- which(quarantine$quara == 0)
  freq.of.infect <- apply(X = inx.non.quara, 2, function(x){determineInfectFreq(city,x)})
  #look at which cities have a high enough level of infection to be put under quarantine
  new.quara.city <- which(freq.of.infect >= start.threshold)
  
  quarantine$quara[quarantine$city.id == new.quara.city] <- TRUE
  
  ###--- Determine bleeding --###
}









