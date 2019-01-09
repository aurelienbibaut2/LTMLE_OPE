# Two helper functions
# Expit and logit
expit <- function(x) { 1 / (1 + exp(-x)) } 
logit_scalar <- Vectorize(function(x) { 
  if(x > 0 & x < 1){
    return(log(x / (1 - x)))
  }else if(x <= 0){
    return(-Inf)
  }else{
    return(Inf)
  }
})
logit <- function(x){
  if(is.null(dim(x))){
    return(logit_scalar(x))
  }else{
    return(array(logit_scalar(x), dim=dim(x)))
  }
}

# Soften weights
soften <- function(a, alpha){
  a^alpha / sum(a^alpha) 
}

# One helper function
distance_to_interval <- function(interval, x){
  lb <- interval[1]; ub <- interval[2]
  if(x >= lb && x <= ub){
    return(0)
  }else{
    return(min(abs(x-lb), abs(x-ub)))
  }
}