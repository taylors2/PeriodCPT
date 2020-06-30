PeriodCPT_extend <- function(object, newiters = 1e4){
  ##--Input check--
  if(class(object) != "pcpt") stop("Unexpected class of `object`.")
  if(length(newiters) != 1 | !is.numeric(newiters)) stop("Argument `newiters` must be a single positive integer.")
  if(newiters <= 0 | floor(newiters) != newiters) stop("Argument `newiters` must be a single positive integer.")

  ##Evaluate inital values as last case from last scenario.
  first.inits <- MCMC.inits(object)
  if(summarised(object) == TRUE){
    MCMC.inits(object) <- MCMC.last(object)
    if(length(MCMC.inits(object)) == 0) stop("Cannot determine last sample to initiate next batch.")
  }else{
    MCMC.inits(object) <- Calc_MCMC_last(object)
  }
  newburn <- n.burn(object) + n.iter(object) ##add last iters to burn (hold in temp variable)

  ##Reset slots
  n.burn(object) <- 0               ##Set burn to zero (ie, treat up to current iteration as burn-in)
  n.iter(object) <- newiters        ##Set n.iter to input value
  MCMC.last(object) <- list()
  pcpt.mode(object) <- numeric()
  param.mode(object) <- array(NA, dim = c(0, 0))
  fit(object) <- numeric()
  results(object) <- initialise.MCMC.list(n.chains(object))
  summarised(object) <- FALSE

  ##Run new batch of iterations
  object <- PeriodCPT.main(object)

  ##Assign back relavent information
  n.burn(object) <- newburn
  MCMC.inits(object) <- first.inits
  return(object)
}
