PeriodCPT_extend <- function(object, newiters = 1e4){
  ##--Input check--
  if(class(object) != "pcpt") stop("Unexpected class of `object`.")
  if(length(object@date) == 0) stop("Argument `object` does not appear to be an output from a PeriodCPT funciton.")
  if(length(newiters) != 1 | !is.numeric(newiters)) stop("Argument `newiters` must be a single positive integer.")
  if(newiters <= 0 | floor(newiters) != newiters) stop("Argument `newiters` must be a single positive integer.")

  ##Evaluate inital values as last case from last scenario.
  first.inits <- MCMC.inits(object)
  MCMC.inits(object) <- MCMC.last(object)
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
