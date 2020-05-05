##pcpt, class_input

class_input <- function(data, periodlength, minseglen, distribution,
                        Mprior, Mhyp, spread, inits, ...){

  ans               = new("pcpt")
  data.set(ans)     = data
  periodlength(ans) = periodlength
  minseglen(ans)    = minseglen;
  distribution(ans) = distribution

  pcpt.prior(ans)   = pcpt.prior.make(Mprior, Mhyp, spread)
  param.prior(ans)  = param.prior.make(distribution, ...)
  MCMC.options(ans) = MCMC.options.make(...)
  MCMC.inits(ans)   = Definie.inits(ans, inits, ...)
  MCMC.chains(ans)  = initialise.MCMC.list(n.chains(ans))


  return(ans)
}

## Create/check within cpt prior info for slot
pcpt.prior.make <- function(Mprior = c("unif","pois"), Mhyp, spread = 1){

  Mprior <- match.arg(Mprior)

  if(Mprior == "unif"){
    Mhyp <- 0  ##This is not used, set as zero for holding space.
  }else if(Mprior == "pois"){
    if(missing(Mhyp)){
      Mhyp <- 1
    }else{
      if(length(Mhyp) != 1 || !is.numeric(Mhyp) || any(Mhyp <= 0))
        stop(paste0("Mhyp specified incorrectly for Mprior `",Mprior,"`"))
    }
  }else{
    stop("Implementation Error: Mprior `",Mprior,"` is not supported.")
  }

  if(!is.numeric(spread) || length(spread) != 1)
    stop("pcpt: spread specified incorrectly.")
  if(spread <= 0) stop("pcpt: spread specified incorrectly.")

  return(list(Mprior = Mprior, Mhyp = Mhyp, spread = spread))
}

## Create/check parameter prior info for slot
param.prior.make <- function(distribution, ...){
  param.prior.make.fn <- get(paste0("PeriodCPT:::param.prior.",distribution))
  return(param.prior.make.fn(...))
}

## Create/check MCMC options for slot
MCMC.options.make <- function(n.iter, n.chain, n.burn, cachesize, quiet){

  if(missing(n.iter)){
    stop("MCMC option - n.iter not specified.")
  }else{
    if(!is.numeric(n.iter) || length(n.iter) != 1 ||
       any(n.iter <= 0) || any(floor(n.iter) != n.iter) )
      stop("MCMC option - n.iter specified incorrectly")
  }

  if(missing(n.chains)){
    n.chains <- 1
  }else{
    if(!is.numeric(n.chains) || length(n.chains) != 1 ||
       any(n.chains <= 0) || any(floor(n.chains) != n.chains) )
      stop("MCMC option - n.chains specified incorrectly")
  }

  if(missing(n.burn)){
    n.burn <- 0
  }else{
    if(!is.numeric(n.burn) || length(n.burn) != 1 ||
       any(n.burn < 0) || any(floor(n.burn) != n.burn) )
      stop("MCMC option - n.burn specified incorrectly")
  }

  if(missing(cachesize)){
    cachesize <- 50
  }else{
    if(!is.numeric(cachesize) || length(cachesize) != 1 ||
       any(cachesize <= 0) || any(floor(cachesize) != cachesize) )
      stop("MCMC option - cachesize specified incorrectly")
  }

  if(missing(quiet)){
    quiet <- FALSE
  }else{
    if(!is.numeric(quiet) || length(quiet) != 1 || !is.logical(quiet) || anyNA(quiet))
      stop("MCMC option - quiet specified incorrectly")
  }

  return(
    list(n.iter = n.iter, n.chain = n.chain,
         n.burn = n.burn, cachesize = cachesize,
         quiet = quiet)
  )
}

## Simulate the number of within period changepoints given pcpt.prior info
rMprior <- function(n, object, ...){
  if(pcpt.prior(object)$Mprior){
    return(sample.int(npcpts.max(object), size = n, replace = TRUE))
  }else if(pcpt.prior(object)$Mprior == "pois"){
    m <- rep(0,n)
    for(i in 1:n){
      while(m[i]>npcpts.max(object)){
        m[i] <- 1 + rpois(n = 1, rate = pcpt.prior(object)$Mhyp[1])
      }
    }
    return(m)
  }else{
    stop("Unrecognised prior for the number of within period changepoints.")
  }
}

## Simulate the within period cpt positions given pcpt.prior info
rpcpt_single <- function(object){
  N <- periodlength(object)
  l <- minseglen(object)
  a <- pcpt.prior(object)$spread
  mMax  <- npcpt.max(object)
  m <- rMprior(1, object)
  if(m == 1){
    tau <- 0
  }else{
    Delta <- N - l*m
    p     <- rgamma(n = m, shape = 1, rate = 1)
    delta <- rmultinom(n = 1, size = Delta, prob = p/sum(p))[,1]
    tau0  <- sample.int(n = N, size = 1) - 1
    tau   <- tau0 + cumsum(c(0,delta[-m] + l))
    tau   <- sort(tau %% N)
  }
  return(tau)
}

## Create/simulate & check chain initial values
Definie.inits <- function(object, inits, ...){

  inits.pcpt <- list()

  if(missing(inits) | is.null(inits)){
    ##Generate inits from prior
    for(i in 1:n.chains(object)) inits.pcpt[[i]] <- rpcpt_single(object)
  }else if(any(inits == "ends")){
    ##Define two chains that are initiated at the ends of the parameter space
    n.chains(object) <- 2
    inits.pcpt[[1]] <- 0
    inits.pcpt[[2]] <- seq(from = 0, by = minseglen(object), len = npcpt.max(object))
  }else if(is.function(inits)){
    ##Generate inits according to provided function
    for(i in 1:n.chains(object)){
      if(any(c("pcpt.object", "chain.index") %in% names(list(...))))
        stop("Cannot pass arguments `pcpt.object` and `chain.index` from `...` into inits().
             These inputs to inits() are managed internally.")
      inits.pcpt[[i]] <- inits(pcpt.object = object, chain.index = i, ...)
    }
  }else if(is.list(inits)){
    ##Inits defined within a list
    if(length(inits) < n.chains(object)){
      stop("Too few inital values provided for specified number of chains")
    }
    inits.pcpt[1:n.chains(object)] <- inits[1:n.chains(object)]
  }else if(is.numeric(inits)){
    ##Init defined in vector (running only one chain)
    if(n.chains(object) != 1)
      stop("Too few inital values provided for specified number of chains.")
    inits.pcpt[[1]] <- inits
  }

  if(length(inits.pcpt) != n.chains(object))
    stop("Incorrect number of initial values for specified number of chains.")
  if(!all(unlist(lapply(tmp, class)) %in% c("numeric","integer")))
    stop("Class of at least one inits is not numeric or interger.")
  m <- unlist(lapply(inits.pcpt, length))
  if(any(m<1 || m>npcpt.max(object)))
    stop("Incorrect number of within period changepoints specified by inits.")
  for(i in 1:n.chains(object)){
    if(any(floor(inits.pcpt[[i]]) != inits.pcpt[[i]]))
      stop("In inits, within period cpts must be whole numbers.")
    if(inits.pcpt[[i]][1] < 0 || inits.pcpt[[i]][1] >= npcpt.max(object))
      stop("In inits, within period cpts must be within [0. period lenght).")
    if(any(diff(c(inits.pcpt[[i]], inits.pcpt[[i]][1]+periodlength(object))) < l))
      stop("In inits, within period cpts does not satisfy minimum segment length condition.")
  }

  return(inits.pcpt)

}

## Create blank list for MCMC output
initialise.MCMC.list <- function(n){
  chains <- vector("list", n)
  names(chains) <- 1:n
  return(chains)
}


populate.chain <- function(object, C.chain.output, blank = -1){
  TAU <- matrix(C.chain.output, ncol = npcpt.max(object), byrow = TRUE)
  TAU[TAU == blank] <- NA
  colnames(TAU) <- paste0("tau",1:npcpt.max(object))
  TAU <- TAU[ , !apply(is.na(TAU), 2, all), drop = FALSE]
  return(TAU)
}


