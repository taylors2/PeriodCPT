##pcpt, class_input

class_input <- function(data, periodlength, minseglen, distribution, nsegparam,
                        Mprior, Mhyp, spread, inits, ...){

  ans               = methods::new("pcpt")
  distribution(ans) = distribution
  data.set(ans)     = data
  ans <- eval(parse(text = paste0("data_value_check.", distribution(ans), "(object = ans)")))

  if(missing(periodlength)) periodlength <- NULL
  if(!is.null(periodlength)){
    periodlength(ans) = periodlength
  }else if(is.ts(data)){
    periodlength(ans) = frequency(data)
  }else{
    stop("Period length is not defined either via data as `ts` object or explicitly given as input.")
  }

  if(length(data.set(ans)) < periodlength(ans))
    stop("Length of data is too short or period length is too long.")

  if(missing(minseglen)){
    minseglen(ans)    = 1
  }else if(is.null(minseglen)){
    minseglen(ans)    = 1
  }else{
    minseglen(ans)    = minseglen;
  }

  nsegparam(ans)    = nsegparam
  pcpt.prior(ans)   = pcpt.prior.make(Mprior, Mhyp, spread)
  param.prior(ans)  = param.prior.make(distribution, ...)
  MCMC.options(ans) = MCMC.options.make(...)
  MCMC.inits(ans)   = Definie.inits(ans, inits, ...)
  results(ans)      = initialise.MCMC.list(n.chains(ans))
  summarised(ans)   = FALSE


  return(ans)
}

## Create/check within cpt prior info for slot
pcpt.prior.make <- function(Mprior = c("pois", "unif"), Mhyp, spread = 1){

  if(is.null(Mprior)) stop("Mprior cannot be NULL.")
  Mprior <- match.arg(Mprior)

  if(Mprior == "unif"){
    Mhyp <- 0  ##This is not used, set as zero for holding space.
  }else if(Mprior == "pois"){
    if(missing(Mhyp) | is.null(Mhyp)){
      Mhyp <- 1
    }else{
      if(length(Mhyp) != 1 | !is.numeric(Mhyp) | any(Mhyp <= 0) | anyNA(Mhyp))
        stop(paste0("Mhyp specified incorrectly for Mprior `",Mprior,"`."))
    }
  }else{
    stop(paste0("Implementation Error: Mprior `",Mprior,"` is not supported."))
  }

  if(missing(spread) | is.null(spread)) spread = 1
  if(!is.numeric(spread) | length(spread) != 1 | anyNA(spread))
    stop("Hyper-parameter `spread` specified incorrectly.")
  if(spread <= 0) stop("Hyper-parameter `spread` specified incorrectly.")

  return(list(Mprior = Mprior, Mhyp = Mhyp, spread = spread))
}

## Create/check parameter prior info for slot
param.prior.make <- function(distribution, ...){

  param.prior.make.fn <- paste0("PeriodCPT:::param.prior.make.",distribution,"(...)")
  return(eval(parse(text = param.prior.make.fn)))
}

## Create/check MCMC options for slot
MCMC.options.make <- function(n.iter, n.chains, n.burn,
                              cachesize, quiet, ...){

  if(missing(n.iter) | is.null(n.iter) | !is.numeric(n.iter)){
    stop("MCMC option - n.iter not specified.")
  }else if( length(n.iter) != 1 | any(n.iter <= 0) | anyNA(n.iter) |
            any(floor(n.iter) != n.iter) ){
    stop("MCMC option - n.iter specified incorrectly.")
  }

  if(missing(n.chains) | is.null(n.chains)){
    n.chains <- 1
  }else if(!is.numeric(n.chains)){
    stop("MCMC option - n.chains specified incorrectly.")
  }else if( length(n.chains) != 1 | anyNA(n.chains) | any(n.chains <= 0) |
            any(floor(n.chains) != n.chains) ){
    stop("MCMC option - n.chains specified incorrectly.")
  }

  if(missing(n.burn) | is.null(n.burn)){
    n.burn <- 0
  }else if(!is.numeric(n.burn)){
    stop("MCMC option - n.burn specified incorrectly.")
  }else if( length(n.burn) != 1 | anyNA(n.burn) | any(n.burn < 0) |
            any(floor(n.burn) != n.burn) ){
    stop("MCMC option - n.burn specified incorrectly.")
  }

  if(missing(cachesize) | is.null(cachesize)){
    cachesize <- 50
  }else if(!is.numeric(cachesize)){
    stop("MCMC option - cachesize specified incorrectly.")
  }else if(length(cachesize) != 1 | anyNA(cachesize) | any(cachesize <= 0) |
           any(floor(cachesize) != cachesize) ){
    stop("MCMC option - cachesize specified incorrectly.")
  }

  if(missing(quiet) | is.null(quiet)){
    quiet <- FALSE
  }else if(!is.logical(quiet)){
    stop("MCMC option - quiet specified incorrectly.")
  }else if(length(quiet) != 1 | anyNA(quiet)){
    stop("MCMC option - quiet specified incorrectly.")
  }

  return(
    list(n.iter = n.iter, n.chains = n.chains,
         n.burn = n.burn, cachesize = cachesize,
         quiet = quiet)
  )
}

## Simulate the number of within period changepoints given pcpt.prior info
rMprior <- function(n, object){
  if(pcpt.prior(object)$Mprior == "unif"){
    return(sample.int(npcpts.max(object), size = n, replace = TRUE))
  }else if(pcpt.prior(object)$Mprior == "pois"){
    lam <- pcpt.prior(object)$Mhyp[1]
    m <- qpois(p = runif(n = n,
                         min = ppois(0,lam),
                         max = ppois(npcpts.max(object),lam) ),
               lambda = lam)
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
  m <- rMprior(1, object)
  if(m == 1){
    tau <- 1
  }else{
    Delta <- N - l*m
    p     <- rgamma(n = m, shape = a, rate = 1)
    delta <- rmultinom(n = 1, size = Delta, prob = p/sum(p))[,1]
    tau0  <- sample.int(n = N, size = 1)
    tau   <- tau0 + cumsum(c(0,delta[-m] + l))
    tau   <- sort(tau %% N) + 1
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
    inits.pcpt[[1]] <- 1
    inits.pcpt[[2]] <- seq(from = 1, by = minseglen(object), len = npcpts.max(object))
  }else if(is.function(inits)){
    ##Generate inits according to provided function
    for(i in 1:n.chains(object)){
      if(any(c("pcpt.object", "chain.index") %in% names(list(...))))
        stop(paste0("Cannot pass arguments `pcpt.object` and `chain.index` from `...` into inits().",
             " These inputs to inits() are managed internally."))
      inits.pcpt[[i]] <- inits(pcpt.object = object, chain.index = i, ...)
    }
  }else if(is.list(inits)){
    ##Inits defined within a list
    if(length(inits) < n.chains(object)){
      stop("Too few inital values provided for specified number of chains.")
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
  if(!all(unlist(lapply(inits.pcpt, is.numeric))))
    stop("Class of at least one inits is not numeric.")
  m <- unlist(lapply(inits.pcpt, length))
  if(any(m<1 | m>npcpts.max(object)))
    stop("Incorrect number of within period changepoints specified by inits.")
  for(i in 1:n.chains(object)){
    if(any(floor(inits.pcpt[[i]]) != inits.pcpt[[i]]))
      stop("In inits, within period cpts must be whole numbers.")
    if((min(inits.pcpt[[i]]) < 1) | (max(inits.pcpt[[i]]) > periodlength(object)))
      stop("In inits, within period cpts must be within [1, period length].")
    if(any(diff(c(inits.pcpt[[i]], inits.pcpt[[i]][1]+periodlength(object))) < minseglen(object)))
      stop("In inits, within period cpts does not satisfy minimum segment length condition.")
  }
  names(inits.pcpt) <- as.character(1:n.chains(object))
  return(inits.pcpt)

}

## Create blank list for MCMC output
initialise.MCMC.list <- function(n){
  chains <- vector("list", n)
  names(chains) <- 1:n
  return(chains)
}


populate.chain <- function(object, C.chain.output, blank = -1, ...){
  TAU <- matrix(C.chain.output, ncol = npcpts.max(object), byrow = TRUE)
  TAU[TAU == blank] <- NA
  colnames(TAU) <- paste0("tau",1:npcpts.max(object))
  TAU <- TAU[ , !apply(is.na(TAU), 2, all), drop = FALSE]
  return(TAU)
}


