##pcpt, class_input

class_input <- function(data, periodlength = NULL, minseglen = 1, distribution = NULL,
                        nsegparam = NULL, Mprior = NULL, Mhyp = NULL, spread = 1,
                        inits = NULL, ...){

  ans               = methods::new("pcpt")
  distribution(ans) = distribution
  data.set(ans)     = data
  ans <- eval(parse(text = paste0("data_value_check.", distribution(ans), "(object = ans)")))

  if(!is.null(periodlength)){
    periodlength(ans) = periodlength
  }else if(is.ts(data)){
    periodlength(ans) = frequency(data)
  }else{
    stop("Period length is not defined either via data as `ts` object or explicitly given as input.")
  }

  if(length(data.set(ans)) < periodlength(ans))
    stop("Length of data is too short or period length is too long.")
  if(periodlength(ans) == 1) stop("Period length must be greater than 1.")

  minseglen(ans)    = minseglen;
  nsegparam(ans)    = nsegparam
  pcpt.prior(ans)   = pcpt.prior.make(Mprior, Mhyp, spread)
  param.prior(ans)  = param.prior.make(distribution, ...)
  ans               = MCMC.options.make(object = ans, ...)
  ans               = Definie.inits(ans, inits, ...)
  results(ans)      = initialise.MCMC.list(n.chains(ans))
  summarised(ans)   = FALSE


  return(ans)
}

## Create/check within cpt prior info for slot
pcpt.prior.make <- function(Mprior = c("pois", "unif"), Mhyp = NULL, spread = 1){

  if(is.null(Mprior)) stop("Mprior cannot be NULL.")
  Mprior <- match.arg(Mprior)

  if(Mprior == "unif"){
    Mhyp <- 0  ##This is not used, set as zero for holding space.
  }else if(Mprior == "pois"){
    if(is.null(Mhyp)){
      Mhyp <- 1
    }else{
      if(length(Mhyp) != 1 | !is.numeric(Mhyp) | any(Mhyp <= 0) | anyNA(Mhyp))
        stop(paste0("Mhyp specified incorrectly for Mprior `",Mprior,"`."))
    }
  }

  if(is.null(spread)) spread = 1
  if(!is.numeric(spread) | length(spread) != 1 | anyNA(spread))
    stop("Hyper-parameter `spread` specified incorrectly.")
  if(spread <= 0) stop("Hyper-parameter `spread` specified incorrectly.")

  return(list(Mprior = Mprior, Mhyp = Mhyp, spread = spread))
}

## Create/check parameter prior info for slot
param.prior.make <- function(distribution = NULL, ...){
  param.prior.make.fn <- paste0("PeriodCPT:::param.prior.make.",distribution,"(...)")
  return(eval(parse(text = param.prior.make.fn)))
}

## Create/check MCMC options for slot
MCMC.options.make <- function(object, n.iter = NULL, n.chains = 1, n.burn = 0,
                              cachesize = 50, quiet = FALSE, ...){
  MCMC.options(object) <- list()

  if(!is.numeric(n.iter)){
    stop("MCMC option - n.iter specified incorrectly.")
  }else if( length(n.iter) != 1 | any(n.iter <= 0) | anyNA(n.iter) |
            any(floor(n.iter) != n.iter) ){
    stop("MCMC option - n.iter specified incorrectly.")
  }
  n.iter(object) <- n.iter

  if(!is.numeric(n.chains)){
    stop("MCMC option - n.chains specified incorrectly.")
  }else if( length(n.chains) != 1 | anyNA(n.chains) | any(n.chains <= 0) |
            any(floor(n.chains) != n.chains) ){
    stop("MCMC option - n.chains specified incorrectly.")
  }
  n.chains(object) <- n.chains

  if(!is.numeric(n.burn)){
    stop("MCMC option - n.burn specified incorrectly.")
  }else if( length(n.burn) != 1 | anyNA(n.burn) | any(n.burn < 0) |
            any(floor(n.burn) != n.burn) ){
    stop("MCMC option - n.burn specified incorrectly.")
  }
  n.burn(object) <- n.burn

  if(!is.numeric(cachesize)){
    stop("MCMC option - cachesize specified incorrectly.")
  }else if(length(cachesize) != 1 | anyNA(cachesize) | any(cachesize <= 0) |
           any(floor(cachesize) != cachesize) ){
    stop("MCMC option - cachesize specified incorrectly.")
  }
  object@MCMC.options[["cachesize"]] <- cachesize

  if(!is.logical(quiet)){
    stop("MCMC option - quiet specified incorrectly.")
  }else if(length(quiet) != 1 | anyNA(quiet)){
    stop("MCMC option - quiet specified incorrectly.")
  }
  quiet(object) <- quiet

  return(object)
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
    return(0)
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
  }else if(m > 1){
    Delta <- N - l*m
    p     <- rgamma(n = m, shape = a, rate = 1)
    delta <- rmultinom(n = 1, size = Delta, prob = p/sum(p))[,1]
    tau0  <- sample.int(n = N, size = 1)
    tau   <- tau0 + cumsum(c(0,delta[-m] + l))
    tau   <- sort(tau %% N) + 1
  }else{
    tau <- 0
  }
  return(tau)
}

## Create/simulate & check chain initial values
Definie.inits <- function(object, inits = NULL, ...){

  inits.pcpt <- list()

  if(is.null(inits)){
    ##Generate inits from prior
    for(i in 1:n.chains(object)) inits.pcpt[[i]] <- rpcpt_single(object)
  }else if(is.function(inits)){
    ##Generate inits according to provided function
    for(i in 1:n.chains(object)){
      if(any(c("pcpt.object", "chain.index") %in% names(list(...))))
        stop("Cannot pass arguments 'pcpt.object' and 'chain.index' from '...' into inits(). These inputs to inits() are managed internally.")
      g_try <- try(inits(pcpt.object = object, chain.index = i, ...), silent = TRUE)
      if(class(g_try) == "try-error"){
        if(grepl("unused argument", g_try)){
          stop("Passing unused arguments to inits function, consider adding '...' to the inputs of your function.")
        }else{
          stop(paste0(strsplit(g_try[1],": ")[[1]][-1], collapse = ": "))
        }
      }else{
        inits.pcpt[[i]] <- g_try
      }
    }
  }else if(any(inits == "ends", na.rm = TRUE)){
    ##Define two chains that are initiated at the ends of the parameter space
    n.chains(object) <- 2
    inits.pcpt[[1]] <- 1
    inits.pcpt[[2]] <- seq(from = 1, by = minseglen(object), len = npcpts.max(object))
  }else if(is.list(inits)){
   if(length(inits) >= n.chains(object)){
      inits.pcpt[1:n.chains(object)] <- inits[1:n.chains(object)]
    }else{
      inits.pcpt <- inits
    }
  }else if(is.numeric(inits)){
    ##Init defined in vector (running only one chain)
    inits.pcpt[[1]] <- inits
  }else{
    stop("Inital values have been specified incorrectly.")
  }

  if(length(inits.pcpt) != n.chains(object))
    stop("Incorrect number of initial values for specified number of chains.")
  if(!all(unlist(lapply(inits.pcpt, is.numeric))))
    stop("Class of at least one inits is not numeric.")
  if(anyNA(unlist(inits.pcpt)))
    stop("Inits must not contain missing NA values.")
  m <- unlist(lapply(inits.pcpt, length))
  if(any(m<1 | m>npcpts.max(object)))
    stop("Incorrect number of within period changepoints specified by inits.")
  for(i in 1:n.chains(object)){
    if(any(floor(inits.pcpt[[i]]) != inits.pcpt[[i]]))
      stop("In inits, within period cpts must be whole numbers.")
    if((min(inits.pcpt[[i]]) < 1) | (max(inits.pcpt[[i]]) > periodlength(object)))
      stop("In inits, within period cpts must be between 1 and period length.")
    if(any(diff(c(inits.pcpt[[i]], inits.pcpt[[i]][1]+periodlength(object))) < minseglen(object)))
      stop("In inits, within period cpts are not ordered or do not satisfy minimum segment length condition.")
  }
  names(inits.pcpt) <- as.character(1:n.chains(object))
  MCMC.inits(object) <- inits.pcpt
  return(object)

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


