##pcpt, class_input

class_input <- function(data, periodlength, minseglen, distribution,
                        Mprior, Mhyp, spread, ...){

  ans = new("pcpt")
  data.set(ans)     = data
  periodlength(ans) = periodlength
  minseglen(ans)    = minseglen;
  distribution(ans) = distribution

  pcpt.prior(ans)   = pcpt.prior.make(Mprior, Mhyp, spread)
  param.prior(ans)  = param.prior.make(distribution, ...)
  MCMC.options(ans) = MCMC.options.make(...)

  for(i in 1:n.chains(ans)){

  }

}


pcpt.prior.make <- function(Mprior = c("unif","pois"), Mhyp, spread = 1){

  Mprior <- match.arg(Mprior)

  if(Mprior == "unif"){
    if(!missing(Mhyp))
      warning("Specified argument `Mhyp` is ignored for Mprior `unif`.")
    Mhyp <- 0
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

#######

param.prior.make <- function(distribution, ...){
  param.prior.make.fn <- get(paste0("PeriodCPT:::param.prior.",distribution))
  return(param.prior.make.fn(...))
}

##########

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







