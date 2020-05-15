##Wrapper function for PeriodCPT
PeriodCPT <- function(data,
                      distribution = c("bern","pois","norm","norm","mean","var"),
                      ...){
  distribution = match.arg(distribution)
  if(missing(data)) stop("data is missing.")
  ans <- eval(parse(text=paste0("PeriodCPT.",distribution,"(data = data, ...)")))
  ans <- eval(parse(text = paste0("SummariseOutput.",distribution,"(ans)")))
  return(ans)
}

##The main Pernd distibution specific functioiodCPT function that calls the C code
PeriodCPT.main <- function(object){

  ##Input check
  if(class(object) != "pcpt") stop("Unexpected class of `object`.")
  #All other formal checks should be done by class specific operators
  # ans


  ##Format inputs to C
  BLANK <- -1
  offset <- start(data.set(object))[2]
  time <- seq(from = offset, by = 1, length.out = length(data.set(object)))
  time2C <- ((time - 1) %% periodlength(object)) + 1
  inits2C <- rep(BLANK, n.chains(object)*npcpts.max(object))
  for(i in 1:n.chains(object)){
    tau <- MCMC.inits(object)[[as.character(i)]]
    inits2C[(n.chains(object))* + (1:length(tau))] <- tau
  }

  ##Run C code
  draw <- .C("PeriodCPT_RJMCMC",
             data    = as.numeric(data.set(object)),
             time    = as.integer(time2C),
             n       = as.integer(length(data.set(object))),
             N       = as.integer(periodlength(object)),
             l       = as.integer(minseglen(object)),
             Mmax    = as.integer(npcpts.max(object)),
             Mdist   = as.character(pcpt.prior(object)$Mprior),
             Mhyp    = as.numeric(pcpt.prior(object)$Mhyp),
             spread  = as.numeric(pcpt.prior(object)$spread),
             Pdist   = as.character(distribution(object)),
             Phyp    = as.numeric(param.prior(object)),
             inits   = as.integer(inits2C),
             nchains = as.integer(n.chains(object)),
             nburn   = as.integer(n.burn(object)),
             niter   = as.integer(n.iter(object)),
             cache   = as.integer(MCMC.options(object)$cachesize),
             quiet   = as.integer(MCMC.options(object)$quiet),
             blank   = as.integer(BLANK),
             error   = as.integer(0),
             draw    = vector("integer", length = n.chains(object) *
                                n.iter(object) * npcpts.max(object))
  )

  ##Error checking
  if(draw$error != 0)
    warning(paste0("PeriodCPT C code had an error. Error code: ",draw$error,"."))

  ##Format C output and put into object
  len.chain.samples <- npcpts.max(object)*n.iter(object)
  for(i in 1:n.chains(object)){
    MCMC.chain(object, i) <- populate.chain(object,
                                            draw$draw[(i-1)*len.chain.samples + (1:len.chain.samples)],
                                            blank = BLANK)
  }

  return(object)
}

##################################################



