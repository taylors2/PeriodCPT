##Wrapper function for PeriodCPT
PeriodCPT <- function(data,
                      distribution = c("bern","pois","norm","mean","var"),
                      ...){
  distribution = match.arg(distribution)
  if(missing(data)) stop("Data is missing.")
  ans <- eval(parse(text=paste0("PeriodCPT.",distribution,"(data = data, ...)")))
  return(ans)
}

##The main Pernd distibution specific functioiodCPT function that calls the C code
PeriodCPT.main <- function(object){

  ##--Input check--
  if(class(object) != "pcpt") stop("Unexpected class of `object`.")
  #All other formal checks should be done by class specific operators

  ##--Format inputs to C--

  Eval_Mode <- TRUE  ##Hard code the evaluation of fits and mode summaries
  BLANK <- -1        ##Internal holding character for blanks

  offset <- start(data.set(object))[2]    ##Evaluate time axis
  time <- seq(from = offset, by = 1, length.out = length(data.set(object)))
  time2C <- ((time - 1) %% periodlength(object)) + 1
  inits2C <- rep(BLANK, n.chains(object)*npcpts.max(object))

  for(i in 1:n.chains(object)){ ##Vectorise inital values
    tau <- MCMC.inits(object)[[as.character(i)]]
    inits2C[(npcpts.max(object))*(i-1) + (1:length(tau))] <- tau
  }

  ##--Run C code--
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
                                n.iter(object) * npcpts.max(object)),

             Eval_Mode   = as.integer(Eval_Mode),
             mode_pcpt   = vector("integer", length = npcpts.max(object)),
             mode_params = vector("numeric", length = npcpts.max(object) * nsegparam(object)),
             nsegparams  = as.integer(nsegparam(object)),
             fits        = vector("numeric", length = 5),
             nfits       = as.integer(5)
  )

  ##Error checking
  PeriodCPT_ErrorBank(draw$error)

  ##Format C output and put into object
  len.chain.samples <- npcpts.max(object)*n.iter(object)
  for(i in 1:n.chains(object)){
    result(object, i) <- populate.chain(object,
                                            draw$draw[(i-1)*len.chain.samples + (1:len.chain.samples)],
                                            blank = BLANK)
  }

  if(Eval_Mode){
    pcpt.mode(object) <- draw$mode_pcpt[draw$mode_pcpt != BLANK]
    mode_param <- matrix(draw$mode_param[draw$mode_param != BLANK], nrow=nsegparam(object))
    rownames(mode_param) <- paste0("Param", 1:nsegparam(object))
    colnames(mode_param) <- paste0("Seg", 1:ncol(mode_param))
    param.mode(object) <- mode_param
    fit(object) <- draw$fit
  }

  return(object)
}

##################################################

PeriodCPT_ErrorBank <- function(errorcode){
  if(errorcode == 1){
    stop("Mprior distribution not recognised, issue with not recognising pcpt number prior distribution.")
  }else if(errorcode == 2){
    stop("Sampling distribution not recognised, issue with not recognising sampling distribution.")
  }else if(errorcode == 3){
    stop("Summary Statistic function not identified, issue with not recognising sampling distribution.")
  }else if(errorcode == 4){
    stop("Sufficient Statistic function not identified, issue with not recognising sampling distribution.")
  }else if(errorcode == 5){
    stop("Segment parameter mode function not identified, issue with not recognising sampling distribution.")
  }else if(errorcode == 6){
    stop("Fit metric calculation function not identified, issue with not recognising sampling distribution.")
  }else if(errorcode == 101){
    warning("Segment param mode is not unique. Consider changing prior hyper-params.")
  }else if(errorcode == 102){
    warning("Small sample size on segment for reliable fit metric. Consider increasing minimum segment length.")
  }
}





