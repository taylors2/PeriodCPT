##Wrapper function for PeriodCPT
PeriodCPT <- function(data,
                      distribution = c("bern","pois","norm","norm","mean","var"),
                      ...){
  distribution = match.arg(distribution)
  if(missing(data)) stop("data is missing.")
  ans <- eval(call(paste0("PeriodCPT.",distribution,"(data = data, ...)")))
  return(ans)
}

##The main PeriodCPT function that calls the C code
PeriodCPT.main <- function(object){

  ##Input check
  ##Create pcpt class object
  ##Set inital values
  ##Perform MCMC (call C code)
  ##Format C output and put into object

  draw <- list(C.chain.output = rep(0, n.chains(object)*npcpts.max(object)*n.iter(object)))

  len.chain.samples <- npcpts.max(object)*n.iter(object)
  for(i in 1:n.chains(object)){
    MCMC.chain(object, i) <- populate.chain(object,
      draw$C.chain.output[(i-1)*len.chain.samples + (1:len.chain.samples)],
      blank = -1)
  }

}

##################################################



