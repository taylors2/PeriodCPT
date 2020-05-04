##Wrapper function for PeriodCPT
PeriodCPT <- function(data,
                      distribution = c("bern","pois","norm","norm","mean","var"),
                      ...){
  distribution = match.arg(distribution)
  if(missing(data)) stop("data is missing.")
  method = paste0("PeriodCPT:::PeriodCPT.",distribution)
  PeriodCPT.FN <- get(paste0("PeriodCPT:::PeriodCPT.",distribution))
  ans <- PeriodCPT.FN(data, distribution, ...)
  return(ans)
}

##The main PeriodCPT function that calls the C code
PeriodCPT.main(pcpt){
  ##Input check
  ##Create pcpt class object
  ##Set inital values
  ##Perform MCMC (call C code)
  ##Format C output and put into object
}



#PeriodCPT.bern
#PeriodCPT.pois
#PeriodCPT.norm
#PeriodCPT.mean
#PeriodCPT.var

if(missing(inits)){
  if(n.chains(ans)==2){
    tau1 <- NULL
    tau2 <- seq(from = 0, by= minseglent(ans), length = npcpts.max(ans)) ##???
    inits.pcpt <- list(tau1, tau2)
  }else{
    ##sample from prior
    inits.pcpt <- list()
    for(i in 1:n.chains(ans)){
      ???
    }
  }
}else if(is.function(inits)){
  init.pcpts <- list()
  for(i in 1:n.chains(ans)){
    inits.pcpts[[i]] <- inits(...)
  }
}else if(is.list(inits)){
  if(length(inits) > n.chains(ans)) n.chains(ans) <- length(inits)
  init.pcpts <- inits[1:n.chains(ans)]
}else if(n.chains(ans)==1 && is.numeric(inits)){

}else if(missing(inits)){
  stop("Invalid inits.")
}


for(i in 1:n.chains(ans)){
  if(is.function(inits)){
    inits.pcpts[i] <- inits()
  }else if(is.list(inits)){
    if
    inits.pcpts[i] <- inits[i]
  }
}

