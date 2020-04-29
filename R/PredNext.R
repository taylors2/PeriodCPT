#########################################################################################
## Evaluate predictive probability of next observation occuring within 'h' steps ahead ##
#########################################################################################

Eval_index <- function(tau, offset, N){
  tau = tau - offset
  j <- 1:length(tau)
  while(tau[1]<=0){
    tau <- c(tau[-1], (tau[1]-1)%%N +1)
    j <- c(j[-1],j[1])
  }
  if(length(tau)<=1){
    nn <- 1:N
    ii <- rep(j,N)
  }else{
    nn <- rep(0,N)
    ii <- rep(j[1],N)
    itmp <- seq(from=1,to=tau[1],by=1)
    nn[itmp] <- itmp
    ist <- length(itmp)
    for(i in 2:length(tau)){
      itmp <- 1:(tau[i]-tau[i-1])
      nn[ist + itmp] <- itmp
      ii[ist + itmp] <- j[i]
      ist <- ist + length(itmp)
    }
    if(ist<N){
      nn[-(1:ist)] <- tau[1]+(1:(N - tau[length(tau)]))
    }
  }
  return(cbind(ii,nn))
}

Pred_Prob_of_Next_Obs <- function(DRAW, offset=0){
  if(!any("Freq"==colnames(DRAW$TAU))){
    TAU = SummariseOutput(DRAW)
  }else{
    TAU = DRAW$TAU
  }
  N <- attributes(DRAW)$N
  time_txt <- attributes(DRAW)$readtime
  Mmax = max(TAU[,1])
  
  if(length(offset)>1) stop("Invalid offset.")
  if(is.numeric(offset)){
    if(!(offset %in% (0:N))) stop("Invalid offset.")  
  }else if(is.character(offset)){
    i <- which(time_txt==offset)[1]
    if(length(i)==0) stop("Invalid offset.")
    offset <- i
  }else{
    stop("Invalid offset.")
  }
  ii <- nn <- matrix(NA,ncol=N,nrow=nrow(TAU))
  for(it in 1:nrow(ii)){
    tmp <- Eval_index(tau=TAU[it,2+(1:TAU[it,1])], offset=offset, N=N)
    ii[it,] <- tmp[,1] ##which segment on for h steps ahead from offset
    nn[it,] <- tmp[,2] ##times been on segment ii[it,h] at h steps ahead from offset
  }
  if(attributes(DRAW)$DIST == "Poisson"){
    AA <- BB <- array(NA,dim=dim(ii))
    for(it in 1:nrow(ii)){
      AA[it,] <- TAU[it,2+Mmax+ii[it,]]
      BB[it,] <- TAU[it,2+2*Mmax+ii[it,]]
    }
    p_0 <- ((BB+nn-1)/(BB+nn))^AA
    p_NonSoFar <- t(apply(cbind(1,p_0[,-ncol(p_0)]),1,cumprod))
    p <- p_NonSoFar*(1-p_0)
    cumP <- t(apply(p,1,cumsum))
    out <- apply(cumP*TAU[,2],2,sum)/sum(TAU[,2])
  }else if(attributes(DRAW)$DIST == "Binomial"){
    AA <- BB <- array(NA,dim=dim(ii))
    for(it in 1:nrow(ii)){
      AA[it,] <- TAU[it,2+Mmax+ii[it,]]
      BB[it,] <- TAU[it,2+2*Mmax+ii[it,]]
    }
    p_0 <- (BB+nn-1)/(AA+BB+nn-1)
    #p_1 <- (AA)/(AA+BB+nn-1)
    p_NonSoFar <- t(apply(cbind(1,p_0[,-ncol(p_0)]),1,cumprod))
    p <- p_NonSoFar*(1-p_0)
    cumP <- t(apply(p,1,cumsum))
    out <- apply(cumP*TAU[,2],2,sum)/sum(TAU[,2])
  }else{
    stop("Distributional assumption not recognised.")
  }
  out <- as.table(out)
  if(offset==0 | offset==N){
    names(out) <- time_txt
  }else{
    names(out) <- c(time_txt[-(1:offset)],time_txt[1:offset])
  }
  return(out)
}
