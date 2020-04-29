##################################
## Evaluate Night-time Interval ##
##################################

SS_by_time <- function(DRAW){
  if(!any("Freq"==colnames(DRAW$TAU))){
    TAU = SummariseOutput(DRAW)
  }else{
    TAU = DRAW$TAU
  }
  N <- attributes(DRAW)$N
  Mmax = max(TAU[,1])
  nSS <- (ncol(TAU)-2)/Mmax - 1 #(M,Freq,tau1--tauM,SS(1)_1--SS(1)_M,...,SS(nSS)_1--SS(nSS)_M)
  Stats <- array(NA,dim=c(nrow(TAU),N,nSS))
  time <- 1:N
  for(i in 1:nrow(TAU)){
    if(TAU[i,1]==1){
      for(s in 1:nSS) Stats[i,,s] <- TAU[i,2+s*Mmax+1]
    }else{
      L <- time<=TAU[i,3] | time>TAU[i,2+TAU[i,1]]
      for(s in 1:nSS) Stats[i,L,s] <- TAU[i,2+s*Mmax+1]
      for(j in 2:TAU[i,1]){
        L <- time <= TAU[i,2+j] & time>TAU[i,1+j]
        for(s in 1:nSS) Stats[i,L,s] <- TAU[i,2+s*Mmax+j]
      }
    }
  }
  STAT <- matrix(NA,nrow=nrow(TAU),ncol=nSS*N + 1)
  for(s in 1:nSS) STAT[,(s-1)*N+(1:N)] <- Stats[,,s]
  STAT[,ncol(STAT)] <- TAU[,2]/sum(TAU[,2])
  ssname <- unlist(strsplit(colnames(TAU)[(1:nSS)*Mmax+3],"1"))
  ssnameN <- paste0(rep(ssname,each=N),gsub(" ","0",format(rep(1:N,nSS),
    zero.print=TRUE)))
  colnames(STAT) <- c(ssnameN,"wt")
  return(STAT)
}



FindBounds <- function(DRAW, prob=0.5){
  if(!is.numeric(prob)) stop("Invalid prob argument.")
  if(length(prob)!=1 | any(prob<=0) | any(prob>=1)) stop("Invalid prob argument.")
  SS <- SS_by_time(DRAW)
  N <- attributes(DRAW)$N
  nSS <- (ncol(SS)-1)/N
  if(attributes(DRAW)$DIST=="Poisson"){
    tmp <- SS[,ncol(SS)]*SS[,1:N]/SS[,N+(1:N)]
    if(!is.matrix(tmp)){tmp <- matrix(tmp,nrow=1)}
    Expectation <- apply(tmp,2,sum)
    fbound <- function(x,sstats,N,prob){
      As <- sstats[,1:N]
      Bs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      tmp2 <- Wt*pgamma(x,As,Bs)
      if(!is.matrix(tmp2)) tmp2 <- matrix(tmp2,nrow=1)
      return(as.numeric(apply(tmp2,2,sum)>prob))
    }
    f <- function(x,sstats,N){
      As <- sstats[,1:N]
      Bs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      return(sum(Wt*dgamma(x,As,Bs)))
    }
    LOWER <- 0
    UPPER <- max(DRAW$data)
  }else if(attributes(DRAW)$DIST=="Binomial"){
    tmp <- SS[,ncol(SS)]*SS[,1:N]/(SS[,1:N] + SS[,N+(1:N)])
    if(!is.matrix(tmp)){tmp <- matrix(tmp,nrow=1)}
    Expectation <- apply(tmp,2,sum)
    fbound <- function(x,sstats,N,prob){
      As <- sstats[,1:N]
      Bs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      tmp2 <- Wt*pbeta(x,As,Bs)
      if(!is.matrix(tmp2)) tmp2 <- matrix(tmp2,nrow=1)
      return(as.numeric(apply(tmp2,2,sum)>prob))
    }
    f <- function(x,sstats,N){
      As <- sstats[,1:N]
      Bs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      return(sum(Wt*dbeta(x,As,Bs)))
    }
    LOWER <- 0
    UPPER <- max(DRAW$data)
  }else if(attributes(DRAW)$DIST=="Normal"){
    ##Evaluate bounds based on variance parameter
    tmp <- SS[,ncol(SS)]*SS[,3*N+(1:N)]/pmax(SS[,2*N+(1:N)] - 1,.Machine$double.eps)
    if(!is.matrix(tmp)){tmp <- matrix(tmp,nrow=1)}
    Expectation <- apply(tmp,2,sum)
    fbound <- function(x,sstats,N,prob){
      As <- sstats[,2*N+(1:N)]
      Bs <- sstats[,3*N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      tmp2 <- Wt*(1-pgamma(1/x,As,Bs))
      if(!is.matrix(tmp2)) tmp2 <- matrix(tmp2,nrow=1)
      return(as.numeric(apply(tmp2,2,sum)>prob))
    }
    f <- function(x,sstats,N){
      As <- sstats[,2*N+(1:N)]
      Bs <- sstats[,3*N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      return(sum(Wt*dgamma(1/x,As,Bs)/x^2))
    }
    LOWER <- .Machine$double.eps
    UPPER <- max(DRAW$data)
  }else if(attributes(DRAW)$DIST=="Mean"){
    tmp <- SS[,ncol(SS)]*SS[,1:N]
    if(!is.matrix(tmp)){tmp <- matrix(tmp,nrow=1)}
    Expectation <- apply(tmp,2,sum)
    fbound <- function(x,sstats,N,prob){
      Ms <- sstats[,1:N]
      Cs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      tmp2 <- Wt*pnorm(x,mean=Ms,sd=sqrt(Cs))
      if(!is.matrix(tmp2)) tmp2 <- matrix(tmp2,nrow=1)
      return(as.numeric(apply(tmp2,2,sum)>prob))
    }
    f <- function(x,sstats,N){
      Ms <- sstats[,1:N]
      Cs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      return(sum(Wt*dnorm(x,mean=Ms,sd=sqrt(Cs))))
    }
    LOWER <- min(DRAW$data)
    UPPER <- max(DRAW$data)
  }else if(attributes(DRAW)$DIST=="Var"){
    tmp <- SS[,ncol(SS)]*SS[,N+(1:N)]/pmax(SS[,(1:N)] - 1,.Machine$double.eps)
    if(!is.matrix(tmp)){tmp <- matrix(tmp,nrow=1)}
    Expectation <- apply(tmp,2,sum)
    fbound <- function(x,sstats,N,prob){
      As <- sstats[,1:N]
      Bs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      tmp2 <- Wt*(1-pgamma(1/x,As,Bs))
      if(!is.matrix(tmp2)) tmp2 <- matrix(tmp2,nrow=1)
      return(as.numeric(apply(tmp2,2,sum)>prob))
    }
    f <- function(x,sstats,N){
      As <- sstats[,1:N]
      Bs <- sstats[,N+(1:N)]
      Wt <- sstats[,ncol(sstats)]
      return(sum(Wt*dgamma(1/x,As,Bs)/x^2))
    }
    LOWER <- .Machine$double.eps
    UPPER <- max(DRAW$data)
  }else{
    stop("Distribution assumption not recognised or supported.")
  }
  E <- c(LOWER,sort(unique(round(Expectation,3))),UPPER)  ####
  Em <- 0.5*(E[-1] + E[-length(E)])   ##??##
  
  Eall <- sort(c(E,Em))
  fEall <- sapply(Eall,f,sstats=SS,N=N)
  L1 <- diff(fEall[-1])>0
  L2 <- diff(fEall[length(Eall):1][-1])[(length(Eall)-2):1]>0
  L <- c(TRUE,L1 & L2,TRUE)
  MINish <- Eall[L]
  
  BOUND <- sapply(MINish,fbound,sstats=SS,N=N, prob=prob)
  BOUND <- cbind(0,BOUND,1)
  L <- rep(TRUE,ncol(BOUND))
  for(i in 2:ncol(BOUND))  L[i] <- !all(BOUND[,i-1]==BOUND[,i])
  BOUND <- BOUND[,L]
  rownames(BOUND) <- attributes(DRAW)$readtime
  return(BOUND)
}

NightBounds <- function(DRAW,NightAt="03:00",prob=0.5){
  if(!is.numeric(prob)) stop("Invalid prob argument.")
  if(length(prob)!=1 | any(prob<=0) | any(prob>=1)) stop("Invalid prob argument.")
  Bounds <- FindBounds(DRAW, prob)
  if(ncol(Bounds)==2){ ##Cannot identify any intermediate bound levels
    return(c(NA,NA,NA))
  }
  N <- attributes(DRAW)$N
  wd <- 24*60/N
  if(length(NightAt)!=1) stop("Can only accept a single time for NightAt.")
  if(is.character(NightAt)){
    time_txt <- attributes(DRAW)$readtime
    i <- which(time_txt==NightAt)
    if(length(i)!=1) stop("Cannot recognise NightAt.")
  }else{
    i <- NightAt
    if(!(i%in%(1:N))) stop("Cannot recognise NightAt.")
  }
  b <- which(Bounds[i,]==1)[1]
  if(b==ncol(Bounds)){
    return(c(NA,NA,NA)) #NightAt happens at the highest level!!!
  }
  
  pts_st <- pts_ed <- NULL
  IN <- FALSE
  for(k in 1:(2*N)){
    kthis <- (k-1)%%N + 1
    kprev <- (k-2)%%N + 1
    if(Bounds[kthis,b]==1 && Bounds[kprev,b]==0 && !IN){
      if(length(pts_st)>0){ if(pts_st[1]==kprev) break }
      pts_st <- c(pts_st,kprev)
      IN <- !IN
    }else if(Bounds[kthis,b]==0 && Bounds[kprev,b]==1 && IN){
      pts_ed <- c(pts_ed,kprev)
      IN <- !IN
    }
  }
  pts_st <- pts_st - N*as.numeric(pts_st>pts_ed)
  k <- pts_st < i & i<=pts_ed
  ngt_st <- pts_st[k]
  ngt_ed <- pts_ed[k]
  if(length(ngt_st)==0 | length(ngt_ed)==0){
    #NightAt is at highest level of activity in period
    return(c(NA,NA,NA))
  }
  if(ngt_st<=0) ngt_st <- ngt_st + N
  ngt <- c(ngt_st,ngt_ed)
  
  ii <- i%%N + 1
  bb <- b
  looped <- found <- FALSE
  while(!looped && !found){
    if(ii==i){
      looped <- TRUE
    }else if(any(Bounds[ii,(1:(bb-1))]==1)){
      found <- TRUE
    }else if(Bounds[ii,bb]==1){
      ii <- ii%%N + 1
    }else if(Bounds[ii,bb]==0){
      bb <- bb+1
    }
  }
  late <- ii-1 + N*as.numeric(ii==1)
  out <- time_txt[c(ngt,late)]
  if(late==ngt[1]) out[3] <- NA
  return(out)
} #ngt_st, ngt_end, (~end morn)

TauProb <- function(tau=NULL, DRAW){
  if(!any("Freq"==colnames(DRAW$TAU))){
    TAU = SummariseOutput(DRAW)
  }else{
    TAU = DRAW$TAU
  }
  if(is.null(tau)){
    count <- DRAW$TAU[DRAW$TAU[,1]==1,2]
  }else{
    count <- 0
    for(i in 1:nrow(DRAW$TAU)){
      L <- all(tau %in% as.numeric(DRAW$TAU[i,2+(1:DRAW$TAU[i,1])]))
      if(L & DRAW$TAU[i,1]!=1){
        count <- count + DRAW$TAU[i,2]
      }
    }
  }
  return(as.numeric(count/sum(DRAW$TAU[,2])))
}

NightProb <- function(ngt,DRAW){
  ngt <- ngt[1:2]
  if(anyNA(ngt)){
    return(NA)
  }else if(any(is.character(ngt))){
    nt1 <- which(attributes(DRAW)$readtime==ngt[1])
    nt2 <- which(attributes(DRAW)$readtime==ngt[2])
    if(length(nt1)!=1) return(NA)
    if(length(nt2)!=1) return(NA)
    if(nt1==nt2) return(NA)
    tau <- c(nt1,nt2)
  }else if(is.numeric(ngt)){
    if(any(floor(ngt)!=ngt)) return(NA)
    if(ngt[1]==ngt[2]) return(NA)
    if(any(ngt<=0) | any(ngt>attributes(DRAW)$N)) return(NA)
    tau <- ngt
  }else{
    stop("Invalid 'ngt' argument.")
  }
  return(TauProb(tau,DRAW))
}


NightIntervals <- function(DRAW,NightAt="03:00",prob=0.5){
  if(!is.numeric(prob)) stop("Invalid prob argument.")
  if(length(prob)!=1 | any(prob<=0) | any(prob>=1)) stop("Invalid prob argument.")

  #Calculate the interval boundary matrix
  Bounds <- FindBounds(DRAW, prob)
  
  if (ncol(Bounds) == 2) {
    #If no temporal variation exists
    return(list(start = NA, end = NA))
  }
  N <- attributes(DRAW)$N
  wd <- 24 * 60/N
  if (length(NightAt) != 1){
    stop("Can only accept a single time for NightAt.")
  }
  
  #Work out the numerical value for NightAt.
  if (is.character(NightAt)) {
    time_txt <- attributes(DRAW)$readtime
    i <- which(time_txt == NightAt)[1]
    if (length(i) != 1) 
      stop("Cannot recognise NightAt.")
  } else {
    i <- NightAt
    if (!(i %in% (1:N))) 
      stop("Cannot recognise NightAt.")
  }
  
  #Find period starts and ends for where Bounds == 1
  STEDv2 <- NULL
  for(j in 2:(ncol(Bounds)-1)){
    #starts: where Bounds has 0->1 transition
    st <- unname(which(diff(c(Bounds[,j],Bounds[1,j]))==1))
    ed <- st
    LOOP <- rep(0,length(st))
    for(m in 1:length(ed)){
      #Keep searching until 1->0 transition as been found
      d <- c(Bounds[,j],Bounds[1,j])
      while(d[st[m]]!=d[ed[m]+1+LOOP[m]*N]){
        ed[m] <- ed[m]+1
        if(LOOP[m]==0 & ed[m]>N) LOOP[m] <- -1 #Interval spans the N->1 period boundary
      }
    }
    #Append found intervals to storage matrix
    STEDv2 <- rbind(STEDv2,cbind(j,-LOOP,st,ed+LOOP*N,ed-st))
  }
  STEDv2 <- rbind(STEDv2,cbind(ncol(Bounds),1,1,1,N)) ##Include full N period
  colnames(STEDv2) <- c("lev", "loop", "st", "ed", "len")
  
  #Identify where intevals contain NightAt
  Lv2 <- (STEDv2[,3] < i & (STEDv2[,4]+STEDv2[,2]*N) > i) | ((STEDv2[,3]-STEDv2[,2]*N) < i & STEDv2[,4] > i)
  if(sum(Lv2)>1){
    #Eliminate repeated intervals
    tmp <- unique(apply(STEDv2[Lv2,-1],1,paste0,collapse = "_"))
    tmp2 <- matrix(as.numeric(sapply(tmp,function(a){strsplit(a,"_")[[1]]})),nrow=length(tmp),byrow=TRUE)
    tmp2 <- tmp2[order(tmp2[,ncol(tmp2)]),]
    STout <- tmp2[,2]
    EDout <- tmp2[,3]
  }else{
    #Should be only the full N period!
    STout <- STEDv2[Lv2,3]
    EDout <- STEDv2[Lv2,4]
  }
  #Return unique period starts and ends that contain provided NightAt
  return(list(start = STout, end = EDout))

}


