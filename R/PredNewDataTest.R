pred_given_tau <- function(newdata, newtime, TAU, N, hyp, Mmax){
  ##Calculates predictive probabilty of 'newdata' at 'newtime' for a given 
  ##  cpt vector ('TAU') with associated sufficient statistics
  
  ##NB: this function is to only be used within PredProb()
  
  M <- TAU[1]
  segj <- rep(1,N)      #seg id as function of period time
  if(M>1) for(j in 2:M) segj[(TAU[1+j]+1):(TAU[2+j])] <- j
  jj <- segj[newtime]   #seg ids for newtime interval
  out <- 0
  for(j in unique(jj)){  ##loop over unique segments within newtime interval
    k <- which(jj==j)
    nj <- sum(k)
    ystar <- newdata[k]
    tstar <- newtime[k]
    Aj <- TAU[Mmax+2+j]
    Bj <- TAU[2*Mmax+2+j]
    mu <- hyp[1]
    ##Conditional Predictive prob: \int P(Y* = y*| \sigma^2, \tau) \pi(\sigma^2|\tau,y) \pi(\tau|y) d(\sigma^2)
    lprobseg <- lgamma(Aj+nj/2) - lgamma(Aj) - 0.5*nj*log(2*pi*Aj)
    lprobseg <- lprobseg - (Aj+nj/2)*log(1 + sum((ystar-mu)^2)/(2*Bj))
    out <- out + lprobseg
  }
  return(unname(out))
}

PredProb <- function(newdata,newtime,DRAW){
  #Calculates how predictible 'newdata' at 'newtime' is given the RJMCMC 'DRAW'
  if(attributes(DRAW)$DIST!="Var") stop("Distribution assumption not supported.")
  if (!any("Freq" == colnames(DRAW$TAU))) {
    TAU = SummariseOutput(DRAW)
  }else {
    TAU = DRAW$TAU
  }
  N <- attributes(DRAW)$N
  Mmax <- max(DRAW$TAU[,1])
  hyp <- attributes(DRAW)$hyp
  #Apply pred_given_tau() to each cpt vector, i.e. each row of DRAW$TAU
  lpredprob_tau <- apply(TAU,1,pred_given_tau,N=N,hyp=hyp,newdata=newdata,newtime=newtime,Mmax=Mmax)
  
  #Marginalise conditional predictive probs over cpt vector posteior
  Wt <- TAU[,2]/sum(TAU[,2])  
  lmax <- max(lpredprob_tau)
  out <- sum(exp(lpredprob_tau-lmax)*Wt)*exp(lmax)
  
  #out = P(Y* = y* | y)
  return(out)
}


SampleNULL <- function(newtime,DRAW,B=1000){
  #For given 'newtime', what is the predictive distribution P(Y*=y*|y) under the 
  #  null hypothesis that the 'newdata' is equivalently distributed as what has been observed
  #Perform a bootstrap sample of size 'B' from available data
  
  #Input checking
  if(length(B)!=1 | !is.numeric(B)) stop("Invalid B.")
  if(B%%1!=0 | B<1) stop("Invaild B.")
  if(!is.list(DRAW)) stop("Invalid DRAW.")
  if(!all(c("TAU","data","time")%in%names(DRAW))) stop("Invalid DRAW.")
  att <- c("names","DIST","hyp","gamma","delta","l","N","control","name","date","readtime","Message")
  if(!all(att%in%names(attributes(DRAW)))) stop("Invalid DRAW.")
  if(!any("Freq" == colnames(DRAW$TAU))){
    TAU = SummariseOutput(DRAW)
  }else {
    TAU = DRAW$TAU
  }
  Mmax <- max(DRAW$TAU[,1])
  N <- attributes(DRAW)$N
  if(missing(newtime)) stop("newtime is missing")
  if(is.character(newtime)){
    tmp <- rep(NA,length(newtime))
    for(i in 1:length(newtime)){
      k <- which(newtime[i]==attributes(DRAW)$readtime)
      if(length(k)!=1) stop("Invalid newtime.")
      tmp[i] <- k
    }
    newtime <- k
  }else{
    if(!is.numeric(newtime)) stop("Invailid newtime,")
    if(!all(newtime%%1==0) | any(newtime<1) | any(newtime>N)) stop("Invaild newtime.")
  }
  
  #sample cpt vector from posterior
  t <- sample(1:nrow(DRAW$TAU),size = B, prob=DRAW$TAU[,2],replace=TRUE)
  m <- unname(DRAW$TAU[t,1])
  
  #generate 'newdata' under emperical distribution
  nulldata <- matrix(NA,nrow=B,ncol=length(newtime))
  for(b in 1:B){
    if(m[b]==1){
      #if sampled single segment model, sample data from any time point
      nulldata[b,] <- sample(DRAW$data,size=length(newtime),replace=FALSE)
    }else{
      #Else, sample data anywhewre from appropriate segment
      segid <- rep(1,N)
      for(j in 2:m[b]) segid[(DRAW$TAU[t[b],1+j]+1):(DRAW$TAU[t[b],2+j])] <- j
      for(i in 1:length(newtime)){
        times <- (1:N)[segid==segid[newtime[i]]]
        nulldata[b,i] <- sample(DRAW$data[DRAW$time %in% times],size=1)
      }
    }
  }
  #apply PredProb to each row of nulldata distribution of statistic under null hypothesis
  nulldist <- apply(nulldata,1,PredProb,newtime=newtime,DRAW=DRAW)
  return(nulldist)
}

PredLogOdd <- function(newdata,newtime,DRAW,B=1000){
  ##Evaluate predictive under null hypothesis & for 'newdata' at 'newtime'
  ##B = size of boostrap sample under null hypothesis
  ##Calculate the log-odds ratio:
  ##  LO = log(P(Y*=y*| y)) - log(P(Y*=y'| y, y'~H0))
  
  #If high, then 'newdata' was highly predictable (near to mu, so lower variance??)
  #If low, then 'newdata' was very difficult to predict (far from mu, so high variance??)
  #If log-odd ~= 0, then little evidence to distinguish from NULL (same variance??)
  
  nulldist <- SampleNULL(newtime = newtime, DRAW = DRAW, B=B)
  alt <- PredProb(newdata=newdata,newtime = newtime,DRAW = DRAW)
  PredlogOdd <- log(alt)-log(nulldist)
  return(PredlogOdd)
}
