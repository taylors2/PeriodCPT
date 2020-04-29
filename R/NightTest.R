################################
##  Night-time activity test  ##
################################

nightTestDist <- function(DRAW, ngt, tol=1e-08){
  if(!any("Freq"==colnames(DRAW$TAU))){
    TAU = SummariseOutput(DRAW)
  }else{
    TAU = DRAW$TAU
  }
  if(missing(ngt)) ngt <- NightBounds(DRAW)
  N <- attributes(DRAW)$N
  time_txt <- attributes(DRAW)$readtime
  ngt_st <- which(ngt[1] == time_txt)
  ngt_ed <- which(ngt[2] == time_txt)
  if(ngt_st==N){
    ngti <- 1:ngt_ed
  }else if(ngt_ed>ngt_st){
    ngti <- (ngt_st+1):ngt_ed
  }else{
    ngti <- c((ngt_st+1):N,1:ngt_ed)
  }
  data <- DRAW$data[DRAW$time %in%ngti]
  binsPerH <- N/24
  nsums <- sum(table(sapply(time_txt[ngti],function(a){strsplit(a,":")[[1]][1]}))==binsPerH)
  if(attributes(DRAW)$DIST=="Poisson"){
    ALPHA <- attributes(DRAW)$hyp[1] + sum(data)   
    BETA  <- attributes(DRAW)$hyp[2] + length(data) 
    cmf_max <- function(rate, m,binsPerH,nsums,ALPHA,BETA){
      return(ppois(m,binsPerH*rate)^nsums * dgamma(rate,ALPHA,BETA))
    }
    summax <- qnbinom(tol,size=ALPHA,prob=BETA/(BETA+length(ngti)),lower.tail = FALSE)
    cmfM <- NULL; cont <- TRUE; m <- 0
    while(cont){
      #Eval P[max <= m] by marginalising over probability
      i <- integrate(f=cmf_max,lower=0,upper=Inf,m=m,binsPerH=binsPerH,
                     nsums=nsums,ALPHA=ALPHA,BETA=BETA,rel.tol = tol)
      cmfM <- c(cmfM,i$value)
      if(cmfM[m+1]>=1-tol && m>=summax) cont <- FALSE
      m <- m + 1
    }
    cmfM[length(cmfM)] <- 1
    PMF <- diff(c(0,cmfM))
    PMF[PMF<tol] <- 0
    summax <- max(summax,length(PMF)-1)
    PMFtotal <- dnbinom(0:summax,size=ALPHA,prob=BETA/(BETA+length(ngti)))
  }else if(attributes(DRAW)$DIST=="Binomial"){
    
    pmf_joint_longest_p <- function(p,T,Pmaster,vec,result,ALPHA,BETA,k=NA){
      #Joint probability for longest seq=k and posterior prob = p
      OUT <- rep(NA,length(p))
      for(i in 1:length(p)){
        P <- Pmaster;P[P==2] <- 1-p[i];P[P==3] <- p[i]
        p_longest <- vec
        for(t in 1:T) p_longest <- p_longest %*% P #Steps for each bin in night
        p_longest <- as.numeric(p_longest %*% result) #Collapse states to max length statistic
        out <- p_longest*dbeta(p[i],ALPHA,BETA)
        OUT[i] <- out[k+1]
      }
      return(OUT)
    }
    ALPHA <- attributes(DRAW)$hyp[1] + sum(data)   
    BETA  <- attributes(DRAW)$hyp[2] + length(data) - sum(data)
    T <- length(ngti)
    
    PMF <- rep(0,T+1)                          #Prob(longest seq == m) -- to be evaluated
    curCMF <- 0
    for(k in 0:(T-1)){    #For each possible maximum length (up to period length)
      NN <- min(k+1,T)
      grid <- matrix(NA,ncol=2,nrow=0.5*NN*(NN+1)+1) #Define the set of states
      M <- 0; l <- 0
      for(i in 1:nrow(grid)){
        grid[i,] <- c(l,M)
        l <- ifelse(l==0,M+1,l-1)
        M <- ifelse(l>M,l,M)
      }
      Pmaster <- matrix(0,ncol=nrow(grid),nrow=nrow(grid))  #Evaluate the transision matrix
      for(i in 1:(nrow(Pmaster)-1)){
        Pmaster[i,which(grid[,1]==0 & grid[,2]==grid[i,2])] <- 2
        Pmaster[i,which(grid[,1]==(1+grid[i,1]) & 
                          grid[,2]==(grid[i,2] + as.numeric(grid[i,1]==grid[i,2])))] <- 3
      }
      Pmaster[nrow(Pmaster),ncol(Pmaster)] <- 1   #last scenario is an absorbant state
      vec <- cbind(1,matrix(0,nrow=1,ncol=nrow(Pmaster)-1)) #initial state vector (start at [0,0])
      result <- matrix(0,nrow=ncol(Pmaster),ncol=1+NN)    #marginalise final state vector by seq. max
      for(M in 0:(ncol(result)-1)) result[which(grid[,2]==M),M+1] <- 1
      
      #Marginalise probability wrt posterior distribution.
      i <- integrate(f=pmf_joint_longest_p,lower=0,upper=1,T=T,ALPHA=ALPHA,BETA=BETA,
                     Pmaster=Pmaster,vec=vec,result=result, k=k,rel.tol=tol)
      PMF[k+1] <- i$value
      curCMF <- curCMF + i$value
      if(curCMF >= 1-tol){
        PMF[k+1] <- 1-sum(PMF[1:k])
        break
      }
      if((k+1)==T){
        PMF[T+1] <- 1 - curCMF
      }
    }
    PMFtotal <- dbbinom(0:T,size=T,alpha=ALPHA,beta=BETA)
  }else{
    stop("Distributional assumption not recognised.")
  }
  
  PMF <- data.frame(value=0:(length(PMF)-1), maxPMF=PMF, sumPMF=PMFtotal)
  attr(PMF,"ngt") <- ngt[1:2]
  attr(PMF,"ngti") <- ngti
  attr(PMF,"binsPerH") <- binsPerH
  attr(PMF,"nsums") <- nsums
  attr(PMF,"ALPHA") <- ALPHA
  attr(PMF,"BETA") <- BETA
  attr(PMF,"DIST") <- attributes(DRAW)$DIST
  return(PMF)
  
}

max_len <- function(x){
  #Evaluate longest sequence of 1
  m <- l <- x[1]
  for(i in 2:length(x)){
    l <- ifelse(x[i]==1,l+1,0)
    m <- max(m,l)
  }
  return(m)
}

nightTest <- function(DRAW, TestPMF, TestData, ...){
  if(missing(TestPMF)) TestPMF <- nightTestDist(DRAW, ...)
  #Collate data for complete nights
  ngti <- attributes(TestPMF)$ngti
  ngt_st <- ngti[1]
  ngt_ed <- ngti[length(ngti)]
  ngtDATA <- TestData[TestData$binid %in% ngti,]
  ist <- which(ngtDATA$binid == ngt_st)[1]
  ied <- max(which(ngtDATA$binid == ngt_ed))
  if(ist>ied) stop("No complete night in TestData.")
  DATE <- as.character(ngtDATA$Time[ngtDATA$binid==ngt_ed])
  DATE <- unname(sapply(DATE, function(a){strsplit(a," ")[[1]][1]}))
  if(ngt_st>ngt_ed) DATE <- DATE[-1]
  if(attributes(TestPMF)$DIST=="Poisson"){
    data <- matrix(ngtDATA$Count[ist:ied],nrow=length(ngti))
    time <- sapply(as.character(ngtDATA$Time),function(a){strsplit(a," ")[[1]][2]})
    time <- matrix(time[ist:ied],nrow=length(ngti))[,1]
    while(!grepl(":00",time[1])){
      data <- as.matrix(data[-1,]); time <- time[-1]
    }
    binsPerH <- attributes(TestPMF)$binsPerH
    nsums <- attributes(TestPMF)$nsums
    data <- as.matrix(data[1:(binsPerH*nsums),])
    probMl <- probSl <- probM <- M <- probS <- S <- rep(0,ncol(data))
    for(i in 1:ncol(data)){
      M[i] <- max(apply(matrix(data[,i],nrow=binsPerH,ncol=nsums),2,sum))
      if(M[i]==0){
        probM[i] <- 1
      }else if(M[i] < nrow(TestPMF)){
        probM[i] <- sum(TestPMF$maxPMF[-(1:M[i])])
      }
      if(M[i]==length(TestPMF$maxPMF)){
        probMl[i] <- 1
      }else{
        probMl[i] <- sum(TestPMF$maxPMF[(1:(1+M[i]))])
      }
      S[i] <- sum(data[,i])
      if(S[i] == 0){
        probS[i] <- 1
      }else if(S[i] < nrow(TestPMF)){
        probS[i] <- sum(TestPMF$sumPMF[-(1:S[i])])
      }
      if(S[i]==length(TestPMF$sumPMF)){
        probSl[i] <- 1
      }else{
        probSl[i] <- sum(TestPMF$sumPMF[(1:(1+S[i]))])
      }
    }
  }else if(attributes(TestPMF)$DIST=="Binomial"){
    data <- matrix(ngtDATA$Active[ist:ied],nrow=length(ngti))
    probMl <- probSl <- probM <- M <- probS <- S <- rep(0,ncol(data))
    for(i in 1:ncol(data)){
      M[i] <- max_len(data[,i])
      if(M[i]==0){
        probM[i] <- 1
      }else if(M[i] < nrow(TestPMF)){
        probM[i] <- sum(TestPMF$maxPMF[-(1:M[i])])
      }
      if(M[i]==length(TestPMF$maxPMF)){
        probMl[i] <- 1
      }else{
        probMl[i] <- sum(TestPMF$maxPMF[(1:(1+M[i]))])
      }
      S[i] <- sum(data[,i])
      if(S[i] == 0){
        probS[i] <- 1
      }else if(S[i] < nrow(TestPMF)){
        probS[i] <- sum(TestPMF$sumPMF[-(1:S[i])])
      }
      if(S[i]==length(TestPMF$sumPMF)){
        probSl[i] <- 1
      }else{
        probSl[i] <- sum(TestPMF$sumPMF[(1:(1+S[i]))])
      }
      
    }
  }else{
    stop("Distributional assumption not recognised.")
  }
  
  out <- data.frame(date = DATE, max=M, probGeqMax=round(probM,5), probLeqMax=round(probMl,5), 
                    sum=S, probGeqSum = round(probS,5), probLeqSum = round(probSl,5))
  return(out)
}

dbbinom <- function(x, size, alpha=1, beta=1, log=FALSE){
  out = lchoose(size, x)
  out = out + lbeta(alpha+x, size-x-beta)
  out = out - lbeta(alpha, beta)
  if(!log) out = exp(out)
  return(out)
}