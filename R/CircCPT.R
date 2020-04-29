#########################################
##  RJMCMC CODE TO ESTIMATE Circ CPTS  ##
#########################################

CircCPT <- function(data, time, DIST=c("Poisson","Binomial","Normal","Mean","Var"),
                              hyp=NULL, gamma=1, delta=1, l=4,
                              init1=NULL, init2=NULL,
                              control=list(), summary=TRUE, MODE = FALSE){
  
  DIST <- match.arg(DIST)
  iDIST <- which(DIST == c("Poisson","Binomial","Normal","Mean","Var"))-1
  if(anyNA(data) | anyNA(time)) stop("data and time connot contain missing values.")
  if(!is.numeric(data) | !is.numeric(time)) stop("data and time must be a numerical vector.")
  if(length(data)!=length(time)) stop("data and time must have the same length,")
  if(any(time != as.integer(time))) stop("time must be an integer.")
  if(any(time<1)) stop("time must be greater than 0.")
  N = max(time)
  Mmax = floor(N/l)

  if(is.null(hyp)){
    ##Default hyp
    if(DIST=="Binomial") hyp=c(1,1)
    if(DIST=="Poisson") hyp=c(1,1)
    if(DIST=="Normal") hyp=c(0,1,1,1)  #m,c,a,b
    if(DIST=="Mean") hyp=c(1,1,1) #sig2, m, c
    if(DIST=="Var") hyp=c(0,1,1)
  }
  chyp = hyp
  if(iDIST==0){##POISSON
    if(any(data<0) | any(data != as.integer(data))) stop("Invalid data for distribution.")
    if(length(chyp)!=2) stop("Invalid prior parameters.")
    if(!all(is.numeric(chyp))) stop("Invalid prior parameters.")
    if(!all(chyp>0)) stop("Invalid prior parameters.")
    nparamPERseg <- 1
  }else if(iDIST==1){##Binary
    if(!all(data %in% c(0,1))) stop("Invalid data for distribution.")
    if(length(chyp)!=2) stop("Invalid prior parameter.")
    if(!all(is.numeric(chyp))) stop("Invalid prior parameters.")
    if(!all(chyp>0)) stop("Invalid prior parameter.")
    nparamPERseg <- 1
  }else if(iDIST==2){ #Normal, both mean and variance/precision unknown
    ##Data can take any real value
    if(length(chyp)!=4) stop("Invalid prior parameter.")
    if(!all(is.numeric(chyp))) stop("Invalid prior parameters.")
    if(!all(chyp[-1]>0)) stop("Invalid prior parameter.")
    chyp[2] <- 1/hyp[2]  #convert variance scale to precision scale
    nparamPERseg <- 2
  }else if(iDIST==3){ # Mean (Normal, known variance/precision)
    ##Data can take any real value
    if(length(chyp)!=3) stop("Invalid prior parameter.")
    if(!all(is.numeric(chyp))) stop("Invalid prior parameters.")
    if(!all(chyp>0)) stop("Invalid prior parameter.")    
    chyp[1] <- 1/hyp[1]  #convert known variance to known precision
    chyp[3] <- 1/hyp[3]  #convert variance scale to precision scale
    nparamPERseg <- 1
  }else if(iDIST==4){ # Var (Normal, known mean)
    ##Data can take any real value
    if(length(chyp)!=3) stop("Invalid prior parameter.")
    if(!all(is.numeric(chyp))) stop("Invalid prior parameters.")
    if(!all(chyp[-1]>0)) stop("Invalid prior parameter.")
    nparamPERseg <- 1
  }else{
    stop("Invalid distribution definition.")
  }
  
  if(length(gamma) >1 | length(delta) >1 ) stop("Invalid prior parameter.")
  if(!is.numeric(gamma) | !is.numeric(delta) ) stop("Invalid prior parameter.")
  if(gamma<=0 | delta<=0 ) stop("Invalid prior parameter.")
  if(length(l)>1 | !is.numeric(l)) stop("l must be a single numerical value.")
  if(l!=as.integer(l) | l<1 | l>N) stop("invalid value for minimum segment length.")
  
  if(!is.null(init1)){
    if(!is.numeric(init1)) stop("Invalid init1 vector.")
    if(length(init1)>Mmax) stop("Invalid init1 vector.")
    if(length(init1)<=1) init1 = 1
    if(any(init1 != as.integer(init1))) stop("Invalid init1 vector.")
    init1 = sort(init1)
    if(init1[1]<1 | init1[length(init1)]>N) stop("Invalid init1 vector.")
    d = diff(c(init1,init1[1]+N))
    if(min(d)<l) stop("Invalid init1 vector.")
  }else{
    init1 <- 1
  }
  
  if(!is.null(init2)){
    if(!is.numeric(init2)) stop("Invalid init2 vector.")
    if(length(init2)>Mmax) stop("Invalid init2 vector.")
    if(length(init2)<=1) init2 = 1
    if(any(init2 != as.integer(init2))) stop("Invalid init2 vector.")
    init2 = sort(init2)
    if(init2[1]<1 | init2[length(init2)]>N) stop("Invalid init2 vector.")
    d = diff(c(init2,init2[1]+N))
    if(min(d)<l) stop("Invalid init2 vector.")
  }else{
    init2 <- seq(from = 1, to = l*Mmax, by = l)
  }
  
  control0=list(itstep=1000,Maxit=10,cachesize=50,typeI=0.05, burn=TRUE, progress=FALSE)
  if(!all(names(control) %in% names(control0)))
    warning("Unrecognised items in control list, these will be ignored.")
  if(!is.list(control)) stop("Invalid control list.")
  if("itstep" %in% names(control)){
    if(!is.numeric(control$itstep) | length(control$itstep)>1 ) stop("Invalid control list -- itstep.")
    if(control$itstep<1 | control$itstep!=as.integer(control$itstep)) stop("Invalid control list -- itstep.")
    control0$itstep = control$itstep
  }
  if("Maxit" %in% names(control)){  
    if(!is.null(control$Maxit)){
      if(!is.numeric(control$Maxit) | length(control$Maxit)>1 ) stop("Invalid control list -- Maxit.")
      if(control$Maxit<1 | control$Maxit!=as.integer(control$Maxit)) stop("Invalid control list -- Maxit.")
      control0$Maxit = control$Maxit
    }else{
      control0$Maxit = Inf
    }
  }
  if("cachesize"%in%names(control)){
    if(!is.numeric(control$cachesize) | length(control$cachesize)>1 ) stop("Invalid control list -- cachesize.")
    if(control$cachesize<1 | control$cachesize!=as.integer(control$cachesize)) stop("Invalid control list -- cachesize.")
    control0$cachesize = control$cachesize
  }
  if("typeI" %in% names(control)){
    if(!is.numeric(control$typeI) | length(control$typeI)>1) stop("Invalid control list -- typeI.")
    if(control$typeI<0 | control$typeI>1) stop("Invalid control list -- typeI.")
    control0$typeI = control$typeI    
  }
  if("burn" %in% names(control)){
    if(length(control$burn)>1 | !is.logical(control$burn) | anyNA(control$burn)) stop("Invalid control list -- burn.")
    control0$burn = control$burn
  }
  if("progress" %in% names(control)){
    if(length(control$progress)>1 | !is.logical(control$progress) | anyNA(control$progress)) stop("Invalid control list -- progress.")
    control0$progress = control$progress
  }
  
  if(length(summary)>1 | !is.logical(summary) | anyNA(summary)) stop("invalid summary flag.")
  
  if(length(MODE)!=1 | any(!is.logical(MODE))) stop("invalid MODE flag.")
  if(iDIST!=2) MODE <- FALSE
  
  ##FUDGE!!!
  if(all(data==data[1])) init2 = 1  ##FORCE m=1 fit if var(data)==0

  Maxit = control0$Maxit
  if(Maxit==Inf) Maxit = 0
  nout = 2*(Mmax+1)*control0$itstep
  
  #dyn.load("CircCPT_RJMCMC.so")
  draw = .C("CircCPT_RJMCMC",
            as.double(data),		                ##data
            as.integer(time),		                ##time
            as.integer(length(data)),	                ##n
            as.integer(N),	        	        ##N
            as.integer(iDIST),		                ##dist(0=pois, 1=binom)
            as.double(chyp),		                ##hypparam
            as.double(gamma),			        ##hypS
            as.double(delta),			        ##hypM
            as.integer(l),			        ##minseglen
            as.integer(control0$itstep),		##itstep
            as.integer(Maxit),		        	##Maxit
            as.double(control0$typeI),	                ##typeI
            as.integer(init1),		                ##init1
            as.integer(length(init1)),		      	##m1
            as.integer(init2),			        ##init2
            as.integer(length(init2)),			##m2
            as.integer(control0$progress),	   	##ProgressBar
            as.integer(control0$cachesize),             ##CACHEMAX
            as.integer(Mmax),                           ##Max model size
            as.integer(control0$burn),                  ##burn
            err=as.integer(0),			        ##err
            draw=vector("integer",nout),                ##draw
            ndraw=as.integer(nout),	                ##ndraw
            conv = as.integer(0),                 ##conv
            as.integer(nparamPERseg),              ##number of segment parameters per segment
            out_mmode = as.integer(0),             ##Mode estimate for the number of changepoints
            out_taumode = vector("integer", Mmax), ##Mode estimate for the changepoint positions
            out_segparammode = vector("double",Mmax*nparamPERseg),   ##Mode estimates for the segment parameters
            out_fits = vector("double", 4)            ##Measures of fit
            )		        
  #dyn.unload("CircCPT_RJMCMC.so")

  if(draw$err!=0) stop(paste0("Function exited early with error code ",draw$err))
  
  TAU <- matrix(NA,ncol=Mmax,nrow=2*(control0$itstep+1))
  i=1
  j=1
  k=1
  while(i<=draw$ndraw){
    if(draw$draw[i]==0){
      i=draw$ndraw
    }else if(draw$draw[i]<0){
      j=j+1
      k=1
    }else{
      TAU[j,k] = draw$draw[i]
      k=k+1
    }
    i=i+1
  }
  TAU <- as.matrix(TAU[,!apply(is.na(TAU),2,all)])
  colnames(TAU) <- paste0("tau",1:ncol(TAU))
  if(!summary){ 
    out <- list(TAU=TAU, data=data, time=time)
  }else{
    out <- list(TAU=SummariseOutput(TAU, data, time, DIST, hyp), data=data, time=time)
  }
  attr(out,"DIST")     <- DIST
  attr(out,"hyp")      <- hyp
  attr(out,"gamma")    <- gamma
  attr(out,"delta")    <- delta
  attr(out,"l")        <- l
  attr(out,"N")        <- N
  attr(out,"DIST")     <- DIST
  attr(out,"control")  <- control0
  attr(out,"name")     <- deparse(substitute(data))
  attr(out,"date")     <- Sys.time()
  attr(out,"readtime") <- format(seq.POSIXt(as.POSIXct("1970-01-01 00:00:00"),
                                            as.POSIXct("1970-01-01 23:59:59"),by = 24*3600/N),"%H:%M")
  if(draw$conv==0){
    warning("Not converged, hit maximum number of iterations.")
    attr(out,"Message") <- paste0("Not converged, hit maximum number of iteration.")
  }else{
    attr(out,"Message") <- paste0("Chain convergence test is satisfied.")
  }
  
  if(MODE){
    ModeM <- draw$out_mmode
    ModeSegParam <- matrix(draw$out_segparammode[1:(ModeM*nparamPERseg)],
                           nrow=nparamPERseg,ncol=ModeM,byrow=TRUE)
    rownames(ModeSegParam) <- paste0("Param",1:nparamPERseg)
    colnames(ModeSegParam) <- paste0("Seg",1:ModeM)
    ModeTau <- draw$out_tau[1:ModeM]
    names(ModeTau) <- paste0("Tau",1:ModeM)
    ModeFits <- draw$out_fits
    names(ModeFits) <- c("MJTmode","JTmode","MLmode","Lmode")
    out$MODE <- list(M = ModeM, SegParam = ModeSegParam, Tau=ModeTau, Fits=ModeFits)
  }else{
    out$MODE <- NULL
  }
  return(out)
}

SummariseOutput <- function(TAU, data, time, DIST, hyp){
  if(is.list(TAU)){
    data = TAU$data
    time = TAU$time
    DIST = attributes(TAU)$DIST
    hyp = attributes(TAU)$hyp
    TAU = TAU$TAU
  }
  cpttxt = apply(TAU,1,paste0,collapse="_")
  tab = table(cpttxt)
  tab = tab[order(tab,decreasing = TRUE)]
  Summary = matrix(NA,ncol=ncol(TAU),nrow=length(tab))
  for(i in 1:length(tab)) Summary[i,] = TAU[which(names(tab)[i]==cpttxt)[1],]
  Summary = cbind(apply(!is.na(Summary),1,sum),as.numeric(tab),Summary)
  colnames(Summary) = c("M","Freq",colnames(TAU))
  ##POSTERIOR Sufficient Statistics Parameters
  SS = EvalSS(Summary, data, time, DIST, hyp)
  return(cbind(Summary,SS))
}

EvalSS = function(summary, data, time, dist, hyp){
  count = sum = sq = matrix(NA,nrow=nrow(summary),ncol=ncol(summary)-2)
  for(i in 1:nrow(summary)){
    m = summary[i,1]
    if(m==1){
      count[i,1] = length(time)
      sum[i,1] = sum(data)
      sq[i,1] = sum(data^2)
    }else{
      L = time<=summary[i,3] | time>summary[i,2+m]
      count[i,1] = sum(L)
      sum[i,1] = sum(data[L])
      sq[i,1] = sum(data[L]^2)
      for(j in 2:m){
        L = time<=summary[i,j+2] & time>summary[i,j+1]
        count[i,j] = sum(L)
        sum[i,j] = sum(data[L])
        sq[i,j] = sum(data[L]^2)
      }
    }
  }
  
  if(dist=="Poisson"){ ##Poisson
    ALPHA = hyp[1] + sum
    BETA = hyp[2] + count
    colnames(ALPHA) = paste0("A",1:ncol(ALPHA))
    colnames(BETA) = paste0("B",1:ncol(ALPHA))
    OUT <- cbind(ALPHA,BETA)
  }else if(dist=="Binomial"){  #Binary
    ALPHA = hyp[1] + sum
    BETA = hyp[2] + count - sum
    colnames(ALPHA) = paste0("A",1:ncol(ALPHA))
    colnames(BETA) = paste0("B",1:ncol(ALPHA))
    OUT <- cbind(ALPHA,BETA)
  }else if(dist=="Normal"){ #Nornal(meanvar)   ###????
    M = (sum + hyp[1]/hyp[2])/(count + 1/hyp[2])
    C = 1/(1/hyp[2] + count)
    A = hyp[3] + count/2
    mean = sum/count
    B = hyp[4] + 0.5*(sq - 2*mean*sum + count*mean^2) + 0.5*count*(1/hyp[2])*(mean-hyp[1])^2/(count+1/hyp[2])

    colnames(M) = paste0("M",1:ncol(M))
    colnames(C) = paste0("C",1:ncol(C))
    colnames(A) = paste0("A",1:ncol(A))
    colnames(B) = paste0("B",1:ncol(B))
    OUT <- cbind(M,C,A,B)
  }else if(dist=="Mean"){  ##hyp=(tau,m,c)   
    M = (sum/hyp[1] + hyp[2]/hyp[3])/(count/hyp[1] + 1/hyp[3])
    C = 1/(1/hyp[3] + count/hyp[1])

    colnames(M) <- paste0("M",1:ncol(M))
    colnames(C) <- paste0("C",1:ncol(C))
    OUT <- cbind(M,C)
  }else if(dist=="Var"){ ##hyp=(mu,a,b)
    A = hyp[2] + count/2
    B = hyp[3] + 0.5*(sq - 2*hyp[1]*sum + count*hyp[1]^2)
    colnames(A) <- paste0("A",1:ncol(A))
    colnames(B) <- paste0("B",1:ncol(B))
    OUT <- cbind(A,B)
  }else{
    stop("Invalid distribution assumption")
  }
  
  return(OUT)
}

