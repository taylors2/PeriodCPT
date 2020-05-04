param.prior.make.norm <- function(){
  ##y ~ Norm(theta[1], theta[2]),
  ##theta[1] ~ Norm(param.m, param.c * theta[2])
  ##theta[2] ~ IGamma(param.a, param.b)

  if(missing(param.m)){
    param.m <- 0
  }else{
    if(!is.numeric(param.m) || length(param.m) != 1)
      stop("Normal-InvGamma `param.m` hyper parameter specified incorrectly.")
  }

  if(missing(param.c)){
    param.c <- 1
  }else{
    if(!is.numeric(param.c) || length(param.c) != 1 || any(param.c <= 0))
      stop("Normal-InvGamma `param.c` hyper parameter specified incorrectly.")
  }

  if(missing(param.a)){
    param.a <- 1
  }else{
    if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0))
      stop("Normal-InvGamma `param.a` hyper parameter specified incorrectly.")
  }

  if(missing(param.b)){
    param.b <- 1
  }else{
    if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0))
      stop("Normal-InvGamma `param.b` hyper parameter specified incorrectly.")
  }

  out <- c("param.m" = param.m, "param.c" = param.c,
           "param.a" = param.a, "param.b" = param.b)
  return(out)
}
