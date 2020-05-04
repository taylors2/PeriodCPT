param.prior.make.bern <- function(param.a, param.b){
  ##y~bern(theta), theta~beta(param.a, param.b)
  if(missing(param.a)){
    param.a <- 1
  }else{
    if(!is.numeric(param.a) || length(param.a) != 1 || any(param.a <= 0))
      stop("Beta `param.a` hyper parameter specified incorrectly.")
  }

  if(missing(param.b)){
    param.b <- 1
  }else{
    if(!is.numeric(param.b) || length(param.b) != 1 || any(param.b <= 0))
      stop("Beta `param.b` hyper parameter specified incorrectly.")
  }
  out <- c("param.a" = param.a, "param.b" = param.b)
  return(out)
}
