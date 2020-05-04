param.prior.make.mean <- function(param.m, param.c, var){
  ##y ~ Norm(theta, var), theta~Norm(param.m, param.c)

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

  if(missing(var)){
    var <- 1
  }else{
    if(!is.numeric(var) || length(var) != 1 || any(var <= 0))
      stop("Known `var` parameter in distribution function is specified incorrectly.")
  }

  out <- c("var" = var, "param.m" = param.m, "param.c" = param.c)
  return(out)

}
