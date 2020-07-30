plot.master <- function(
  object, fn = NULL, main = NULL, rlim = NULL, type = "p",
  draw.est = TRUE, probs.est = 0.5, param.est = 1, draw.bar = FALSE,
  at.axis = NULL, labels.axis = NULL, rotate.axis = FALSE, ...){

  if(!is.logical(draw.est) | length(draw.est) != 1 | anyNA(draw.est))
    stop("Invalid `draw.est` value.")
  if(!is.logical(draw.bar) | length(draw.bar) != 1 | anyNA(draw.bar))
    stop("Invalid `draw.bar` value.")
  if(is.null(main)) main <- paste0("Period plot of ", deparse(substitute(object)))
  if(!is.character(type) | length(type) != 1) stop("Invalid `type` value.")
  if(max(nchar(type)) != 1) stop("Invalid `type` value.")
  if(!(type %in% c("p","l","o","n")))
    stop("Unknown `type`, must be one of \"p\", \"l\", \"o\" or \"n\".")

  #Format axis information
  if(is.null(at.axis) | length(at.axis) == 0) at.axis <- 1:periodlength(object)
  if(!is.numeric(at.axis)) stop("Inalid 'at.axis' value.")
  L <- at.axis %in% 1:periodlength(object)
  if(!all(L)) stop("Invalid 'at.axis' value, can only be integers from 1 to period length.")
  if(is.null(labels.axis) | length(labels.axis) == 0) labels.axis <- at.axis
  if(is.numeric(labels.axis)) labels.axis <- as.character(labels.axis)
  if(!is.character(labels.axis)) stop("Invalid 'tlables', can only be character or numeric.")
  labels.axis <- labels.axis[rep_len(1:length(labels.axis), length.out = length(at.axis))]
  if(!is.logical(rotate.axis) | is.na(rotate.axis) | length(rotate.axis) == 0) rotate.axis <- FALSE
  if(length(rotate.axis) > 1) rotate.axis <- rotate.axis[1]

  if(!draw.est & !draw.bar & type == "n")
    stop("Nothing to plot with `draw.est = FALSE`, `draw.bar = FALSE` and `type = \"n\"`.")

  ## Grab plotting values (times, data, quantiles)
  ##  -- times and data
  draw.data <- type != "n"
  time_data <- plot.data.format(object, draw.data, fn, ...)

  ##  -- requested quantiles from fit and pcpt probabilities
  if(length(object@date) == 0 & (draw.est | draw.bar)){
    draw.est <- FALSE
    draw.bar <- FALSE
    warning("There are no estimates to add to the plot.")
  }
  Q <- plot.est.format(object, draw.est, probs.est, param.est)
  P <- plot.bar.format(object, draw.bar)

  ##Format the radius limits
  rtick <- plot.evaluate.rtick(rlim, time_data, Q)
  rlim <- range(rtick)

  ##Transform time/value to plotting angle/radius
  RAD_LIM <- c(0.5, 1.5)
  frame_half_span <- RAD_LIM[2] + 0.1*diff(RAD_LIM)

  time_data <- plot.transform.data(time_data, periodlength(object), rlim, RAD_LIM)
  Q         <- plot.transform.data(Q,         periodlength(object), rlim, RAD_LIM)
  P         <- plot.transform.data(P,         periodlength(object), c(0,1), RAD_LIM)
  rtick_std <- plot.transform.data(rtick,     periodlength(object), rlim, RAD_LIM)

  #Create a blank plotting space
  plot.blank.frame(main, frame_half_span, ...)

  #add bars
  if(draw.bar) plot.draw.bars(P, periodlength(object), RAD_LIM, ...)

  #add grid
  plot.grid(N = periodlength(object), rtick = rtick_std, ...)


  #add estimates
  if(draw.est) plot.draw.estimates(Q, ...)

  #add data
  plot.draw.data(time_data, type, join = !is.null(fn), ...)

  #add axes
  plot.draw.axis(periodlength(object), rtick, rtick_std, RAD_LIM, at.axis,
      labels.axis, rotate.axis, draw.est | type != "n", draw.bar, ...)

}


plot.draw.axis <- function(N, rtick, rtick_std, RAD_LIM, at, labels, rotate,
                           draw.EstOrData, draw.bar,...){
  dots <- list(...)
  dots.axis <- dots[names(dots) %in% names(par())]
  dots.axis["cex"] <- ifelse("cex.axis" %in% names(dots.axis), dots.axis["cex.axis"], par("cex.axis"))
  dots.axis["col"] <- ifelse("col.axis" %in% names(dots.axis), dots.axis["col.axis"], par("col.axis"))
  if(!("font" %in% names(dots.axis))){
    dots.axis["font"] <- ifelse("font.axis" %in% names(dots.axis), dots.axis["font.axis"], par("font.axis"))
  }
  if(!draw.EstOrData & !draw.bar) draw.EstOrData <- TRUE

  #Add top axis
  dots.axis$x = rep(0,length(rtick))
  dots.axis$y = rtick_std
  if(draw.EstOrData){
    dots.axis$labels = rtick
  }else{
    dots.axis$labels = paste0(0:(length(rtick)-1),"/",(length(rtick)-1))
  }
  dots.axis$pos = 2
  do.call(text, dots.axis)

  #Add bottom axis
  dots.axis$y = -dots.axis$y
  dots.axis$pos = 4
  if(draw.bar){
    dots.axis$labels = paste0(0:(length(rtick)-1),"/",(length(rtick)-1))
  }else{
    dots.axis$labels = rtick
  }
  do.call(text, dots.axis)

  #add rotation axis
  ax_rad <- RAD_LIM[2] + 0.15*diff(RAD_LIM)
  dots.axis <- dots.axis[names(dots.axis) != "pos"]
  dots.axis$adj = c(0.5, 0.5)
  for(i in 1:length(at)){
    dots.axis$x = ax_rad * cos(pi/2 - (at[i]/N)*2*pi)
    dots.axis$y = ax_rad * sin(pi/2 - (at[i]/N)*2*pi)
    dots.axis$labels = labels[i]
    if(rotate) dots.axis$srt = 360 * (1 - (at[i]/N))
    do.call(text, dots.axis)
  }
}




plot.data.format <- function(object, draw.data = TRUE, fn = NULL, ...){
  N <- periodlength(object)
  if(length(N) == 0) stop("Period length is not defined within pcpt-class object.")
  if(!draw.data){
    return(list(time = 1:N, data = rep(NA, N)))
  }

  datats <- data.set(object)
  if(length(datats) == 0) stop("No data to plot from pcpt-class object.")
  offset <- start(datats)[2]    ##Evaluate time axis
  time <- seq(from = offset, by = 1, length.out = length(datats))
  time <- ((time - 1) %% N) + 1
  data <- as.numeric(datats)
  if(is.function(fn)){
    pts2 <- rep(0, N)
    for(t in 1:N){
      timept_summary <- evaluate_fn(fn = fn, x = data[time == t], ...)
      if(!(is.numeric(timept_summary) | anyNA(timept_summary)) | length(timept_summary) != 1){
        stop("Function `fn` does not return a single numerical or NA value.")
      }
      pts2[t] <- timept_summary
    }
    data <- pts2
    time <- 1:N
  }
  return(list(time = time, data = data))
}


plot.est.format <- function(object, draw.est = TRUE, probs = 0.5, param = 0){
  if(draw.est){
    data <- eval(parse(text = paste0("plot.est.format.",distribution(object),
       "(object = object, probs = probs, param = param)")))
    time <- 1:periodlength(object)
  }else{
    data <- time <- NA
  }
  return(list(time = time, data = data))
}

evaluate_fn <- function(fn, x, ...){
  if(!is.function(fn)) return(NULL)
  fn_string <- paste0(deparse(args(fn)), collapse = "")
  if(grepl("\\.\\.\\.",fn_string)) return(fn(x, ...))

  fn_argstr <- substr(fn_string, 11, nchar(fn_string) - 6)
  commas <- NULL
  count_bracksts <- 0
  for(i in 1:nchar(fn_argstr)){
    if(substr(fn_argstr,i,i) == "," & count_bracksts == 0){
      commas <- c(commas, i)
    }else if(substr(fn_argstr,i,i) == "("){
      count_bracksts <- count_bracksts + 1
    }else if(substr(fn_argstr,i,i) == ")"){
      count_bracksts <- count_bracksts - 1
    }
  }
  fn_argstr_each <- substring(fn_argstr, c(1,commas+1), c(commas-1,nchar(fn_argstr)))
  for(i in 1:length(fn_argstr_each)){
    fn_argstr_each[i] <- gsub(" ","",strsplit(fn_argstr_each[i], "=")[[1]][1])
  }
  dots <- list(...)
  inputs <- dots
  inputs <- inputs[names(inputs) %in% fn_argstr_each]
  L <- fn_argstr_each %in% names(inputs)
  if(all(L)){
    inputs[[(fn_argstr_each[!L])[1]]] <- x
  }else{
    inputs[[fn_argstr_each[1]]] <- x
  }
  return(do.call(fn, inputs))
}

plot.evaluate.rtick <- function(rlim = NULL, time_data = NA, Q = NA){
  if(is.null(rlim)){
    tmp <- c(time_data$data, Q$data)
    tmp <- tmp[is.finite(tmp)]
    if(length(tmp) == 0){
      rlim <- c(0, 1)
    }else{
      rlim <- range(tmp)
    }
  }

  if(!is.numeric(rlim) | length(rlim) != 2) stop("Invalid `rlim` value.")
  if(any(is.infinite(rlim))) stop("Need finite `rlim` values.")
  rlim <- sort(rlim)
  if(rlim[1] == rlim[2]) rlim <- rlim + c(-0.8,+0.8)
  rtick <- pretty(rlim)
  return(rtick)
}

plot.bar.format <- function(object, draw.bar = TRUE){
  if(draw.bar){
    data <- table_pcpt(object) / n.iter(object)
    data <- data[grepl("tau",names(data))]
    time <- as.numeric(sub("tau","",names(data)))
    data <- unname(data)
    data <- data[order(time)]
    time <- time[order(time)]
  }else{
    time <- 1:periodlength(object)
    data <- rep(NA, length(time))
  }
  return(list(time = time, data = data))
}

plot.transform.data <- function(info, N, rlim, RAD_LIM){
  if(!is.list(info)){
    info <- pmax(0, RAD_LIM[1] + (info - rlim[1])*diff(RAD_LIM)/diff(rlim))
    WARN_RANGE <- any(info < 0.5*RAD_LIM[1] |
                        info > (RAD_LIM[2] + 0.1*diff(RAD_LIM)), na.rm = TRUE)
  }else{
    info$time <- pi / 2 - (2 * pi * info$time / N)
    if(is.matrix(info$data)){
      for(i in 1:ncol(info$data)){
        info$data[,i] <- pmax(0, RAD_LIM[1] + (info$data[,i] - rlim[1])*diff(RAD_LIM)/diff(rlim))
      }
    }else{
      info$data <- pmax(0, RAD_LIM[1] + (info$data - rlim[1])*diff(RAD_LIM)/diff(rlim))
    }
    WARN_RANGE <- any(info$data < 0.5*RAD_LIM[1] |
                        info$data > (RAD_LIM[2] + 0.1*diff(RAD_LIM)), na.rm = TRUE)
  }
  if(WARN_RANGE)
    warning("Plotted values are outside the `rlim` radius range.")
  return(info)
}

plot.blank.frame <- function(main, frame_half_span, ...){
  plot.dots <- list(...)
  plot.dots <- plot.dots[ names(plot.dots) %in% names(par()) ]
  plot.dots$xlab  = ""
  plot.dots$ylab  = ""
  plot.dots$frame = FALSE
  plot.dots$type  = "n"
  plot.dots$xaxt  = "n"
  plot.dots$yaxt  = "n"
  plot.dots$asp   = TRUE
  plot.dots$xlim  = c(-1, 1) * frame_half_span
  plot.dots$ylim  = c(-1, 1) * frame_half_span
  plot.dots$main  = main
  plot.dots$x  = 0
  plot.dots$y  = 0
  do.call(plot, args = plot.dots)
}


plot.grid <- function(N = 1, rtick, ...){
  dots <- list(...)
  col <- "grey60"
  if("col.grid" %in% names(dots)) col <- dots$col.grid[1]
  lwd <- par("lwd")
  if("lwd.grid" %in% names(dots)) lwd <- dots$lwd.grid[1]
  lty <- "dotted"
  if("lty.grid" %in% names(dots)) lty <- dots$lty.grid[1]
  theta <- seq(pi/2,-3*pi/2,len = 101)
  for(r in rtick){
    lines(r*cos(theta), r*sin(theta), col = col, lty = lty, lwd = lwd)
  }
  rtheta <- pi/2 - (1:N)*2*pi/N
  segments(min(rtick) * cos(rtheta), min(rtick) * sin(rtheta),
           max(rtick) * cos(rtheta), max(rtick) * sin(rtheta),
           col = col, lty = lty, lwd = lwd)
  segments(0, 0, min(rtick) * cos(pi * c(1,0.5,0,-0.5)),
           min(rtick) * sin(pi * c(1,0.5,0,-0.5)),
           col = col, lty = lty, lwd = lwd)
}

plot.draw.bars <- function(P, N, RAD_LIM, ...){
  dots <- list(...)
  bar.dots <- list(col = "grey80", lty = "solid",
                   lwd = par("lwd"), border = "black")
  for(case in grep(".bar", names(dots), value = TRUE)){
    bar.dots[[sub(".bar","",case)]] <- dots[[case]]
  }
  theta <- seq(from = -pi/N, to = pi/N, length.out = 101)
  for(i in 1:N){
    r <- c(rep(P$data[i],length(theta)), rep(RAD_LIM[1], length(theta)))
    bar.dots$x = r*cos(P$time[i] + c(theta,rev(theta)))
    bar.dots$y = r*sin(P$time[i] + c(theta,rev(theta)))
    do.call(polygon, bar.dots)
  }
}

plot.draw.estimates <- function(Q, ...){
  dots.est <- list(...)
  col <- "red"
  if("col.est" %in% names(dots.est)) col <- dots.est$col.est
  col <- rep_len(col, ncol(Q$data))
  lwd <- par("lwd")
  if("lwd.est" %in% names(dots.est)) lwd <- dots.est$lwd.est
  lwd <- rep_len(lwd, ncol(Q$data))
  lty <- "solid"
  if("lty.est" %in% names(dots.est)) lty <- dots.est$lty.est
  lty <- rep_len(lty, ncol(Q$data))
  for(i in 1:ncol(Q$data)){
    ang <- c(Q$time,Q$time[1])
    rad <- unname(c(Q$data[,i],Q$data[1,i]))
    lines(x = rad * cos(ang), y = rad  *sin(ang),
          col = col[i], lty = lty[i], lwd = lwd[i])
  }
}

plot.draw.data <- function(time_data, type, join, ...){
  points.dots <- list(...)
  points.dots <- points.dots[names(points.dots) %in% names(par())]
  points.dots$x = time_data$data*cos(time_data$time)
  points.dots$y = time_data$data*sin(time_data$time)
  points.dots$type = type
  do.call(points, points.dots)
  if(join & (type == "l" | type == "o")){
    ends <- c(length(time_data$data), 1)
    points.dots$type <- "l"
    points.dots$x <- points.dots$x[ends]
    points.dots$y <- points.dots$y[ends]
    do.call(lines, points.dots)
  }
}

