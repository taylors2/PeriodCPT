\name{CircPlot}
\alias{CircPlot}
\alias{Circ_plot}
\alias{Circ_rline}
\alias{Circ_arc}
\alias{Circ_line}
\title{
  Circular plotting functions
}
\description{
  The command Circ_plot() creates a radial plot whilst the 
remaining functions adds different types of lines to the 
created image. Circ_rline() adds a radial segment at a given 
time point, Circ_arc() creates an arc between two time points 
at the same radial distance, and Circ_line() adds a line that 
connect a provided sequence of points.
}

\usage{
Circ_plot(data, time=NULL, rlim=NULL, tmax=NULL, dot=TRUE,
  clocktime=TRUE, main=NULL, tlab="Time", rlab=NULL, cex.lab=1,
  cex.axis=1, cex.main=1, ...)
Circ_rline(t, r1=NULL, r2=NULL, plotinfo, ...)
Circ_arc(r, t1 = NULL, t2 = NULL, plotinfo, ...)
Circ_line(radius, time, close=TRUE, plotinfo, ...)
}

\arguments{
  \item{data}{Radial values at which to draw the points vector.}
  \item{time}{Time vector in terms of integer modulo of periodiod 
length. If NULL then data is taken as a single complete period.}
  \item{rlim}{Bound to the radial axis. If NULL, derived 
from data.}
  \item{tmax}{Upper bound to the periodic time axis. If NULL, 
derived from time.}
  \item{dot}{Logical, should a dot be drawn at the centre to 
offset the zero position of the radial axis. Default, TRUE.}
  \item{main}{Plot title}
  \item{tlab}{Label for the time axis}
  \item{rlab}{Label for the radial axis}
  \item{t}{Time bin defining the angle to which to draw the 
radial line at.}
  \item{r1, r2}{Radial endpoints for the line. If NULL then 
taken as r1=0 and r2=rmax.}
  \item{plotinfo}{plotting arguments returned by Circ_plot()}
  \item{r}{Radial level at which to draw an arc.}
  \item{t1, t2}{Time start and when drawing the arc. If either 
are NULL, then a circle is drawn.}
  \item{radius}{Radial values at which to draw connecting lines.}
  \item{close}{Logical, if TRUE then the end point should be 
connected to the start point.}
  \item{cex.lab, cex.axis, cex.main}{Size of axis labels and title.}
  \item{clocktime}{Logical flag indicating if the values on the 
radial axis should be printed in terms of clock time, TRUE by default.
Otherwise, axis values are printed as fraction of period.}
  \item{...}{Additional plotting arguments.}
}

\details{
  time, t, t1 and t2 can be provided in "hh:mm" format.
}

\value{
  Circ_plot() returns a list of plotting arguments that should 
be provided to any further plotting function.
}

\examples{
  set.seed(10)
  rate = c(rep(1,12),rep(10,12),rep(1,12),rep(10,12))
  info <- Circ_plot(data=rpois(48,rate), time=c(1:24,1:24),
    main="plot", rlab="poisson", rlim=c(0,15))
  
  #Add radial bars
  Circ_rline(t=12, r1=8, r2=10, plotinfo=info, col=4, lwd=2)
  Circ_rline(t=24, r1=8, r2=10, plotinfo=info, col=4, lwd=2)
  
  #Add arc
  Circ_arc(r=9, t1=12, t2=24, plotinfo=info, col=4, lwd=2)
  
  #Add line
  Circ_line(radius=c(rep(1,12),rep(10,12)), time=1:24,
    close=TRUE, plotinfo=info, col=2, lwd=2)
}

\keyword{CircPlot}