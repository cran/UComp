% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UCS3functions.R
\name{getp0}
\alias{getp0}
\title{getp0}
\usage{
getp0(y, model = "llt/equal/arma(0,0)", periods = NA)
}
\arguments{
\item{y}{a time series to forecast.}

\item{model}{any valid \code{UComp} model without any ?.}

\item{periods}{vector of fundamental period and harmonics required.}
}
\value{
A set of parameters p0 of an object of class \code{UComp}
to use as input to \code{\link{UC}}, \code{\link{UCforecast}} or \code{\link{UCsetup}}.
}
\description{
Get initial conditions for parameters of \code{UComp} object
}
\details{
Provides initial parameters of a given model for the time series.
They may be changed arbitrarily by the user to include as an input \code{p0} to
\code{UC} or \code{UCforecast} functions (see example below).
There is no guarantee that the model will converge and selecting initial conditions
should be used with care.
}
\examples{
\dontrun{
p0 <- getp0(log(AirPassengers), model = "llt/equal/arma(0,0)")
p0[1] <- 0  # p0[1] <- NA
m <- UCforecast(log(AirPassengers), model = "llt/equal/arma(0,0)", p0 = p0)
}
}
\seealso{
\code{\link{UC}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
         \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
         \code{\link{UChp}}
}
\author{
Diego J. Pedregal
}
