% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UCS3functions.R
\name{AIC.UComp}
\alias{AIC.UComp}
\title{AIC.UComp}
\usage{
\method{AIC}{UComp}(object, ..., k = 2)
}
\arguments{
\item{object}{Object of class \dQuote{UComp}.}

\item{...}{Additional inputs to function.}

\item{k}{The penalty per parameter to be used.}
}
\description{
Extract AIC value of UComp object
}
\details{
Selection criteria for models with different number of 
parameters, the smaller AIC the better. The formula used here is
\eqn{AIC=-2 (ln(L) - k) / n}, where \eqn{ln(L)} is the log-likelihood
at the optimum, \eqn{k} is the number of parameters plus
non-stationary states and \eqn{n} is the number of observations.
Mind that this formulation differs from the usual definition that
does not divide by \eqn{n}. This makes that AIC(m) and AIC(logLik(m))
give different results, being m an UComp object.
}
\examples{
\dontrun{
y <- log(AirPassengers)
m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
AIC(m1)
}
}
\seealso{
\code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
         \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
}
\author{
Diego J. Pedregal
}
