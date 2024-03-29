% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UCS3functions.R
\name{BIC.UComp}
\alias{BIC.UComp}
\title{BIC.UComp}
\usage{
\method{BIC}{UComp}(object, ...)
}
\arguments{
\item{object}{Object of class \dQuote{UComp}.}

\item{...}{Additional inputs to function.}
}
\description{
Extract BIC (or SBC) value of UComp object
}
\details{
Selection criteria for models with different number of 
parameters, the smaller BIC the better. The formula used here is
\eqn{BIC=(-2 ln(L) + k ln(n)) / n}, where \eqn{ln(L)} is the log-likelihood
at the optimum, \eqn{k} is the number of parameters plus
non-stationary states and \eqn{n} is the number of observations.
Mind that this formulation differs from the usual definition that
does not divide by \eqn{n}. This makes that BIC(m) and BIC(logLik(m))
give different results, being m an UComp object.
}
\examples{
\dontrun{
y <- log(AirPassengers)
m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
BIC(m1)
}
}
\seealso{
\code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
         \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
}
\author{
Diego J. Pedregal
}
