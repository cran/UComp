% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ARIMAestim.R
\name{ARIMAestim}
\alias{ARIMAestim}
\title{ARIMAestim}
\usage{
ARIMAestim(m)
}
\arguments{
\item{m}{an object of type \code{ARIMA} created with \code{ARIMAmodel}}
}
\value{
The same input object with the appropriate fields 
filled in, in particular:
\item{p}{Estimated parameters}
\item{yFor}{Forecasted values of output}
\item{yForV}{Variance of forecasted values of output}
\item{ySimul}{Bootstrap simulations for forecasting distribution evaluation}
}
\description{
Estimates and forecasts ARIMA models
}
\details{
\code{ARIMAestim} estimates and forecasts a time series using 
an ARIMA model
}
\examples{
\dontrun{
m1 <- ARIMAsetup(log(gdp))
m1 <- ARIMAestim(m1)
}
}
\seealso{
\code{\link{ARIMA}}, \code{\link{ARIMAmodel}}, \code{\link{ARIMAvalidate}},
}
\author{
Diego J. Pedregal
}
