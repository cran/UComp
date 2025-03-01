% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ETSS3functions.R, R/TETSS3functions.R
\name{plot.ETS}
\alias{plot.ETS}
\alias{plot.TETS}
\title{plot.ETS}
\usage{
\method{plot}{ETS}(x, ...)

\method{plot}{TETS}(x, ...)
}
\arguments{
\item{x}{Object of class \dQuote{TETS}.}

\item{...}{Additional inputs to function.}
}
\description{
Plot components of ETS object

Plot components of TETS object
}
\details{
See help of \code{ETS}.

See help of \code{TETS}.
}
\examples{
\dontrun{
m1 <- ETSforecast(log(gdp))
plot(m1)
}
\dontrun{
m1 <- TETSforecast(log(gdp))
plot(m1)
}
}
\seealso{
\code{\link{ETS}}, \code{\link{ETSforecast}}, \code{\link{ETSvalidate}},
         \code{\link{ETScomponents}}, \code{\link{ETSestim}}

\code{\link{TETS}}, \code{\link{TETSforecast}}, \code{\link{TETSvalidate}},
         \code{\link{TETScomponents}}, \code{\link{TETSestim}}
}
\author{
Diego J. Pedregal
}
