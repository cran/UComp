% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ETSS3functions.R, R/PTSS3functions.R,
%   R/TETSS3functions.R
\name{fitted.ETS}
\alias{fitted.ETS}
\alias{fitted.PTS}
\alias{fitted.TETS}
\title{fitted.ETS}
\usage{
\method{fitted}{ETS}(object, ...)

\method{fitted}{PTS}(object, ...)

\method{fitted}{TETS}(object, ...)
}
\arguments{
\item{object}{Object of class \dQuote{TETS}.}

\item{...}{Additional inputs to function.}
}
\description{
Fitted output values of ETS object

Fitted output values of PTS object

Fitted output values of TETS object
}
\details{
See help of \code{ETS}.

See help of \code{PTS}.

See help of \code{TETS}.
}
\examples{
\dontrun{
m1 <- ETSmodel(log(gdp))
fitted(m1)
}
\dontrun{
m1 <- PTSmodel(log(AirPassengers))
fitted(m1)
}
\dontrun{
m1 <- TETSmodel(log(gdp))
fitted(m1)
}
}
\seealso{
\code{\link{ETS}}, \code{\link{ETSmodel}}, \code{\link{ETSvalidate}},
         \code{\link{ETScomponents}}, \code{\link{ETSestim}}

\code{\link{PTS}}, \code{\link{PTSmodel}}, \code{\link{PTSvalidate}},
         \code{\link{PTScomponents}}, \code{\link{PTSestim}}

\code{\link{TETS}}, \code{\link{TETSmodel}}, \code{\link{TETSvalidate}},
         \code{\link{TETScomponents}}, \code{\link{TETSestim}}
}
\author{
Diego J. Pedregal
}
