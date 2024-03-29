% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{roots}
\alias{roots}
\title{roots}
\usage{
roots(x)
}
\arguments{
\item{x}{coefficients of polynomial in descending order}
}
\value{
Roots of polynomial
}
\description{
Roots of polynomial
}
\examples{
roots(c(1, -2 ,1))
roots(conv(c(1, -1), c(1, 0.8)))
}
\seealso{
\code{\link{colMedians}}, \code{\link{rowMedians}}, \code{\link{tests}},
         \code{\link{sumStats}}, \code{\link{gaussTest}}, \code{\link{ident}},
         \code{\link{cusum}}, \code{\link{varTest}}, \code{\link{conv}},
         \code{\link{armaFilter}}, \code{\link{dif}},
         \code{\link{zplane}}, \code{\link{acft}}, \code{\link{slide}},
         \code{\link{plotSlide}}, \code{\link{Accuracy}}, \code{\link{tsDisplay}},
         \code{\link{size}}
}
\author{
Diego J. Pedregal
}
