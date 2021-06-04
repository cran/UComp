## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(UComp)

## ----load data----------------------------------------------------------------
load(url("https://github.com/djpedregal/djpedregal.github.io/blob/master/ipi.Rdata?raw=true"))
ipi = ts(y[1 : 208], start = 2000, frequency = 12)
plot(ipi)

## ----run model----------------------------------------------------------------
m = UC(log(ipi))

## ----plot components----------------------------------------------------------
plot(m)

## ----run model 2--------------------------------------------------------------
m = UC(log(ipi), model = "rw/equal/arma(0,0)", h = 24, verbose = FALSE)

## ----outliers-----------------------------------------------------------------
m = UC(log(ipi), u = u, outlier = 4)

## ----plot outliers------------------------------------------------------------
plot(m)

## ----forecasts----------------------------------------------------------------
plot(m$yFor, ylim = c(4.3, 4.85))
lines(m$yFor + 2 * sqrt(m$yForV), col = "red")
lines(m$yFor - 2 * sqrt(m$yForV), col = "red")

