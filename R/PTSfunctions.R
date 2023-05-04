#' @title PTSsetup
#' @description Run up PTS general univariate MSOE models
#'
#' @details See help of \code{PTS}.
#'
#' @param y a time series to forecast (it may be either a numerical vector or
#' a time series object). This is the only input required. If a vector, the additional
#' input \code{s} should be supplied compulsorily (see below).
#' @param u a matrix of input time series. If
#' the output wanted to be forecast, matrix \code{u} should contain future values for inputs.
#' @param model the model to estimate. It is a single string indicating the type of
#' model for each component with one or two letters:
#' \itemize{
#' \item Error: ? / N / A
#'
#' \item Trend: ? / N / A / Ad / L
#'
#' \item Seasonal: ? / N / A / D (trigonometric with different variances)
#'
#' }
#' @param s seasonal period of time series (1 for annual, 4 for quarterly, ...)
#' @param h forecast horizon. If the model includes inputs h is not used, the lenght of u is used instead.
#' @param criterion information criterion for identification ("aic", "bic" or "aicc").
#' @param lambda Box-Cox lambda parameter (NULL: estimate)
#' @param armaIdent check for arma models for error component (TRUE / FALSE).
#' @param verbose intermediate estimation output (TRUE / FALSE)
#'
#' @author Diego J. Pedregal
#'
#' @return An object of class \code{PTS}. It is a list with fields including all the inputs and
#'         the fields listed below as outputs. All the functions in this package fill in
#'         part of the fields of any \code{PTS} object as specified in what follows (function
#'         \code{PTS} fills in all of them at once):
#'
#' After running \code{PTSmodel} or \code{PTSestim}:
#' \itemize{
#' \item p0:       Initial values for parameter search
#' \item p:        Estimated parameters
#' \item lambda:   Estimated Box-Cox lambda parameter
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecasted values of output
#' \item yForV:    Variance of forecasted values of output
#' }
#'
#' After running \code{PTSvalidate}:
#' \itemize{
#' \item table: Estimation and validation table
#' }
#'
#' After running \code{PTScomponents}:
#' \itemize{
#' \item comp:  Estimated components in matrix form

#' }
#'
#' Standard methods applicable to PTS objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @seealso \code{\link{PTS}}, \code{\link{PTSmodel}}, \code{\link{PTSvalidate}},
#'          \code{\link{PTScomponents}}, \code{\link{PTSestim}}
#'
#' @examples
#' \dontrun{
#' m1 <- PTSsetup(log(AirPassengers))
#' }
#' @rdname PTSsetup
#' @export
PTSsetup <- function(y, u = NULL, model="???", s = frequency(y), h = 2 * s, 
                criterion = "aic", lambda = 1, armaIdent = FALSE, 
                verbose = FALSE){
        if (s <= 1)
                s = 0
        # modelU = PTS2modelUC(model)
        out = list(
                y = y,
                u = u,
                model = model,
                s = s,
                h = h,
                p0 = NA,
                criterion = criterion,
                lambda = lambda,
                verbose = verbose,
                armaIdent = armaIdent,
                armaOrders = c(0,0),
                yFor = NA,
                yForV = NA,
                comp = NA,
                table = "",
                p = NA,
                v = NA,
                modelUC = NA
        )
        return(structure(out, class = "PTS"))
}
#' @title PTSmodel
#' @description Estimates and forecasts PTS general univariate models
#'
#' @details \code{PTSmodel} is a function for modelling and forecasting univariate
#' time series according to Power-Trend-Seasonal (PTS).
#' It sets up the model with a number of control variables that
#' govern the way the rest of functions in the package work. It also estimates
#' the model parameters by Maximum Likelihood and forecasts the data.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @inheritParams PTSsetup
#'
#' @return An object of class \code{PTS}. It is a list with fields including all the inputs and
#'         the fields listed below as outputs. All the functions in this package fill in
#'         part of the fields of any \code{PTS} object as specified in what follows (function
#'         \code{PTS} fills in all of them at once):
#'
#' After running \code{PTSmodel} or \code{PTSestim}:
#' \itemize{
#' \item p0:       Initial values for parameter search
#' \item p:        Estimated parameters
#' \item lambda:   Estimated Box-Cox lambda parameter
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecasted values of output
#' \item yForV:    Variance of forecasted values of output
#' }
#'
#' After running \code{PTSvalidate}:
#' \itemize{
#' \item table: Estimation and validation table
#' }
#'
#' After running \code{PTScomponents}:
#' \itemize{
#' \item comp:  Estimated components in matrix form
#' }
#'
#' @author Diego J. Pedregal
#'
#' @seealso \code{\link{PTS}}, \code{\link{PTSsetup}}, \code{\link{PTSvalidate}},
#'          \code{\link{PTScomponents}}, \code{\link{PTSestim}}
#'
#' @examples
#' \dontrun{
#' m1 <- PTSmodel(log(AirPassengers))
#' }
#' @rdname PTSmodel
#' @export
PTSmodel <- function(y, u = NULL, model="???", s = frequency(y), h = 2 * s, 
                     criterion = "aicc", lambda = 1, armaIdent = FALSE,
                     verbose = FALSE){
        if (s <= 1)
                s = 0
        m = PTSsetup(y, u, model, s, h, criterion, lambda, armaIdent, verbose)
        m = PTSestim(m)
        return(m)
}
#' @title PTS
#' @description Estimates, forecasts and smooth PTS general univariate models
#'
#' @details \code{PTS} is a function for modelling and forecasting univariate
#' time series according to Power-Trend-Seasonal (PTS).
#' It sets up the model with a number of control variables that
#' govern the way the rest of functions in the package work. It also estimates
#' the model parameters by Maximum Likelihood and forecasts the data.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @inheritParams PTSsetup
#'
#' @return An object of class \code{PTS}. It is a list with fields including all the inputs and
#'         the fields listed below as outputs. All the functions in this package fill in
#'         part of the fields of any \code{PTS} object as specified in what follows (function
#'         \code{PTS} fills in all of them at once):
#'
#' After running \code{PTSmodel} or \code{PTSestim}:
#' \itemize{
#' \item p0:       Initial values for parameter search
#' \item p:        Estimated parameters
#' \item lambda:   Estimated Box-Cox lambda parameter
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecasted values of output
#' \item yForV:    Variance of forecasted values of output
#' }
#'
#' After running \code{PTSvalidate}:
#' \itemize{
#' \item table: Estimation and validation table
#' }
#'
#' After running \code{PTScomponents}:
#' \itemize{
#' \item comp:  Estimated components in matrix form
#' }
#'
#' @author Diego J. Pedregal
#'
#' @seealso \code{\link{PTSmodel}}, \code{\link{PTSsetup}}, \code{\link{PTSvalidate}},
#'          \code{\link{PTScomponents}}, \code{\link{PTSestim}}
#'
#' @examples
#' \dontrun{
#' m1 <- PTS(log(AirPassengers))
#' }
#' @rdname PTS
#' @export
PTS <- function(y, u = NULL, model="???", s = frequency(y), h = 2 * s,
                criterion = "aicc", lambda = 1, armaIdent = FALSE, 
                verbose = FALSE){
        if (s <= 1)
                s = 0
        m = PTSsetup(y, u, model, s, h, criterion, lambda, armaIdent, verbose)
        m = PTSestim(m)
        m = PTSvalidate(m, verbose)
        m = PTScomponents(m)
        return(m)
}
#' @title PTSestim
#' @description Estimates and forecasts PTS models
#'
#' @details \code{PTSestim} estimates and forecasts a time series using an
#' a PTS model
#'
#' @param m an object of type \code{PTS} created with \code{PTSmodel}
#'
#' @return The same input object with the appropriate fields
#' filled in, in particular:
#' \itemize{
#' \item p0:       Initial values for parameter search
#' \item p:        Estimated parameters
#' \item lambda:   Estimated Box-Cox lambda parameter
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecasted values of output
#' \item yForV:    Variance of forecasted values of output
#' }
#'
#' @author Diego J. Pedregal
#'
#' @seealso \code{\link{PTSmodel}}, \code{\link{PTSsetup}}, \code{\link{PTSvalidate}},
#'          \code{\link{PTScomponents}}, \code{\link{PTS}}
#'
#' @examples
#' \dontrun{
#' m1 <- PTSsetup(log(AirPassengers))
#' m1 <- PTSestim(m1)
#' }
#' @rdname PTSestim
#' @export
PTSestim <- function(m){
        modelUC = PTS2modelUC(m$model, m$armaOrders)
        if (m$s < 2)
                periods = 1
        else
                periods = m$s / (1 : floor(m$s / 2))
        mUC = UCsetup(m$y, m$u, modelUC, m$h, m$lambda, 0, FALSE, m$criterion,
                      periods, m$verbose, FALSE, -9999.9, m$armaIdent, NULL,
                      "rw/llt/srw", "none/linear/different",
                      "none/arma(0,0)")
        mUC$hidden$MSOE = TRUE
        mUC$hidden$PTSnames = TRUE
        mUC = UCestim(mUC)
        if (mUC$model == "error")
                return(m)
        m$armaOrders = modelUC2arma(mUC$model)
        m$model = modelUC2PTS(mUC$model)
        m$p0 = mUC$p0
        m$lambda = mUC$lambda
        m$yFor = mUC$yFor
        m$yForV = mUC$yForV
        m$p = mUC$p
        m$modelUC = mUC
        return(m)
}
#' @title PTSvalidate
#' @description Shows a table of estimation and diagnostics results for PTS models
#'
#' @param m an object of type \code{PTS} created with \code{PTSmodel}
#' @param verbose verbose mode TRUE/FALSE
#'
#' @return The same input object with the appropriate fields
#' filled in, in particular:
#' \itemize{
#' \item table: Estimation and validation table
#' }
#'
#' @author Diego J. Pedregal
#'
#' @seealso \code{\link{PTSmodel}}, \code{\link{PTSsetup}}, \code{\link{PTSestim}},
#'          \code{\link{PTScomponents}}, \code{\link{PTS}}
#'
#' @examples
#' \dontrun{
#' m1 <- PTSmodel(log(AirPassengers))
#' m1 <- PTSvalidate(m1)
#' }
#' @rdname PTSvalidate
#' @export
PTSvalidate <- function(m, verbose = TRUE){
        m$modelUC = UCvalidate(m$modelUC, FALSE)
        m$table = m$modelUC$table
        # Buscando test heterocedasticidad con valor p NaN
        ind = which(grepl("nan", m$table))
        if (any(ind)){
                for (i in 1 : length(ind)){
                        line = m$table[ind[i]]
                        df = as.numeric(substr(line, 9, 12))
                        Fstat = as.numeric(substr(line, 15, 31))
                        pval = round(pf(Fstat, df, df), 4)
                        line = gsub("   nan", pval, line)
                        m$table[ind[i]] = line
                }
        }
        if (verbose)
                cat(m$table)
        m$v = m$modelUC$v
        return(m)
}
#' @title PTScomponents
#' @description Estimates components of PTS models
#'
#' @param m an object of type \code{PTS} created with \code{PTSmodel}
#'
#' @return The same input object with the appropriate fields
#' filled in, in particular:
#' \itemize{
#' \item comp:  Estimated components in matrix form
#' }
#'
#' @author Diego J. Pedregal
#'
#' @seealso \code{\link{PTSmodel}}, \code{\link{PTSsetup}}, \code{\link{PTSestim}},
#'          \code{\link{PTSvalidate}}, \code{\link{PTS}}
#'
#' @examples
#' \dontrun{
#' m1 <- PTS(log(AirPassengers))
#' m1 <- PTScomponents(m1)
#' }
#' @rdname PTScomponents
#' @export
PTScomponents <- function(m){
        m$modelUC = UCcomponents(m$modelUC)
        names = colnames(m$modelUC$comp)
        nComp = length(names)
        ind = c(1, which(names == "Seasonal"), which(names == "Slope"))
        pos = which(names == "Irregular")
        if (length(pos) > 0 && pos < length(names)){
                ind = c(ind, (pos + 1) : length(names))
        }
        m$comp = cbind(m$modelUC$v, m$modelUC$comp[, 1], m$modelUC$comp[, ind])
        if (ncol(m$comp) == 3)
                m$comp[, 2] = m$comp[, 3]
        else
                m$comp[, 2] = rowSums(m$comp[, 3 : ncol(m$comp)])
        colnames(m$comp) = c("Error", "Fit", names[ind])
        return(m)
}
#' @title modelUC2arma
#' @description Extracts arma part of modelUC model
#'
#' @param model a UC model
#'
#' @return arma orders
#'
#' @author Diego J. Pedregal
#'
#' @rdname modelUC2arma
#' @export
modelUC2arma <- function(model){
        # ARMA model orders
        pos1 = regexpr("arma", model, fixed = TRUE)[1]
        armaModel = substr(model, pos1, nchar(model))
        pos1 = pos1 + 5
        pos2 = regexpr(",", model, fixed = TRUE)[1] - 1
        pos3 = regexpr(")", model, fixed = TRUE)[1] - 1
        ar = as.numeric(substr(model, pos1, pos2))
        ma = as.numeric(substr(model, pos2 + 2, pos3))
        if (is.na(ar)) ar = 0
        if (is.na(ma)) ma = 0
        return(c(ar, ma))
}
#' @title modelUC2PTS
#' @description Translates modelUC to model PTS
#'
#' @param modelUC a UC model
#'
#' @return a PTS model
#'
#' @author Diego J. Pedregal
#'
#' @rdname modelUC2PTS
#' @export
modelUC2PTS <- function(modelUC){
        # removing cycle from UC model
        modelUC = sub("/none/", "/", modelUC, fixed = TRUE)
        # extracting components
        aux = gregexpr("/", modelUC)[[1]]
        trend = substr(modelUC, 1, aux[1] - 1)
        seasonal = substr(modelUC, aux[1] + 1, aux[2] - 1)
        noise = substr(modelUC, aux[2] + 1, nchar(modelUC))
        # noise
        model = "A"
        if (noise == "none")
                model = "N"
        # trend
        if (trend == "rw")
                model = paste0(model, "N")
        else if (trend == "srw")
                model = paste0(model, "Ad")
        else if (trend == "llt")
                model = paste0(model, "A")
        else if (trend == "td")
                model = paste0(model, "L")
        # seasonal
        if (seasonal == "none")
                model = paste0(model, "N")
        else if (seasonal == "equal")
                model = paste0(model, "E")
        else if (seasonal == "different")
                model = paste0(model, "D")
        else if (seasonal == "linear")
                model = paste0(model, "L")
        return(model)
}
#' @title PTS2modelUC
#' @description Translates PTS model to UC model
#'
#' @param model a PTS model
#' @param armaOrders arma orders of noise model
#'
#' @return a UC model
#'
#' @author Diego J. Pedregal
#'
#' @rdname PTS2modelUC
#' @export
 PTS2modelUC <- function(model, armaOrders = c(0, 0)){
        modelU = ""
        n = nchar(model)
        # noise
        aux = tolower(substr(model, 1, 1))
        if (aux == "?")
                modelU = "/?"
        else if (aux == "n")
                modelU = paste0("/none")
        else if (aux == "a")
                modelU = paste0("/arma(", armaOrders[1], ",", armaOrders[2], ")")
        else
                stop("ERROR: incorrect error model!!")
        # seasonal
        aux = tolower(substr(model, n, n))
        if (aux == "?")
                modelU = paste0("/?", modelU)
        else if (aux == "n")
                modelU = paste0("/none", modelU)
        else if (aux == "a")
                modelU = paste0("/linear", modelU)
        else if (aux == "e")
                modelU = paste0("/equal", modelU)
        else if (aux == "d")
                modelU = paste0("/different", modelU)
        else
                stop("ERROR: incorrect seasonal model!!")
        # trend
        aux = tolower(substr(model, 2, n - 1))
        if (aux == "?")
                modelU = paste0("?/none", modelU)
        else if (aux == "n")
                modelU = paste0("rw/none", modelU)
        else if (aux == "a")
                modelU = paste0("llt/none", modelU)
        else if (aux == "ad")
                modelU = paste0("srw/none", modelU)
        else if (aux == "l")
                modelU = paste0("td/none", modelU)
        else
                stop("ERROR: incorrect trend model!!")
        return(modelU)
}


