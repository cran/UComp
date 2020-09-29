#' @title print.UComp
#' @description Prints an UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param x Object of class \dQuote{UComp}.
#' @param ... Additional inputs to handle the way to print output.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' print(m1)
#' @rdname print.UComp
#' @export
print.UComp = function(x, ...){
    x = UCvalidate(x, TRUE)
}
#' @title summary.UComp
#' @description Prints an UComp object on screen
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @details See help of \code{UC}.
#'
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' summary(m1)
#' @rdname summary.UComp
#' @export
summary.UComp = function(object, ...){
    print(object)
}
#' @title plot.UComp
#' @description Plot components of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param x Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' plot(m1)
#' @rdname plot.UComp
#' @export
plot.UComp = function(x, ...){
    if (length(x$comp) < 2){
        x = UCcomponents(x)
    }
    if (is.ts(x$comp)){
        plot(x$comp, main = "Time Series Decomposition")
    } else {
        plot(ts(x$comp, frequency = x$periods[1]),
             main = "Time Series Decomposition")
    }
}
#' #' @title autoplot.UComp
#' #' @description Plot components of UComp object
#' #'
#' #' @details See help of \code{UC}.
#' #'
#' #' @author Diego J. Pedregal
#' #' 
#' #' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#' #'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#' #'          
#' #' @examples
#' #' y <- log(AirPassengers)
#' #' m1 <- UCmodel(y)
#' #' autoplot(m1)
#' #' @rdname autoplot.UComp
#' #' @export
#' autoplot.UComp = function(comp){
#'     # Smoothing if necessary
#'     if (any(class(comp) == "UComp")){
#'         if (length(comp$comp) < 2){
#'             comp = UCcomponents(comp)
#'         }
#'         comp = comp$comp
#'     }
#'     autoplot(comp[, colnames(comp)], facets = TRUE) +
#'         xlab("") + ylab("") +
#'         ggtitle("Time Series Components")
#' }
#' @title fitted.UComp
#' @description Fitted output values of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' fitted(m1)
#' @rdname fitted.UComp
#' @export
fitted.UComp = function(object, ...){
    if (length(object$yFit) < 2){
        object = UCsmooth(object)
    }
    return(object$yFit)
}
#' @title residuals.UComp
#' @description Standardised innovations values of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' residuals(m1)
#' @rdname residuals.UComp
#' @export
residuals.UComp = function(object, ...){
    if (length(object$yFit) < 2){
        object = UCfilter(object)
    }
    return(object$v)
}
#' @title logLik.UComp
#' @description Extract log Likelihood value of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' logLik(m1)
#' @rdname logLik.UComp
#' @export
logLik.UComp = function(object, ...){
    out = object$criteria[1]
    class(out) = "logLik"
    attr(out, "df") = length(object$p) - 1
    return(out)
}
#' @title coef.UComp
#' @description Extract model coefficients of UComp object
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' coef(m1)
#' @rdname coef.UComp
#' @export
coef.UComp = function(object, ...){
    if (length(object$table) < 2){
        object = UCvalidate(object, FALSE)
    }
    # Number of u's
    nu = dim(object$u)[1]
    if (dim(object$u)[2] == 2){
        nu = 0
    }
    parameters = matrix(NA, length(object$p) + nu, 1)
    parametersNames = matrix("", length(parameters))
    # Looking for parameters in output table
    hyphen = 1
    rowM = 1
    i = 1
    while (hyphen < 4){
        rowM = rowM + 1
        if (hyphen > 2){
            lineI = object$table[rowM]
            if (substr(object$table[rowM], 1, 1) != "-"){
                # Parameter name
                namePar = substr(lineI, 1, gregexpr(pattern =':', lineI))
                colonN = nchar(namePar) + 1
                namePar = gsub(" ", "", substr(namePar, 1, colonN - 2), fixed = TRUE)
                # Parameter value
                blanks = gregexpr(pattern =' ', substr(lineI, colonN, nchar(lineI)))
                aux = which(diff(blanks[[1]]) > 1)[1] + 1
                aux2 = gsub("*", "", substr(lineI, colonN, colonN + blanks[[1]][aux]), fixed = TRUE)
                parameters[i, 1] = as.numeric(gsub(" ", "", aux2))
                parametersNames[i] = namePar
                i = i + 1
            } 
            # else {
            #     hyphen = hyphen + 1
            # }
        }
        if (substr(object$table[rowM], 1, 1) == "-"){
            hyphen = hyphen + 1
        }
    }
    rownames(parameters) = parametersNames
    return(parameters)
}
#' @title Forecast
#' @description Forecasting using structural Unobseved Components models
#'
#' @details See help of \code{UC}.
#'
#' @param UCompObject UComp object
#' @param level Confidence level for prediction intervals
#' @param fan Set level to seq(51, 99, 3) when fan = TRUE
#' 
#' @author Diego J. Pedregal
#' 
#' @return An object of class \code{forecast}.
#' 
#' Function \code{summary} produces a summary of results, while the function \code{plot} produces a
#' plot of the forecast intervals.
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' \dontrun{
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/eq/arma(0,0)")
#' f1 <- Forecast(m1)
#' plot(f1)
#' }
#' @rdname Forecast
#' @export
Forecast = function(UCompObject,
                    level = c(80, 95),
                    fan = FALSE){
    # Creates forecast object from UComp object
    if (fan){
        level = seq(51, 99, by=3)
    }
    if (length(UCompObject$table) < 2){
        UCompObject = UCvalidate(UCompObject, FALSE)
    }
    if (length(UCompObject$yFit) < 2){
        UCompObject = UCsmooth(UCompObject)
    }
    # Number of u's
    nu = dim(UCompObject$u)[1]
    if (dim(UCompObject$u)[2] == 2){
        nu = 0
    }
    # Number of parameters
    nPar = length(UCompObject$p) - 1 + nu
    h = length(UCompObject$yFor)
    # 2 sided confidence levels
    if (is.null(level)){
        upper = matrix(NA, h, 1)
        lower = upper
    } else {
        level = matrix(level, length(level), 1)
        level2 = level / 100 + (1 - level / 100) / 2
        coefs = qt(level2, df = length(UCompObject$y) - nPar)
        nCoefs = length(coefs)
        upper = matrix(UCompObject$yFor, h, nCoefs) + t(matrix(coefs, nCoefs, h)) * 
                matrix(sqrt(UCompObject$yForV), h, nCoefs)
        lower = matrix(UCompObject$yFor, h, nCoefs) - t(matrix(coefs, nCoefs, h)) * 
                matrix(sqrt(UCompObject$yForV), h, nCoefs)
        colnames(lower) = level
        colnames(upper) = level
    }
    if (is.ts(UCompObject$yFor)){
        freq = frequency(UCompObject$yFor)
        upper = ts(upper, start(UCompObject$yFor, frequency = freq), frequency = freq)
        lower = ts(lower, start(UCompObject$yFor, frequency = freq), frequency = freq)
    } else {
        nn = length(UCompObject$y) + 1
        upper = ts(upper, nn, frequency = 1)
        lower = ts(lower, nn, frequency = 1)
        UCompObject$yFor = ts(UCompObject$yFor, nn, frequency = 1)
        UCompObject$y = ts(UCompObject$y, 1, nn - 1)
    }
    # Finding model
    i = 2
    nPos = list(-1)
    while (nPos[[1]] == -1){
        lineI = UCompObject$table[i]
        nPos = gregexpr(pattern ='Model: ', lineI)
        i = i + 1
    }
    method = substr(lineI, nPos[[1]] + 7, nchar(lineI) -1)
    # Creating object
    fObject = list(model = method,
                   mean = UCompObject$yFor,
                   level = level,
                   x = UCompObject$y,
                   upper = upper,
                   lower = lower,
                   fitted = UCompObject$yFit,
                   method = method,
                   series = "Time Series",
                   residuals = UCompObject$v)
    return(structure(fObject, class = "forecast"))
}
#' @title predict.UComp
#' @description Forecasting using structural Unobseved Components models
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @return A list with components \code{pred} for the predictions and 
#'         \code{se} for standard errors (if \code{se.fit = TRUE})
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/eq/arma(0,0)")
#' f1 <- predict(m1)
#' @rdname predict.UComp
#' @export
predict.UComp = function(object, ...){
    out = list(pred = object$yFor,
               se = sqrt(object$yForV))
    return(out)
}
#' @title tsdiag.UComp
#' @description Diagnostic plots for UComp objects
#'
#' @details See help of \code{UC}.
#'
#' @param object Object of class \dQuote{UComp}.
#' @param gof.lag Maximum number of lags for pormanteau Ljung-Box test
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}}
#'          
#' @examples
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y, model = "llt/eq/arma(0,0)")
#' tsdiag(m1)
#' @rdname tsdiag.UComp
#' @export
tsdiag.UComp = function(object, gof.lag = NULL, ...){
    if (length(object$v) < 2){
        object = UCfilter(object)
    }
    # Number of u's
    nu = dim(object$u)[1]
    if (dim(object$u)[2] == 2){
        nu = 0
    }
    # Number of parameters
    nPar = length(object$p) - 1 + nu
    aux = list(residuals = object$v,
               sigma2 = 1,
               nobs = length(object$y) - nPar,
               coef = object$p,
               x = object$y,
               fitted = object$yFit)
    if (is.null(gof.lag)){
        tsdiag(structure(aux, class = "Arima"))
    } else {
        tsdiag(structure(aux, class = "Arima"), gof.lag)
    }
}

