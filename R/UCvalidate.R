#' @title UCvalidate
#' @description Shows a table of estimation and diagnostics results for UC models.
#' Equivalent to print or summary.
#' The table shows information in four sections:
#' Firstly, information about the model estimated, the relevant 
#' periods of the seasonal component included, and further information about
#' convergence.
#' Secondly, parameters with their names are provided, the asymptotic standard errors, 
#' the ratio of the two, and the gradient at the optimum. One asterisk indicates 
#' concentrated-out parameters and two asterisks signals parameters constrained during estimation.
#' Thirdly, information criteria and the value of the log-likelihood.
#' Finally, diagnostic statistics about innovations, namely, the Ljung-Box Q test of absense
#' of autocorrelation statistic for several lags, the Jarque-Bera gaussianity test, and a
#' standard ratio of variances test.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' @param printScreen print to screen or just return output table
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{table}{Estimation and validation table}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' m1 <- UC(log(AirPassengers))
#' m1 <- UCvalidate(m1)
#' @rdname UCvalidate
#' @export
UCvalidate = function(sys, printScreen = TRUE){
    if (is.ts(sys$y)){
        y = as.numeric(sys$y)
    } else {
        y = sys$y
    }
    if (is.ts(sys$u)){
        u = as.numeric(sys$u)
    } else {
        u = sys$u
    }
    if (any(is.na(sys$hidden$truePar))){
        if (printScreen){
            print(sys$table)
        }
        return(sys)
    }
    # Convert to R list
    #sys$periods = sys$hidden$periods0
    #sys$rhos = sys$hidden$rhos0
    rubbish = c(sys$hidden$d_t, sys$hidden$innVariance, sys$hidden$objFunValue, TRUE, sys$outlier, sys$arma, sys$iter)
    rubbish2 = cbind(sys$grad, sys$hidden$constPar, sys$hidden$typePar)
    rubbish3 = cbind(sys$hidden$ns, sys$hidden$nPar)
    output = UCompC("validate", y, u, sys$model, sys$periods, sys$rhos,
                    sys$h, sys$tTest, sys$criterion, sys$hidden$truePar, rubbish2, rubbish, sys$verbose, 
                    sys$stepwise, sys$hidden$estimOk, sys$p0, sys$v, sys$yFitV,
                    sys$hidden$nonStationaryTerms, rubbish3, sys$hidden$harmonics,
                    as.vector(sys$criteria), sys$hidden$cycleLimits, 
                    cbind(sys$hidden$beta, sys$hidden$betaV), sys$hidden$typeOutliers)
    sys$table = output$table
    if (is.ts(sys$y)){
        fY = frequency(sys$y)
        sY = start(sys$y, frequency = fY)
        aux = ts(matrix(NA, length(sys$y) - length(output$v) + 1, 1), sY, frequency = fY)
        sys$v = ts(output$v, end(aux), frequency = fY)
    } else {
        sys$v = output$v
    }
    if (printScreen){
        cat(output$table)
    }
    sys$covp = output$covp
    sys$p = as.vector(output$coef)
    # Parameter names from table
    nPar = length(sys$p)
    parNames = rep("", nPar)
    rowM = 2
    hyphen = 1
    i = 1
    while (hyphen < 4){
        lineI = sys$table[rowM]
        if (substr(lineI, 1, 1) == "-"){
            hyphen = hyphen + 1
        }
        if (hyphen > 2 && substr(lineI, 1, 1) != "-"){
            parNames[i] = substr(lineI, 1, gregexpr(pattern =':', lineI))
            i = i + 1
        }
        rowM = rowM + 1
    }
    rownames(sys$covp) = parNames[1 : dim(sys$covp)[1]]
    colnames(sys$covp) = parNames[1 : dim(sys$covp)[1]]
    names(sys$p) = parNames[1 : nPar]
    return(sys)
}
    