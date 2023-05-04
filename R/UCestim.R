#' @title UCestim
#' @description Estimates and forecasts UC models
#'
#' @details \code{UCestim} estimates and forecasts a time series using an
#' UC model.
#' The optimization method is a BFGS quasi-Newton algorithm with a 
#' backtracking line search using Armijo conditions.
#' Parameter names in output table are the following:
#' \itemize{
#' \item Damping:   Damping factor for DT trend.
#' \item Level:     Variance of level disturbance.
#' \item Slope:     Variance of slope disturbance.
#' \item Rho(#):    Damping factor of cycle #.
#' \item Period(#): Estimated period of cycle #.
#' \item Var(#):    Variance of cycle #.
#' \item Seas(#):   Seasonal harmonic with period #.
#' \item Irregular: Variance of irregular component.
#' \item AR(#):     AR parameter of lag #.
#' \item MA(#):     MA parameter of lag #.
#' \item AO#:       Additive outlier in observation #.
#' \item LS#:       Level shift outlier in observation #.
#' \item SC#:       Slope change outlier in observation #.
#' \item Beta(#):   Beta parameter of input #.
#' \item Cnst:      Constant.
#' }
#' 
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @param sys an object of type \code{UComp} created with \code{UC}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \itemize{
#' \item p:        Estimated transformed parameters
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecast values of output
#' \item yForV:    Forecasted values variance
#' \item criteria: Value of criteria for estimated model
#' \item covp:     Covariance matrix of estimated transformed parameters
#' \item grad:     Gradient of log-likelihood at the optimum
#' \item iter:     Estimation iterations
#' }
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, 
#'          \code{\link{UCsmooth}}, \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' \dontrun{
#' m1 <- UCsetup(log(AirPassengers))
#' m1 <- UCestim(m1)
#' }
#' @rdname UCestim
#' @export
UCestim = function(sys){
    sys$table = NA
    sys$hidden$constPar = NA
    # Estimation
    rubbish = c(sys$hidden$d_t, sys$hidden$innVariance, sys$hidden$objFunValue, TRUE, 
                sys$outlier, sys$arma, sys$iter, sys$hidden$seas, sys$lambda,
                sys$hidden$MSOE, sys$hidden$PTSnames)
    rubbish2 = cbind(sys$grad, sys$hidden$constPar, sys$hidden$typePar)
    rubbish3 = cbind(sys$hidden$ns, sys$hidden$nPar)
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
    nu = dim(u)[2]
    kInitial = dim(u)[1]
    if (nu == 2){
        nu = length(sys$y) + sys$h
        kInitial = 0
    }
    output = UCompC("estimate", y, u, sys$model, sys$periods, sys$rhos,
                    sys$h, sys$tTest, sys$criterion, sys$hidden$truePar, rubbish2, rubbish, sys$verbose, 
                    sys$stepwise, sys$hidden$estimOk, sys$p0, sys$v, sys$yFitV,
                    sys$hidden$nonStationaryTerms, rubbish3, sys$hidden$harmonics,
                    as.vector(sys$criteria), sys$hidden$cycleLimits, 
                    cbind(sys$hidden$beta, sys$hidden$betaV), sys$hidden$typeOutliers,
                    sys$TVP, sys$trendOptions, sys$seasonalOptions, sys$irregularOptions)
    if (output$model == "error"){
        sys$model = "error"
        return(sys);
    }
    if (is.ts(sys$y)){
        fY = frequency(sys$y)
        sY = start(sys$y, frequency = fY)
        aux = ts(matrix(NA, length(sys$y) + 1, 1), sY, frequency = fY)
        if (length(output$yFor > 0)){
            sys$yFor = ts(output$yFor, end(aux), frequency = fY)
            sys$yForV = ts(output$yForV, end(aux), frequency = fY)
        }
    } else {
        if (length(output$yFor > 0)){
            sys$yFor = output$yFor
            sys$yForV = output$yForV
        }
    }
    # Convert to R list
    #sys$p = output$p[, 2]
    sys$hidden$truePar = output$p[, 1]
    sys$p0 = output$p0
    if (grepl("?", sys$model, fixed = TRUE)){
        sys$model = output$model
    }
    n = length(sys$hidden$truePar)
    rubbish2 = matrix(output$rubbish2, n, 3)
    sys$grad = rubbish2[, 1]
    sys$hidden$constPar = rubbish2[, 2]
    sys$hidden$typePar = rubbish2[, 3]
    sys$hidden$cycleLimits = matrix(output$cycleLimits, 
                                    length(output$cycleLimits) / 2, 2)
    sys$hidden$d_t = output$rubbish[1]
    sys$hidden$innVariance = output$rubbish[2]
    sys$hidden$objFunValue = output$rubbish[3]
    sys$iter = output$rubbish[6]
    sys$h = output$rubbish[7]
    sys$lambda = output$rubbish[8]
    betas = matrix(output$betas, length(output$betas) / 2, 2)
    sys$hidden$beta = betas[, 1]
    sys$hidden$betaV = betas[, 2]
    sys$periods = output$periods
    sys$rhos = output$rhos
    sys$hidden$estimOk = output$estimOk
    sys$hidden$nonStationaryTerms = output$nonStationaryTerms
    rubbish3 = matrix(output$rubbish3, 7, 2)
    sys$hidden$ns = rubbish3[, 1]
    sys$hidden$nPar = rubbish3[, 2]
    sys$hidden$harmonics = output$harmonics
    criteria = output$criteria;
    sys$criteria = matrix(criteria, 1, 4)
    colnames(sys$criteria) = c("LLIK", "AIC", "BIC", "AICc")
    sys$u = output$u
    if (!is.na(sys$outlier) && !is.null(u)){
        nu = length(sys$y) + sys$h;
        k = length(output$u) / nu
        nOut = k - kInitial
        if (nOut > 0){
            sys$u = matrix(output$u, k, nu)
            sys$hidden$typeOutliers = output$typeOutliers
        }
    }
    return(sys)
}
    