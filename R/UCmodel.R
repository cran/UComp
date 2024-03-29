#' @title UCsetup
#' @description Sets up UC general univariate models
#'
#' @details See help of \code{UC}.
#'
#' @param y a time series to forecast (it may be either a numerical vector or
#' a time series object). This is the only input required. If a vector, the additional
#' input \code{periods} should be supplied compulsorily (see below).
#' @param u a matrix of external regressors included only in the observation equation. 
#' (it may be either a numerical vector or a time series object). If the output wanted 
#' to be forecast, matrix \code{u} should contain future values for inputs.
#' @param model the model to estimate. It is a single string indicating the type of 
#' model for each component. It allows two formats "trend/seasonal/irregular" or 
#' "trend/cycle/seasonal/irregular". The possibilities available for each component are:
#' \itemize{
#' \item Trend: ? / none / rw / irw / llt / dt / td; 
#' 
#' \item Seasonal: ? / none / equal / different;
#' 
#' \item Irregular: ? / none / arma(0, 0) / arma(p, q) - with p and q integer positive orders;
#'     
#' \item Cycles: ? / none / combination of positive or negative numbers. Positive numbers fix
#' the period of the cycle while negative values estimate the period taking as initial
#' condition the absolute value of the period supplied. Several cycles with positive or negative values are possible
#' and if a question mark is included, the model test for the existence of the cycles
#' specified. The following are valid examples with different meanings: 48, 48?, -48, -48?,
#' 48+60, -48+60, -48-60, 48-60, 48+60?, -48+60?, -48-60?, 48-60?.
#' }
#' @param outlier critical level of outlier tests. If NA it does not carry out any 
#' outlier detection (default). A positive value indicates the critical minimum
#' t test for outlier detection in any model during identification. Three types of outliers are
#' identified, namely Additive Outliers (AO), Level Shifts (LS) and Slope Change (SC).
#' @param stepwise stepwise identification procedure (TRUE / FALSE).
#' @param tTest augmented Dickey Fuller test for unit roots used in stepwise algorithm (TRUE / FALSE). 
#' The number of models to search for is reduced, depending on the result of this test.
#' @param p0 initial parameter vector for optimisation search.
#' @param h forecast horizon. If the model includes inputs h is not used, the lenght of u is used instead.
#' @param lambda Box-Cox transformation lambda, NULL for automatic estimation
#' @param criterion information criterion for identification ("aic", "bic" or "aicc").
#' @param periods vector of fundamental period and harmonics required.
#' @param verbose intermediate results shown about progress of estimation (TRUE / FALSE).
#' @param arma check for arma models for irregular components (TRUE / FALSE).
#' @param TVP vector of zeros and ones to indicate TVP parameters.
#' @param trendOptions trend models to select amongst (e.g., "rw/llt").
#' @param seasonalOptions seasonal models to select amongst (e.g., "none/differentt").
#' @param irregularOptions irregular models to select amongst (e.g., "none/arma(0,1)").
#' 
#' @author Diego J. Pedregal
#' 
#' @return An object of class \code{UComp}. It is a list with fields including all the inputs and
#'         the fields listed below as outputs. All the functions in this package fill in
#'         part of the fields of any \code{UComp} object as specified in what follows (function 
#'         \code{UC} fills in all of them at once):
#' 
#' After running \code{UCmodel} or \code{UCestim}:
#' \itemize{
#' \item p:        Estimated parameters
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecasted values of output
#' \item yForV:    Variance of forecasts
#' \item criteria: Value of criteria for estimated model
#' \item iter:     Number of iterations in estimation
#' \item grad:     Gradient at estimated parameters
#' \item covp:     Covariance matrix of parameters
#' }
#' 
#' After running \code{UCvalidate}:
#' \itemize{
#' \item table: Estimation and validation table
#' }
#' 
#' After running \code{UCcomponents}:
#' \itemize{
#' \item comp:      Estimated components in matrix form
#' \item compV:     Estimated components variance in matrix form
#' }
#' 
#' After running \code{UCfilter}, \code{UCsmooth} or  \code{UCdisturb}:
#' \itemize{
#' \item yFit:  Fitted values of output
#' \item yFitV: Estimated fitted values variance
#' \item a:     State estimates
#' \item P:     Variance of state estimates
#' \item aFor:  Forecasts of states
#' \item PFor:  Forecasts of states variances
#' }
#' 
#' After running \code{UCdisturb}:
#' \itemize{
#' \item eta: State perturbations estimates
#' \item eps: Observed perturbations estimates
#' }
#' 
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCmodel}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' \dontrun{
#' y <- log(AirPassengers)
#' m1 <- UCsetup(y)
#' m1 <- UCsetup(y, outlier = 4)
#' m1 <- UCsetup(y, model = "llt/equal/arma(0,0)")
#' m1 <- UCsetup(y, model = "?/?/?/?")
#' m1 <- UCsetup(y, model = "llt/?/equal/?", outlier = 4)
#' }
#' @rdname UCsetup
#' @export
UCsetup = function(y, u = NULL, model = "?/none/?/?", h = 9999, lambda = 1, outlier = 9999, tTest = FALSE, criterion = "aic",
                   periods = NA, verbose = FALSE, stepwise = FALSE, p0 = -9999.9, arma = FALSE,
                   TVP = NULL, trendOptions = "none/rw/llt/dt", seasonalOptions = "none/equal/different", irregularOptions = "none/arma(0,0)"){
    # Converting u vector to matrix
    if (length(size(u)) == 1 && size(u) > 0){
        u = matrix(u, 1, length(u))
    }
    if (length(size(u)) == 2 && size(u)[1] > size(u)[2])
        u = t(u)
    if (is.null(TVP) && !is.null(u))
        TVP = matrix(0, nrow(u))
    if (is.null(TVP))
        TVP = -9999.99
    if (is.null(u)){
        u = matrix(0, 1, 2)
    }
    if (is.null(lambda))
        lambda = 9999.9
    # Checking periods
    if (is.ts(y) && any(is.na(periods)) && frequency(y) > 1){
        periods = frequency(y) / (1 : floor(frequency(y) / 2))
    } else if (is.ts(y) && any(is.na(periods))){
        periods = 1
    } else if (!is.ts(y) && any(is.na(periods))){
        stop("Input \"periods\" should be supplied!!")
    }
    rhos = rep(1, length(periods))
    out = list(y = y,
               u = u,
               model = model,
               h = h,
               lambda = lambda,
               # Other less important
               arma =  arma,
               outlier = -abs(outlier),
               tTest = tTest,
               criterion = criterion,
               periods = periods,
               rhos = rhos,
               verbose = verbose,
               stepwise = stepwise,
               p0 = p0,
               criteria = NA,
               TVP = TVP,
               trendOptions = trendOptions,
               seasonalOptions = seasonalOptions,
               irregularOptions = irregularOptions,
               # Outputs
               comp = NA,
               compPlus = NA,
               compMinus = NA,
               p = NA,
               covp = NA,
               grad = NA,
               v= NA,
               yFit = NA,
               yFor = NA,
               yFitV = NA,
               yForV = NA,
               a = NA,
               P = NA,
               eta = NA,
               eps = NA,
               table = NA,
               iter = 0,
               # Other variables
               hidden = list(d_t = 0,
                             estimOk = "Not estimated",
                             objFunValue = 0,
                             innVariance = 1,
                             nonStationaryTerms = NA,
                             ns = NA,
                             nPar = NA,
                             harmonics = 0,
                             constPar = NA,
                             typePar = NA,
                             cycleLimits = NA,
                             typeOutliers = matrix(-1, 1, 2),
                             truePar = NA,
                             beta = NA,
                             betaV = NA,
                             seas = frequency(y),
                             MSOE = FALSE,
                             PTSnames = FALSE))
    return(structure(out, class = "UComp"))
}

#' @title UCmodel
#' @description Estimates and forecasts UC general univariate models
#'
#' @details \code{UCmodel} is a function for modelling and forecasting univariate
#' time series according to Unobserved Components models (UC). 
#' It sets up the model with a number of control variables that
#' govern the way the rest of functions in the package work. It also estimates 
#' the model parameters by Maximum Likelihood and forecasts the data.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @inheritParams UCsetup
#' 
#' @return An object of class \code{UComp}. It is a list with fields including all the inputs and
#'         the fields listed below as outputs. All the functions in this package fill in
#'         part of the fields of any \code{UComp} object as specified in what follows (function 
#'         \code{UC} fills in all of them at once):
#' 
#' After running \code{UCmodel} or \code{UCestim}:
#' \itemize{
#' \item p:        Estimated parameters
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecasted values of output
#' \item yForV:    Forecasted values +- one standard error
#' \item criteria: Value of criteria for estimated model
#' \item iter:     Number of iterations in estimation
#' \item grad:     Gradient at estimated parameters
#' \item covp:     Covariance matrix of parameters
#' }
#' 
#' After running \code{UCvalidate}:
#' \itemize{
#' \item table: Estimation and validation table
#' }
#' 
#' After running \code{UCcomponents}:
#' \itemize{
#' \item comp:  Estimated components in matrix form
#' \item compV: Estimated components variance in matrix form
#' }
#' 
#' After running \code{UCfilter}, \code{UCsmooth} or  \code{UCdisturb}:
#' \itemize{
#' \item yFit:  Fitted values of output
#' \item yFitV: Variance of fitted values of output
#' \item a:     State estimates
#' \item P:     Variance of state estimates
#' \item aFor:  Forecasts of states
#' \item PFor:  Forecasts of states variances
#' }
#' 
#' After running \code{UCdisturb}:
#' \itemize{
#' \item eta: State perturbations estimates
#' \item eps: Observed perturbations estimates
#' }
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' \dontrun{
#' y <- log(AirPassengers)
#' m1 <- UCmodel(y)
#' m1 <- UCmodel(y, model = "llt/equal/arma(0,0)")
#' }
#' @rdname UCmodel
#' @export
UCmodel = function(y, u = NULL, model = "?/none/?/?", h = 9999, lambda = 1, outlier = 9999, tTest = FALSE, criterion = "aic",
                   periods = NA, verbose = FALSE, stepwise = FALSE, p0 = -9999.9, arma = TRUE,
                   TVP = NULL, trendOptions = "none/rw/llt/dt", seasonalOptions = "none/equal/different", irregularOptions = "none/arma(0,0)"){
    m = UCsetup(y, u, model, h, lambda, outlier, tTest, criterion, 
                periods, verbose, stepwise, p0, arma,
                TVP, trendOptions, seasonalOptions, irregularOptions)
    m = UCestim(m)
    return(m)
}
#' @title UC
#' @description Runs all relevant functions for UC modelling
#'
#' @details \code{UC} is a function for modelling and forecasting univariate
#' time series according to Unobserved Components models (UC). 
#' It sets up the model with a number of control variables that
#' govern the way the rest of functions in the package work. It also estimates 
#' the model parameters by Maximum Likelihood, forecasts the data, performs smoothing,
#' estimates model disturbances, estimates components and shows statistical diagnostics.
#' Standard methods applicable to UComp objects are print, summary, plot,
#' fitted, residuals, logLik, AIC, BIC, coef, predict, tsdiag.
#'
#' @inheritParams UCsetup
#' 
#' @author Diego J. Pedregal
#' 
#' @return An object of class \code{UComp}. It is a list with fields including all the inputs and
#'         the fields listed below as outputs. All the functions in this package fill in
#'         part of the fields of any \code{UComp} object as specified in what follows (function 
#'         \code{UC} fills in all of them at once):
#' 
#' After running \code{UCmodel} or \code{UCestim}:
#' \itemize{
#' \item p:        Estimated parameters
#' \item v:        Estimated innovations (white noise in correctly specified models)
#' \item yFor:     Forecasted values of output
#' \item yForV:    Forecasted values +- one standard error
#' \item criteria: Value of criteria for estimated model
#' \item iter:     Number of iterations in estimation
#' \item grad:     Gradient at estimated parameters
#' \item covp:     Covariance matrix of parameters
#' }
#' 
#' After running \code{UCvalidate}:
#' \itemize{
#' \item table: Estimation and validation table
#' }
#' 
#' After running \code{UCcomponents}:
#' \itemize{
#' \item comp:  Estimated components in matrix form
#' \item compV: Estimated components variance in matrix form
#' }
#' 
#' After running \code{UCfilter}, \code{UCsmooth} or  \code{UCdisturb}:
#' \itemize{
#' \item yFit:  Fitted values of output
#' \item yFitV: Variance of fitted values of output
#' \item a:     State estimates
#' \item P:     Variance of state estimates
#' \item aFor:  Forecasts of states
#' \item PFor:  Forecasts of states variances
#' }
#' 
#' After running \code{UCdisturb}:
#' \itemize{
#' \item eta: State perturbations estimates
#' \item eps: Observed perturbations estimates
#' }
#' 
#' @seealso \code{\link{UC}}, \code{\link{UCvalidate}}, \code{\link{UCfilter}}, \code{\link{UCsmooth}}, 
#'          \code{\link{UCdisturb}}, \code{\link{UCcomponents}},
#'          \code{\link{UChp}}
#'          
#' @examples
#' \dontrun{
#' y <- log(AirPassengers)
#' m1 <- UC(y)
#' m1 <- UC(y, model = "llt/different/arma(0,0)")
#' }
#' @rdname UC
#' @export
UC = function(y, u = NULL, model = "?/none/?/?", h = 9999, lambda = 1, outlier = 9999, tTest = FALSE, criterion = "aic",
              periods = NA, verbose = FALSE, stepwise = FALSE, p0 = -9999.9, arma = TRUE,
              TVP = NULL, trendOptions = "none/rw/llt/dt", seasonalOptions = "none/equal/different", irregularOptions = "none/arma(0,0)"){
    m = UCsetup(y, u, model, h, lambda, outlier, tTest, criterion, 
                periods, verbose, stepwise, p0, arma,
                TVP, trendOptions, seasonalOptions, irregularOptions)
    m = UCestim(m)
    if (m$model == "error")
        return(m)
    m = UCvalidate(m, verbose)
    m = UCdisturb(m)
    m = UCsmooth(m)
    m = UCcomponents(m)
    return(m)
}
