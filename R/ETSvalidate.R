#' @title ETSvalidate
#' @description Shows a table of estimation and diagnostics results for UC models
#'
#' @param m an object of type \code{UComp} created with \code{UCmodel}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{table}{Estimation and validation table}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{ETS}}, \code{\link{ETSmodel}}, \code{\link{ETSvalidate}},
#'          \code{\link{ETScomponents}}
#'          
#' @examples
#' \dontrun{
#' m1 <- ETSmodel(log(gdp))
#' m1 <- ETSvalidate(m1)
#' }
#' @rdname ETSvalidate
#' @export
ETSvalidate = function(m){
    if (is.null(m$u))
        u = m$u
    else {
        if (is.vector(m$u)){
            u = matrix(m$u, 1, length(m$u))
        } else {
            nu = dim(m$u)
            u = as.numeric(m$u);
            u = matrix(u, nu[1], nu[2])
        }
    }
    output = ETSc("validate", as.numeric(m$y), u, m$model, m$s, m$h,
                  m$criterion, m$armaIdent, m$identAll, m$forIntervals,
                  m$bootstrap, m$nSimul, FALSE, m$lambda,
                  m$alphaL, m$betaL, m$gammaL, m$phiL, m$p0)
    if (is.ts(m$y))
        m$comp = ts(output$comp, start = start(m$y), frequency = frequency(m$y))
    else
        m$comp = output$comp
    m$table = output$table
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
    if (m$verbose)
        cat(m$table)
    colnames(m$comp) = strsplit(output$compNames, split = "/")[[1]]
    return(m)
}
    