#' @title TETSvalidate
#' @description Shows a table of estimation and diagnostics results for TOBIT TETS models
#'
#' @param m an object of type \code{TETS} created with \code{TETSmodel}
#' 
#' @return The same input object with the appropriate fields 
#' filled in, in particular:
#' \item{table}{Estimation and validation table}
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{TETS}}, \code{\link{TETSmodel}}, \code{\link{TETSvalidate}},
#'          \code{\link{TETScomponents}}
#'          
#' @examples
#' \dontrun{
#' m1 <- TETSmodel(log(gdp))
#' m1 <- TETSvalidate(m1)
#' }
#' @rdname TETSvalidate
#' @export
TETSvalidate = function(m){
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
    output = TETSc("validate", as.numeric(m$y), u, m$model, m$s, m$h,
                  m$criterion, m$armaIdent, m$identAll, m$forIntervals,
                  m$bootstrap, m$nSimul, FALSE, m$lambda,
                  m$alphaL, m$betaL, m$gammaL, m$phiL, m$p0, m$Ymin, m$Ymax)
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
    