#' @title print.TETS
#' @description Prints a TOBIT TETS object
#'
#' @details See help of \code{TETS}.
#'
#' @param x Object of class \dQuote{TETS}.
#' @param ... Additional inputs to handle the way to print output.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{TETS}}, \code{\link{TETSmodel}}, \code{\link{TETSvalidate}},
#'          \code{\link{TETScomponents}}, \code{\link{TETSestim}}
#'          
#' @examples
#' \dontrun{
#' m1 <- TETSmodel(log(gdp))
#' print(m1)
#' }
#' @rdname print
#' @export 
print.TETS = function(x, ...){
    if (length(x$table) == 0){
        x = TETSvalidate(x)
    } else {
        cat(x$table)
    }
}
#' @title summary.TETS
#' @description Prints a TOBIT TETS object on screen
#'
#' @param object Object of class \dQuote{TETS}.
#' @param ... Additional inputs to function.
#' 
#' @details See help of \code{TETS}.
#'
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{TETS}}, \code{\link{TETSmodel}}, \code{\link{TETSvalidate}},
#'          \code{\link{TETScomponents}}, \code{\link{TETSestim}}
#'          
#' @examples
#' \dontrun{
#' m1 <- TETSmodel(log(gdp))
#' summary(m1)
#' }
#' @rdname summary.TETS
#' @export 
summary.TETS = function(object, ...){
    print(object)
}
#' @title plot.TETS
#' @description Plot components of TETS object
#'
#' @details See help of \code{TETS}.
#'
#' @param x Object of class \dQuote{TETS}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{TETS}}, \code{\link{TETSmodel}}, \code{\link{TETSvalidate}},
#'          \code{\link{TETScomponents}}, \code{\link{TETSestim}}
#'          
#' @examples
#' \dontrun{
#' m1 <- TETSmodel(log(gdp))
#' plot(m1)
#' }
#' @rdname plot
#' @export 
plot.TETS = function(x, ...){
    if (length(x$comp) < 2){
        x = TETScomponents(x)
    }
    if (is.ts(x$comp)){
        plot(x$comp, main = "Time Series Decomposition")
    } else {
        plot(ts(x$comp, frequency = x$s),
             main = "Time Series Decomposition")
    }
}
#' @title fitted.TETS
#' @description Fitted output values of TETS object
#'
#' @details See help of \code{TETS}.
#'
#' @param object Object of class \dQuote{TETS}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{TETS}}, \code{\link{TETSmodel}}, \code{\link{TETSvalidate}},
#'          \code{\link{TETScomponents}}, \code{\link{TETSestim}}
#'          
#' @examples
#' \dontrun{
#' m1 <- TETSmodel(log(gdp))
#' fitted(m1)
#' }
#' @rdname fitted
#' @export 
fitted.TETS = function(object, ...){
    if (length(object$comp) < 2){
        object = TETScomponents(object)
    }
    return(object$com[, 2])
}
#' @title residuals.TETS
#' @description Residuals of TETS object
#'
#' @details See help of \code{TETS}.
#'
#' @param object Object of class \dQuote{TETS}.
#' @param ... Additional inputs to function.
#' 
#' @author Diego J. Pedregal
#' 
#' @seealso \code{\link{TETS}}, \code{\link{TETSmodel}}, \code{\link{TETSvalidate}},
#'          \code{\link{TETScomponents}}, \code{\link{TETSestim}}
#'          
#' @examples
#' \dontrun{
#' m1 <- TETSmodel(log(gdp))
#' residuals(m1)
#' }
#' @rdname residuals
#' @export 
residuals.TETS = function(object, ...){
    if (length(object$comp) < 2){
        object = TETScomponents(object)
    }
    return(object$com[, 1])
}
