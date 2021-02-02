#' Print method for paramedic objects
#' 
#' The \code{print} method for paramedic objects displays a summary of the 
#' parameter estimates from the fitted model.
#' 
#' @param x The fitted object
#' @param ... Ignored
#' @return Returns \code{x}, invisibly.
#' @export
print.paramedic <- function(x, ...) {
    print(x$summary)
    invisible(x)
}

#' Summary method for paramedic objects
#' 
#' The \code{summary} method for paramedic objects displays a summary of the 
#' parameter estimates from the fitted model.
#' 
#' @param object The fitted object
#' @param ... Other arguments to \code{rstan::summary}
#' @return Returns a summary of \code{x}.
#' @export
summary.paramedic <- function(object, ...) {
    rstan::summary(object$stan_fit)
}