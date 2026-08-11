#' Summary method for a fitted \code{evgmrf} object
#'
#' @param object a fitted \code{evgmrf} object
#' @param ... not used
#'
#' @details
#' 
#' The main purpose of \code{summary.evgmrf} is to make clear the model that
#' has been fitted. It will give details of fixed effect specifications,
#' which parameters vary according to GMRFs, and the specifications of those
#' GMRFs, where relevant.
#'
#' @return A \code{summary.evgmrf} object
#'
#' @examples
#'
#' data(COtop5prcp)
#' COmxprcp <- COtop5prcp$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' summary(m_gev)
#'
#' @name summary.evgmrf
#'
#' @export
#' 
summary.evgmrf <- function(object, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  out <- list(family = toupper(object$family))
  out$fixed <- object$formula
  got_gmrf <- object$gmrf
  gmrf_type <- toupper(object$model)
  gmrf_order <- object$order
  names(got_gmrf) <- names(gmrf_type) <- names(gmrf_order) <- object$names$response
  out$GMRF <- list(GMRF = got_gmrf, GMRF_type = gmrf_type, GMRF_order = gmrf_order)
  out$call <- object$call
  out$logLik <- object$logLik
  class(out) <- "summary.evgmrf"
  out
}

#' @param x a \code{summary.evgmrf} object
#'
#' @rdname summary.evgmrf
#' 
#' @export
#' 
print.summary.evgmrf <- function(x, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  cat('\n Call: ', deparse(x$call), '\n\n')
  cat(' ** Parametric terms **\n')
  sapply(seq_along(x$fixed), function(i) cat(' - ', names(x$fixed[i]), ': ', gsub('~', '~ ', x$fixed[i]), '\n', sep = ''))
  cat('\n')
  cat(" ** GMRF terms **")
  GMRF_titles <- paste('- GMRF', c('present', 'type', 'order'))
  for (i in seq_along(x$GMRF)) {
    cat('\n', GMRF_titles[i], '\n')
    names(x$GMRF[[i]])[1] <- paste0(' -- ', names(x$GMRF[[i]])[1])
    cat(paste0(paste(names(x$GMRF[[i]]), ': ', x$GMRF[[i]], sep = ''), collapse = ', '), '')
  }
  cat("\n\n ** log-likelihood **\n")
  sapply(seq_along(x$logLik), function(i) cat(' - ', names(x$logLik)[i], ': ', x$logLik[[i]], '\n', sep = ''))
  cat('\n')
  invisible(x)
}

#' Print a fitted \code{evgmrf} object
#'
#' @param x a fitted \code{evgmrf} object
#' @param ... not used
#'
#' @return The call of the \code{evgmrf} object
#'
#' @examples
#'
#' data(COtop5prcp)
#' COmxprcp <- COtop5prcp$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' print(m_gev)
#'
#' @export
#' 
print.evgmrf <- function(x, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  print(x$call)
  invisible(x)
}

#' Log-likelihood from a fitted \code{evgmrf} object
#'
#' @param object a fitted \code{evgmrf} object
#' @param ... not used
#'
#' @return A scalar
#'
#' @examples
#' 
#' data(COtop5prcp)
#' COmxprcp <- COtop5prcp$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' logLik(m_gev)
#'
#' @export
#' 
logLik.evgmrf <- function(object, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  out <- object$logLik
  class(out) <- "logLik"
  out
}

#' Extract Model Fitted Values
#'
#' @param object a fitted \code{evgmrf} object
#' @param ... not used
#'
#' @examples
#'
#' data(COtop5prcp)
#' COmxprcp <- COtop5prcp$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' fitted(m_gev)
#'
#' @return Fitted values extracted from the object `object'.
#' 
#' @export
#' 
fitted.evgmrf <- function(object, ...) {
  predict(object)
}

#' Extract Model Coefficients
#'
#' @param object a fitted \code{evgmrf} object
#' @param ... not used
#'
#' @examples
#'
#' data(COtop5prcp)
#' COmxprcp <- COtop5prcp$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' coef(m_gev)
#'
#' @return Fitted values extracted from the object `object'.
#' 
#' @export
#' 
coef.evgmrf <- function(object, ...) {
  out <- split(object$beta, attr(object$beta, 'split'))
  names(out) <- object$names$response
  out
}
