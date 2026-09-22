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
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' summary(m_gev)
#'
#' @name summary.evgmrf
#'
#' @export
summary.evgmrf <- function(object, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  out <- list(family = toupper(object$family))
  out$fixed_formula <- object$formula
  got_gmrf <- object$gmrf
  gmrf_type <- toupper(object$model)
  gmrf_order <- object$order
  names(got_gmrf) <- names(gmrf_type) <- names(gmrf_order) <- object$names$response
  out$GMRF <- list(GMRF = got_gmrf, GMRF_type = gmrf_type, GMRF_order = gmrf_order)
  out$call <- object$call
  out$logLik <- object$logLik
  beta_table <- lapply(seq_along(object$fixed), function(.) 'None')
  if (length(unlist(object$fixed_id)) > 0) {
    beta_fixed_split <- rep(seq_along(object$fixed_id), sapply(object$fixed_id, length))
    temp_mat <- matrix(0, nrow = length(object$beta), ncol = length(unlist(object$fixed_id)))
    temp_ind <- cbind(unlist(object$fixed_id), 1:ncol(temp_mat))
    temp_mat[temp_ind] <- 1
    beta_fixed_ese <- sqrt(as.matrix(Matrix::solve(object$cholprecondHessian, temp_mat))[temp_ind])
    beta_fixed_ese <- beta_fixed_ese * Matrix::diag(object$diagHessian)[unlist(object$fixed_id)]
    for (i in seq_along(object$fixed_id)) {
      if (length(object$fixed[[i]]) > 0) 
        attr(object$fixed[[i]], 'e.s.e') <- beta_fixed_ese[beta_fixed_split == i]
    }
    out$beta_fixed <- object$fixed
    do_mat <- function(x) {
      out <- cbind(x, attr(x, 'e.s.e'))
      colnames(out) <- c('Estimate', 'Std. Error')
      out
    }
    for (i in seq_along(beta_table)) {
      if (length(out$beta_fixed[[i]]) > 0) {
        beta_table[[i]] <- do_mat(out$beta_fixed[[i]])
      }
    }
  }
  out$fixed_table <- beta_table
  out$names <- object$names
  out$hyper <- .pars2hyper(object$par, object$Qd$hyper, object$Qd$hyper_swap)
  names(out$hyper) <- object$names$response[out$GMRF[[1]]]
  for (i in seq_along(out$hyper)) for (j in seq_along(out$hyper[[i]])) {
    if (length(out$hyper[[i]][[j]]) > 1) {
      nij <- length(out$hyper[[i]][[j]])
      names(out$hyper[[i]][[j]]) <- paste0(names(out$hyper[[i]])[j], '_', 1:nij)
    }
  }
  class(out) <- "summary.evgmrf"
  out
}

#' @param x a \code{summary.evgmrf} object
#'
#' @rdname summary.evgmrf
#' 
#' @export
print.summary.evgmrf <- function(x, ...) {
  if (!missing(...)) warning("extra arguments discarded")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(' ** Parametric terms **\n')
  for (i in seq_along(x$fixed_formula)) {
    coef_name <- names(x$fixed_formula)[i]
    if (!is.matrix(x$fixed_table[[i]])) {
      cat(sprintf(" - %s: None\n", coef_name))
    } else {
      cat(sprintf(" - %s:\n", coef_name))
      xi <- x$fixed_table[[i]]
      xi_formatted <- format(xi, digits = 3)
      rownames(xi_formatted) <- paste0("    ", rownames(xi_formatted))
      print(xi_formatted, quote = FALSE)
      cat("\n")
    }
  }
  cat('\n')
  cat(" ** GMRF terms **")
  GMRF_titles <- paste('- GMRF', c('present', 'type', 'order', 'hyper-parameter estimates'))
  for (i in seq_along(x$GMRF)) {
    cat('\n', GMRF_titles[i], '\n')
    names(x$GMRF[[i]])[1] <- paste0(' -- ', names(x$GMRF[[i]])[1])
    cat(paste0(paste(names(x$GMRF[[i]]), ': ', x$GMRF[[i]], sep = ''), collapse = ', '), '')
  }
  cat('\n', tail(GMRF_titles, 1), '\n')
  for (i in seq_along(x$hyper)) {
    cat(paste0(' -- ', names(x$hyper)[i], ':\n', sep = ''))
    hyper_formatted <- unlist(lapply(x$hyper[[i]], format, digits = 4))
    hyper_vec <- paste0(unlist(lapply(x$hyper[[i]], names)), ' = ', hyper_formatted)
    cat(' ---', paste0(hyper_vec, collapse = ', '), '\n')
  }
  cat("\n ** log-likelihood **\n")
  sapply(seq_along(x$logLik), function(i) cat(' - ', names(x$logLik)[i], ':', c(' ', '   ', '  ')[i], x$logLik[[i]], '\n', sep = ''))
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
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' print(m_gev)
#'
#' @export
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
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' logLik(m_gev)
#'
#' @export
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
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' fitted(m_gev)
#'
#' @return Fitted values extracted from the object `object'.
#' 
#' @export
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
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = 'gev')
#' coef(m_gev)
#'
#' @return Fitted values extracted from the object `object'.
#' 
#' @export
coef.evgmrf <- function(object, ...) {
  out <- split(object$beta, attr(object$beta, 'split'))
  names(out) <- object$names$response
  out
}
