#' Simulations from a fitted \code{evgmrf} object
#'
#' Simulates parameters or quantiles from the posterior or predictive 
#' distribution of a fitted \code{evgmrf} model.
#'
#' @param object A fitted \code{evgmrf} object.
#' @param nsim An integer specifying the number of simulation paths to generate. 
#'   Defaults to `1`.
#' @param seed An integer giving the seed for simulations, passed directly to 
#'   \code{\link{set.seed}}. If `NULL`, the current random number generator state 
#'   is preserved.
#' @param type A character string giving the type of simulation sought. Supported 
#'   options are \code{"link"} (the linear predictor scale) or \code{"response"} 
#'   (the inverse-link parameter scale). If \code{prob} is supplied, this parameter 
#'   is overridden to compute values on the quantile scale. Defaults to \code{"link"}.
#' @param prob A scalar or vector of probabilities mapping to quantiles to be 
#'   estimated. If supplied, the simulations will evaluate predictions on the 
#'   quantile scale. Defaults to `NULL`.
#' @param simplify2array Logical; if `TRUE` (the default), simulated paths are 
#'   coerced and dropped into spatial arrays dimensioned to match the grid layout 
#'   (\code{nx} by \code{ny} by \code{nsim}).
#' @param decompose Logical; if `TRUE`, structural additive terms within Besag-York-Mollié 
#'   (BYM) models are returned broken down into their individual spatial and non-spatial 
#'   components. Defaults to `FALSE`.
#' @param supernodal Logical; activates CHOLMOD supernodal sparse matrix factorization settings. 
#'   Defaults to `FALSE`.
#' @param ... Unused arguments. Passed along for generic compatibility with 
#'   \code{\link[stats]{simulate}}.
#' 
#' @details
#' The simulation routine draws multivariate normal samples using a precision-based 
#' Cholesky solver applied to the preconditioned Hessian matrix evaluated at the 
#' maximum likelihood estimates. These coefficient variations are then mapped 
#' through design matrices and inverse-link functions according to the target 
#' extreme value family definitions.
#' 
#' @references 
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Models. Journal of Statistical Software. \doi{10.18637/jss.v103.i03}
#'
#' @seealso \code{\link{evgmrf}}, \code{\link{predict.evgmrf}}
#'
#' @return A \code{list} or \code{array} of simulated paths structured according 
#'   to the choice of \code{simplify2array} and \code{decompose}.
#' 
#' @examples
#' \dontrun{
#' # Assuming a fitted model object 'm_gev' exists:
#' sim_paths <- simulate(m_gev, nsim = 100, type = "response")
#' }
#' 
#' @export
simulate.evgmrf <- function(object, nsim = 1, seed = NULL, type = 'link', prob = NULL,
                            simplify2array = TRUE, decompose = FALSE, supernodal = FALSE,
                            ...) {
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1) # initialize the RNG if necessary
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  type0 <- type
  if (!is.null(prob))
    type <- 'quantile'
  if (type == 'quantile')
    type0 <- 'response'
  if (decompose && type != 'link')
    stop("Decomposed simulations only available for type = 'link'.")
  dH <- Matrix::diag(object$diagHessian)
  nH <- length(dH)
  z <- matrix(rnorm(nsim * nH), ncol = nsim)
  if (!supernodal) {
    chol_factor <- .chol_analyze_simplicial(object$precondHessian)
    placeholder <- .chol_factorize_simplicial(chol_factor, object$precondHessian)
  } else {
    chol_factor <- .chol_analyze_supernodal(object$precondHessian)
    placeholder <- .chol_factorize_supernodal(chol_factor, object$precondHessian)
  }
  # draws with covariance H^{-1} = D A^{-1} D: half solve, not a full solve
  mat <- object$beta + dH * .chol_solve_Lt(chol_factor, z, supernodal)
  # one list element per component: per parameter, or per term if decompose
  lst <- list()
  for (i in 1:object$np) {
    bi <- mat[attr(object$beta, 'split') == i, , drop = FALSE]
    if (decompose) {
      Xlc <- object$likdata$Xlc[[i]]
      Xlc <- Xlc[sapply(Xlc, ncol) > 0]
      spl <- rep(seq_along(Xlc), sapply(Xlc, ncol))
      lst[[i]] <- lapply(seq_along(Xlc), function(j)
        as.matrix(Xlc[[j]] %*% bi[spl == j, , drop = FALSE]))
    } else {
      lst[[i]] <- list(as.matrix(object$X[[i]] %*% bi))
    }
  }
  lst <- unlist(lst, recursive = FALSE)
  nms <- as.list(object$names[[type0]])
  if (decompose) {
    for (i in seq_along(nms))
      nms[[i]] <- paste(nms[[i]], object$par_type[[i]], sep = ': ')
  }
  if (type %in% c('response', 'quantile')) {
    for (i in 1:object$np) {
      lst[[i]] <- object$unlink[[i]](lst[[i]])
    }
  }
  if (simplify2array) {
    for (i in seq_along(lst))
      lst[[i]] <- drop(array(lst[[i]], c(object$nx, object$ny, nsim)))
  }
  names(lst) <- unlist(nms)
  if (type == 'quantile') {
    nprob <- length(prob)
    if (nprob == 1) {
      lst$p <- prob[1]
      lst <- do.call(object$quantile, lst)
    } else {
      temp <- list()
      for (i in 1:nprob) {
        lst$p <- prob[i]
        temp[[i]] <- do.call(object$quantile, lst)
      }
      lst <- temp
    }
    names(lst) <- paste('q', prob, sep = '_')
  }
  lst
}
                            
.solve_pchol <- function(L_perm, b) {
  L <- as(L_perm, "Matrix")
  P <- as(L_perm, 'pMatrix')
  y <- Matrix::solve(Matrix::t(L), b, system = 'L')
  as.matrix(Matrix::crossprod(P, y))
}