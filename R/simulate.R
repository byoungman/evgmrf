#' Simulations from a fitted \code{evgmrf} object.
#'
#' Produces plots of parameters or quantiles from a \code{evgmrf} fit.
#'
#' @param object a fitted \code{evgmrf} object
#' @param nsim to do
#' @param seed an integer giving the seed for simulations
#' @param type a character string giving the type of prediction sought; see Details. Defaults to \code{"link"}
#' @param prob a scalar or vector of probabilities for quantiles to be estimated if \code{type == "quantile"}; defaults to 0.5
#' @param simplify2array to do
#' @param ... unused
#' 
#' @details
#' 
#' To do.
#' 
#' @references 
#' 
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Models. Journal of Statistical Software. \doi{10.18637/jss.v103.i03}
#'
#' @examples
#'
#' # To follow
#'
#' @seealso \link{evgmrf} \link{predict.evgmrf}
#'
#' @return A \code{list} or \code{array}
#' 
#' @export
#' 
simulate.evgmrf <- function(object, nsim = 1, seed = NULL, type = 'link', prob = NULL,
                            simplify2array = TRUE, ...) {
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
  nH <- length(object$diagHessian)
  z <- matrix(rnorm(nsim * nH), ncol = nsim)
  mat <- .solve_pchol(object$cholprecondHessian, z)
  mat <- object$beta + object$diagHessian * mat
  lst <- list()
  for (i in 1:object$np) {
    lst[[i]] <- mat[attr(object$beta, 'split') == i, , drop = FALSE]
    lst[[i]] <- object$X[[i]] %*% lst[[i]]
  }
  if (type %in% c('response', 'quantile')) {
    for (i in 1:object$np) {
      lst[[i]] <- object$unlink[[i]](lst[[i]])
    }
  }
  if (simplify2array) {
    for (i in 1:object$np)
      lst[[i]] <- drop(array(lst[[i]], c(object$nx, object$ny, nsim)))
  }
  names(lst) <- unlist(object$names[type0])
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