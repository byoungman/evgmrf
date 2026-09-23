#' Predictions and standard errors from a fitted \code{evgmrf} object
#'
#' Obtains predicted values, extreme value quantiles, and associated standard 
#' errors from a fitted Gaussian Markov random field spatial extremes model.
#'
#' @param object A fitted \code{evgmrf} object.
#' @param type A character string specifying the prediction scale. Supported options 
#'   are \code{"link"} (the linear predictor scale) or \code{"response"} 
#'   (the inverse-link parameter scale). If \code{prob} is supplied, this is 
#'   automatically overridden to evaluate on the quantile scale. Defaults to \code{"link"}.
#' @param se.fit Logical; if `TRUE`, calculates associated standard errors for the 
#'   predictions alongside point estimates. Defaults to `FALSE`.
#' @param prob A scalar or vector of probabilities mapping to the target extreme value 
#'   quantiles to be estimated. Defaults to `NULL`.
#' @param nx,ny Positive integers giving the numbers of grid points in the
#'   first (\code{nx}) and second (\code{ny}) dimensions of the spatial grid. 
#'   Defaults to \code{object$nx} and \code{object$nx}, respectively. 
#' @param index A matrix containing spatial row and column coordinate index maps 
#'   used to reconstruct missing grid layouts. Defaults to \code{object$index}.
#' @param set2NA Set parameters to `NA` where there are no data. Defaults to `FALSE`.
#' @param simplify2array Logical; if `TRUE`, coerces and binds the underlying 
#'   prediction surfaces into a single multidimensional array matrix. Defaults to `FALSE`.
#' @param xid,yid Integer index vectors identifying localized grid subsets to extract. 
#'   Defaults to the entire structural dimensions tracked in \code{object$xid} and \code{object$yid}.
#' @param loop Logical; if `TRUE`, forces internal sparse solvers to iterate via 
#'   memory-conserving blocks rather than loading full dense matrices. Defaults to `TRUE`.
#' @param progress Logical; if `TRUE`, renders active text progress bars to the R console 
#'   tracking standard error computation loops. Defaults to the value of \code{loop}.
#' @param chunksize An integer determining the slice sizing processed concurrently 
#'   per standard error iterative loop. Defaults to `100`.
#' @param se.method A character string selecting the error propagation method. 
#'   Supported options are \code{"direct"} (analytical delta method solving) or 
#'   \code{"simulation"} (stochastic Monte Carlo sampling). Defaults to \code{"direct"}.
#' @param nsim An integer tracking random sampling pathways drawn when evaluating 
#'   stochastic errors under \code{se.method = "simulation"}. Defaults to `1000`.
#' @param decompose Logical; if `TRUE`, random additive terms inside Besag-York-Mollié 
#'   (BYM) formulations are broken down into individual spatial and random elements. 
#'   Defaults to `FALSE`.
#' @param random2zero Logical vector; first component corresponds to fitted values
#'   and second to standard errors. When `TRUE`, random additive terms or corresponding 
#'   standard errors inside Besag-York-Mollié (BYM) formulations are set to zero. 
#'   Defaults to `c(FALSE, FALSE)`.
#' @param drop.parametric Logical; if `TRUE` and \code{decompose = TRUE}, deletes 
#'   fixed parametric background terms from the final evaluated list. Defaults to `FALSE`.
#' @param openmp Logical; switches whether analytical solvers exploit shared memory 
#'   multi-core CPU parallelism extensions. Defaults to \code{object$control$openmp}.
#' @param threads An integer controlling maximum system CPU threads allocated 
#'   if \code{openmp = TRUE}. Defaults to \code{object$control$threads}.
#' @param ... Unused auxiliary flags passed along for generic matching alignment with 
#'   \code{\link[stats]{predict}}.
#' 
#' @details
#' Analytical standard errors are solved using Cholesky factorizations of the preconditioned 
#' Hessian matrix evaluated at the maximum likelihood convergence point. If Besag-York-Mollié 
#' (BYM2) configurations are detected, structural projections adjust the scaling matrix properties 
#' automatically. When calculating quantile-scale standard errors under \code{se.method = "direct"}, 
#' the function calls analytical gradients attached as a functional derivative attribute to 
#' \code{object$quantile0}.
#' 
#' See \code{evgmrf} for more details on \code{nx}, \code{ny} and \code{index}.
#' 
#' @references 
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Models. Journal of Statistical Software. \doi{10.18637/jss.v103.i03}
#'
#' @seealso \code{\link{evgmrf}}, \code{\link{simulate.evgmrf}}, \code{\link{family.evgmrf}}
#'
#' @return If \code{se.fit = FALSE}, the function returns a \code{list} or \code{array} of 
#'   point predictions. If \code{se.fit = TRUE}, it returns a nested \code{list} containing:
#'   \itemize{
#'     \item \code{fitted}: The point predictions or quantiles mapped to the grid layout.
#'     \item \code{se}: The associated standard errors matched element-for-element to the 
#'       shapes inside \code{fitted}.
#'   }
#' 
#' @examples
#' \dontrun{
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp)
#' 
#' # Evaluate location and scale predictions
#' link_preds <- predict(m_gev, type = "link")
#' 
#' # Calculate point estimates and standard errors for the 95th percentile
#' q_preds <- predict(m_gev, prob = 0.95, se.fit = TRUE, se.method = "direct")
#' }
#' 
#' @export
predict.evgmrf <- function(object, type = 'link', se.fit = FALSE, prob = NULL, index = NULL, nx = NULL,
                           ny = NULL, set2NA = FALSE, 
                           simplify2array = FALSE, xid = NULL, yid = NULL,
                           loop = TRUE, progress = FALSE, chunksize = 1e2, se.method = 'direct',
                           nsim = 1e3, decompose = FALSE, random2zero = c(FALSE, FALSE), drop.parametric = TRUE, 
                           openmp = object$control$openmp, threads = object$control$threads,  ...) {
  if (type != "link" & decompose)
    stop("Decomposed parameters only available for type = 'link'.")
  type0 <- type
  if (!is.null(prob))
    type <- 'quantile'
  if (type == 'quantile')
    type0 <- 'response'
  if (is.null(index))
    index <- object$index
  openmp <- object$likdata$openmp
  if (random2zero[1]) {
    if (decompose)
      stop("Can't have decompose = TRUE and random2zero = TRUE.")
    id_bym2 <- unlist(object$likdata$id_bym2)
    if (!any(id_bym2))
      stop('No random effects to set to zero.')
    object$beta[id_bym2] <- 0
  }
  out <- .fitted_values(object$beta, object$likdata, decompose)
  if (set2NA)
    out[, object$no_data] <- NA
  np <- nrow(out)
  if (!is.null(index)) {
    object$nx <- nx
    if (is.null(object$nx)) 
      object$nx <- max(index[, 1])
    object$ny <- ny
    if (is.null(object$ny))
      object$ny <- max(index[, 2])
    object$index <- index
    object$holes <- TRUE
    no_data <- matrix(TRUE, object$nx, object$ny)
    no_data[index] <- FALSE
    object$no_data <- as.vector(no_data)
  }
  if (is.null(yid))
    yid <- 1:object$ny
  if (is.null(xid))
    xid <- 1:object$nx
  if (!object$holes) {
    out <- lapply(1:np, function(i) matrix(out[i, ], object$nx))
  } else {
    if (!is.null(index)) {
      outm <- out
      out <- list()
      for (i in 1:np) {
        temp <- matrix(NA, object$nx, object$ny)
        temp[index] <- outm[i, ]
        out[[i]] <- temp
      }
    }
  }
  nms <- as.list(object$names[[type0]])
  if (decompose) {
    for (i in 1:length(nms))
      nms[[i]] <- paste(nms[[i]], object$par_type[[i]], sep = ': ')
  }
  names(out) <- unlist(nms)
  if (type %in% c('response', 'quantile')) {
    if (se.fit) {
      out0 <- out
      names(out0) <- unlist(object$names['link'])
    }
    for (i in 1:object$np) {
      out[[i]] <- object$unlink[[i]](out[[i]])
    }
  }
  out <- lapply(out, function(x) x[xid , yid, drop = FALSE])
  if (simplify2array) {
    out <- array(unlist(out), dim = c(dim(out[[1]]), length(out)))
  }
  if (type != 'link' & decompose) {
    names(out) <- gsub(': spatial', '', names(out))
    names(out) <- gsub(': random', '', names(out))
  }
  if (type == 'quantile') {
    nprob <- length(prob)
    temp <- list()
    for (i in 1:nprob) {
      out$p <- prob[i]
      temp[[i]] <- do.call(object$quantile, out)
    }
    out <- temp
    names(out) <- paste('q', prob, sep = '_')
    if (simplify2array)
      out <- array(out, dim = c(length(prob), object$nx, object$ny))
  } else {
    if (simplify2array) {
      out <- array(unlist(out), dim = c(object$np, object$nx, object$ny))
    }
  }
  if (set2NA)
    out <- lapply(out, function(x) {x[object$no_data] <- NA; x})
  if (se.fit) {
    if (progress) 
      cat('Calculating standard errors...\n')
    nv <- nrow(object$Hessian)
    dH <- Matrix::diag(object$diagHessian)
    cpH <- object$cholprecondHessian
    if (type %in% c('link', 'response')) {
      if (se.method == 'simulation') {
        if (progress) 
          pb <- txtProgressBar(min = 0, max = nsim / chunksize, style = 3)
        spl <- split(1:nsim, c(0:(nsim - 1)) %/% chunksize)
        se <- numeric(nv)
        for (j in 1:length(spl)) {
          z <- matrix(sample(c(-1, 1), length(spl[[j]]) * nv, replace = TRUE), nv)
          mat <- .solve_pchol(cpH, z)
          se <- se + rowSums(mat * mat)
          if (progress) setTxtProgressBar(pb, j)
        }
        se <- object$diagHessian * sqrt(se / nsim)
      } else {
        if (!openmp) {
          X <- Matrix::t(object$likdata$X)
            Xlc <- object$likdata$Xlc
            se_id <- rep(seq_along(Xlc), each = object$likdata$n)
            Xlc <- lapply(Xlc, function(x) x[sapply(x, ncol) > 0])
            if (!decompose) {
              Xlc <- lapply(Xlc, function(x) do.call(cbind, x))
              out_i <- seq_along(Xlc)
            } else {
              out_i <- rep(seq_along(Xlc), sapply(Xlc, length))
              Xlc <- unlist(Xlc, recursive = FALSE)
            }
            reps <- sapply(Xlc, ncol)
            X_col_id <- rep(seq_along(reps), reps)
            splitter <- rep(seq_along(reps), each = object$likdata$n)
            se <- rep(NA, length(splitter))
            if (progress)
              pb <- txtProgressBar(min = 0, max = length(se), style = 3)
            X0 <- object$likdata$X
            n_par <- ncol(X0)
            n_obs <- object$likdata$n
            ind0 <- seq_len(n_obs)

            for (i in seq_along(reps)) {
              
              cols_i <- which(X_col_id == i)
              X0i <- X0[, cols_i, drop = FALSE]
              n_i <- length(cols_i)
              rows_i <- se_id == out_i[i]
              ind <- ind0 + (i - 1L) * n_obs
              Ei <- Matrix::sparseMatrix(i = cols_i, j = seq_len(n_i), x = 1,
                dims = c(n_par, n_i), giveCsparse = TRUE)
              rows <- se_id == out_i[i]
              ind <- ind0 + (i - 1L) * n_obs
              X0r <- X0i[rows, , drop = FALSE]
              Bi <- dH * Matrix::tcrossprod(Ei, X0r)
              
              if (loop && (n_obs > chunksize)) {
                
                for (first in seq.int(1L, n_obs, by = chunksize)) {
                  
                  last <- min(first + chunksize - 1L, n_obs)
                  block <- first:last
                  cols_block <- cols_i[block]
                  t1 <- Matrix::solve(cpH, Bi[, block, drop = FALSE])
                  se[ind[block]] <- sqrt(Matrix::colSums(t1 * Bi[, block, drop = FALSE]))
                  
                  if (progress)
                    setTxtProgressBar(pb, tail(ind[block], 1))
                  
                }
              
              } else {
                
                se[ind] <- sqrt(Matrix::colSums(Bi * Matrix::solve(cpH, Bi)))
                
                if (progress)
                  setTxtProgressBar(pb, tail(ind, 1))
                
              }
              
            }

        } else {
          stop("Can't do openmp yet")
        }
      }
      se <- split(se, splitter)
      if (type == 'response') {
        for (i in 1:object$np) 
          se[[i]] <- se[[i]] * attr(object$unlink, 'deriv')[[i]](out0[[i]])
      }
      se <- lapply(se, matrix, object$nx, object$ny)
      names(se) <- names(out)
    }
    if (type == 'quantile') {
      if (progress)
        pb <- txtProgressBar(min = 0, max = length(prob) * object$n, style = 3)
      se <- list()
      ind0 <- (seq_len(object$np) - 1) * object$n
      for (k in 1:length(prob)) {
        sek <- numeric(object$n)
        if (se.method == 'simulation') {
          spl <- split(1:nsim, c(0:(nsim - 1)) %/% chunksize)
          for (j in 1:length(spl)) {
            ind <- spl[[j]]
            z <- matrix(rnorm(length(ind) * nv), ncol = length(ind))
            lst <- list()
            mat <- .solve_pchol(cpH, z)
            mat <- object$beta + dH * mat
            for (i in 1:object$np) {
              lst[[i]] <- mat[attr(object$beta, 'split') == i, , drop = FALSE]
              lst[[i]] <- object$X[[i]] %*% lst[[i]]
              lst[[i]] <- object$unlink[[i]](lst[[i]])
            }
            lst$p <- prob[k]
            lst <- as.matrix(do.call(object$quantile, lst))
            sek <- sek + rowSums((lst - as.vector(out[[k]]))^2)
            if (progress) 
              setTxtProgressBar(pb, k * object$n + max(ind))
          }
          sek <- sqrt(sek / nsim)
        } else {
          np0 <- object$likdata$np0
          outv <- lapply(out0, c)
          outv$p <- prob[k]
          J <- do.call(attr(object$quantile0, 'deriv'), outv)
          Xl <- object$likdata$Xl
          X <- Matrix::bdiag(Xl)
          dHX <- dH * Matrix::t(X)
          B <- Matrix::sparseMatrix(i = seq_len(np0 * object$n),
                                    j = rep(seq_len(object$n), np0),
                                    x = as.vector(J),
                                    dims = c(np0 * object$n, object$n))
          W   <- dHX %*% B                                  # p x n, column i = w_i
          if (!loop) {
            sek <- sqrt(Matrix::colSums(W * Matrix::solve(cpH, W)))
            if (progress)
              setTxtProgressBar(pb, k * object$n)
          } else {
            spl <- split(1:object$n, c(0:(object$n - 1)) %/% chunksize)
            sek <- numeric(object$n)
            for (j in seq_along(spl)) {
              ind <- spl[[j]]
              Wi <- W[ , ind, drop = FALSE]
              sek[ind] <- sqrt(Matrix::colSums(Wi * Matrix::solve(cpH, Wi)))
              if (progress) 
                setTxtProgressBar(pb, k * object$n + max(ind))
            }
          }
        }
        se[[k]] <- matrix(sek, object$nx, object$ny)[xid, yid, drop = FALSE]
      }
      names(se) <- paste('q', prob, sep = '_')
    }
    if (progress) {
      close(pb)
      cat('Done.\n')
    }
    if (set2NA)
      se <- lapply(se, function(x) {x[object$no_data] <- NA; x})
  }
  if (drop.parametric) {
    gonner <- replace(logical(length(out)), grep('parametric', names(out)), TRUE)
    out <- out[!gonner]
    if (se.fit) 
      se <- se[!gonner]
  }
  if (se.fit)
    out <- list(fitted = out, se = se)
  out
}

.fitted_values <- function(pars, likdata, decompose = FALSE) {
  pl <- split(pars, likdata$psplit)
  if (!decompose) {
    out <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  } else {
    Xlc <- likdata$Xlc
    Xlc <- lapply(likdata$Xlc, function(x) x[sapply(x, ncol) > 0])
    plc <- lapply(seq_along(pl), function(i) split(pl[[i]], rep(seq_along(Xlc[[i]]), sapply(Xlc[[i]], ncol))))
    Xlc <- unlist(Xlc, recursive = FALSE)
    plc <- unlist(plc, recursive = FALSE)
    out <- t(sapply(seq_along(plc), function(i) as.vector(Xlc[[i]] %*% plc[[i]])))
  }
  out
}
