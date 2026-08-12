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
#' @param index A matrix containing spatial row and column coordinate index maps 
#'   used to reconstruct missing grid layouts. Defaults to \code{object$index}.
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
#' @param drop.parametric Logical; if `TRUE` and \code{decompose = TRUE}, deletes 
#'   fixed parametric background terms from the final evaluated list. Defaults to `FALSE`.
#' @param sdif A numeric standard deviation inflation factor modifier scaled across 
#'   the linear predictor error structures. Defaults to `1`.
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
#' data(COtop5prcp)
#' COmxprcp <- COtop5prcp$prcp[, 1, , ]
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
predict.evgmrf <- function(object, type = 'link', se.fit = FALSE, prob = NULL, index = object$index, 
                           simplify2array = FALSE, xid = object$xid, yid = object$yid,
                           loop = TRUE, progress = loop, chunksize = 1e2, se.method = 'direct',
                           nsim = 1e3, decompose = FALSE, drop.parametric = FALSE, sdif = 1, 
                           openmp = object$control$openmp, threads = object$control$threads, ...) {
  if (se.fit & decompose)
    stop("Decomposed parameters can only be plotted for type = 'link'.")
  type0 <- type
  if (!is.null(prob))
    type <- 'quantile'
  if (type == 'quantile')
    type0 <- 'response'
  openmp <- object$likdata$openmp
  out <- .fitted_values(object$beta, object$likdata, decompose)
  np <- nrow(out)
  if (!object$holes) {
    out <- lapply(1:np, function(i) matrix(out[i, ], object$nx))
  } else {
    if (!is.null(index)) {
      outm <- out
      out <- list()
      for (i in 1:np) {
        nx <- ifelse(is.null(object$nx), max(index[, 1]), object$nx)
        ny <- ifelse(is.null(object$ny), max(index[, 2]), object$ny)
        temp <- matrix(NA, nx, ny)
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
  if (decompose & drop.parametric)
    out <- out[unlist(object$par_type) != 'parametric']
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
  if (se.fit) {
    nv <- nrow(object$Hessian)
    dH <- object$diagHessian
    cpH <- object$cholprecondHessian
    if (progress) cat('Calculating standard errors...\n')
    id_bym2 <- object$likdata$id_bym2
    if (any(unlist(id_bym2))) {
      id_bym2_mat <- lapply(id_bym2, function(x) length(x) / object$n)
      XX <- lapply(id_bym2_mat, function(i) lapply(1:i, function(.) Matrix::Diagonal(n = object$n, x = 1)))
      XX <- lapply(XX, function(x) do.call(rbind, x))
      XX <- Matrix::bdiag(XX)
      H <- Matrix::crossprod(XX, object$Hessian) %*% XX
      nv <- nrow(H)
      dH <- pmax(Matrix::diag(H), 1e-8)
      D <- Matrix::Diagonal(nrow(H), 1 / sqrt(dH))
      H <- D %*% H %*% D
      cpH <- Matrix::Cholesky(H, super = object$likdata$control$super, LDL = FALSE)
    }
    if (type %in% c('link', 'response')) {
      if (se.method == 'simulation') {
        if (progress) 
          pb <- txtProgressBar(min = 0, max = nsim / chunksize, style = 3)
        spl <- split(1:nsim, c(0:(nsim - 1)) %/% chunksize)
        se <- numeric(nv)
        for (j in 1:length(spl)) {
          z <- matrix(sample(c(-sdif, sdif), length(spl[[j]]) * nv, replace = TRUE), nv)
          mat <- .solve_pchol(cpH, z)
          se <- se + rowSums(mat * mat)
          if (progress) setTxtProgressBar(pb, j)
        }
        se <- object$diagHessian * sqrt(se / nsim)
      } else {
        if (!openmp) {
        if (progress)
          pb <- txtProgressBar(min = 0, max = nv / chunksize, style = 3)
        if (!loop) {
          t1 <- Matrix::solve(cpH, Matrix::Diagonal(nv))
          se <- dH * sqrt(Matrix::diag(t1))
        } else {
          se <- rep(NA, nv)
          spl <- split(1:nv, c(0:(nv - 1)) %/% chunksize)
          for (j in 1:length(spl)) {
            ind <- spl[[j]]
            Ej <- Matrix::sparseMatrix(i = ind, j = seq_along(ind), x = 1, dims = c(nv, length(ind)))
            v <- Matrix::solve(cpH, Ej)
            se[ind] <- sdif * dH[ind] * sqrt(v[cbind(ind, seq_along(ind))])
            if (progress) setTxtProgressBar(pb, j)
          }
        }
      } else {
        se <- dH * .chol_idiag_omp(object$precondHessian, object$control$threads)
      }
      }
      se <- split(se, rep(1:object$np, each = object$n))
      if (type == 'response') {
        for (i in 1:object$np) 
          se[[i]] <- se[[i]] * attr(object$unlink, 'deriv')[[i]](out0[[i]])
      }
      se <- lapply(se, matrix, object$nx, object$ny)
      names(se) <- names(out)
    }
    if (type == 'quantile') {
      if (progress)
        pb <- txtProgressBar(min = 0, max = object$n / chunksize, style = 3)
      se <- list()
      ind0 <- (seq_len(object$np) - 1) * object$n
      for (k in 1:length(prob)) {
        sek <- numeric(object$n)
        if (se.method == 'simulation') {
          if (progress) 
            pb <- txtProgressBar(min = 0, max = nsim / chunksize, style = 3)
          spl <- split(1:nsim, c(0:(nsim - 1)) %/% chunksize)
          for (j in 1:length(spl)) {
            z <- matrix(rnorm(length(spl[[j]]) * nv), ncol = length(spl[[j]]))
            lst <- list()
            mat <- .solve_pchol(cpH, sdif * z)
            mat <- object$beta + dH * mat
            for (i in 1:object$np) {
              lst[[i]] <- mat[attr(object$beta, 'split') == i, , drop = FALSE]
              lst[[i]] <- object$X[[i]] %*% lst[[i]]
              lst[[i]] <- object$unlink[[i]](lst[[i]])
            }
            lst$p <- prob[k]
            lst <- as.matrix(do.call(object$quantile, lst))
            sek <- sek + rowSums((lst - as.vector(out[[k]]))^2)
            if (progress) setTxtProgressBar(pb, j)
          }
          sek <- sqrt(sek / nsim)
        } else {
          out0$p <- prob[k]
          J <- do.call(attr(object$quantile0, 'deriv'), out0)
          J <- dH * matrix(J, ncol = object$np)
          if (loop) {
            spl <- split(1:object$n, c(0:(object$n - 1)) %/% chunksize)
            for (j in 1:length(spl)) {
              ind <- spl[[j]]
              Jj <- as.vector(t(J[ind, ]))
              indr <- as.integer(outer(ind0, ind, FUN = '+'))
              indc <- rep(seq_along(ind), each = object$np)
              Ej <- Matrix::sparseMatrix(i = indr, j = indc, x = Jj, dims = c(nv, length(ind)))
              temp <- Matrix::solve(cpH, Ej)[indr, , drop = FALSE]
              sek[ind] <- sdif * sqrt(Matrix::colSums(Ej[indr, , drop = FALSE ] * temp))
              if (progress) setTxtProgressBar(pb, j)
            }
          } else {
            indr <- as.integer(outer(ind0, 1:object$n, FUN = '+'))
            indc <- rep(1:object$n, each = object$np)
            EJ <- Matrix::sparseMatrix(i = indr, j = indc, x = as.vector(t(J)), dims = c(nv, object$n))
            sek <- sdif * sqrt(Matrix::colSums(EJ * Matrix::solve(cpH, EJ)))
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
    out <- list(fitted = out, se = se)  
  }
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
