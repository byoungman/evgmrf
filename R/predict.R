#' Plots from a fitted \code{evgmrf} object.
#'
#' Produces plots of parameters or quantiles from a \code{evgmrf} fit.
#'
#' @param object a fitted \code{evgmrf} object
#' @param type a character string giving the type of prediction sought; see Details. Defaults to \code{"link"}
#' @param se.fit to do
#' @param prob a scalar or vector of probabilities for quantiles to be estimated if \code{type == "quantile"}; defaults to 0.5
#' @param index to do
#' @param simplify2array to do
#' @param loop to do
#' @param progress to do
#' @param chunksize to do
#' @param se.method to do
#' @param nsim to do
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
predict.evgmrf <- function(object, type = 'link', se.fit = FALSE, prob = NULL, index = object$index, 
                           simplify2array = FALSE,
                           loop = TRUE, progress = loop, chunksize = 1e2, se.method = 'direct',
                           nsim = 1e3, ...) {
  type0 <- type
  if (!is.null(prob))
    type <- 'quantile'
  if (type == 'quantile')
    type0 <- 'response'
  out <- .fitted_values(object$beta, object$likdata)
  np <- object$likdata$np0
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
  # if (type == 'response') {
  #   for (i in 1:np) out[[i]] <- object$likfns$trans[[i]](out[[i]])
  # }
  if (type %in% c('response', 'quantile')) {
    if (se.fit) {
      out0 <- out
      names(out0) <- unlist(object$names['link'])
    }
    for (i in 1:object$np) {
      out[[i]] <- object$unlink[[i]](out[[i]])
    }
  }
  names(out) <- unlist(object$names[type0])
  if (simplify2array) {
    out <- array(unlist(out), dim = c(dim(out[[1]]), length(out)))
    # for (i in 1:object$np)
    #   out[[i]] <- drop(array(lst[[i]], c(object$nx, object$ny, nsim)))
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
    if (progress) cat('Calculating standard errors...\n')
    nv <- nrow(object$Hessian)
    if (type %in% c('link', 'response')) {
      if (se.method == 'simulation') {
        if (progress) 
          pb <- txtProgressBar(min = 0, max = nsim / chunksize, style = 3)
        spl <- split(1:nsim, c(0:(nsim - 1)) %/% chunksize)
        se <- numeric(nv)
        for (j in 1:length(spl)) {
          z <- matrix(sample(c(-1, 1), length(spl[[j]]) * nv, replace = TRUE), nv)
          mat <- .solve_pchol(object$cholprecondHessian, z)
          se <- se + rowSums(mat * mat)
          if (progress) setTxtProgressBar(pb, j)
        }
        se <- object$diagHessian * sqrt(se / nsim)
      } else {
        if (progress)
          pb <- txtProgressBar(min = 0, max = nv / chunksize, style = 3)
        if (!loop) {
          t1 <- Matrix::solve(object$cholprecondHessian, Matrix::Diagonal(nv))
          se <- object$diagHessian * sqrt(Matrix::diag(t1))
        } else {
          se <- rep(NA, nv)
        # for (j in 1:nv) {
        #   ej <- Matrix::sparseVector(i = j, x = 1, length = nv)
        #   v <- as.vector(Matrix::solve(object$cholprecondHessian, ej))
        #   se[j] <- object$diagHessian[j] * sqrt(v[j])
        #   if (progress) setTxtProgressBar(pb, j)
        # }
          spl <- split(1:nv, c(0:(nv - 1)) %/% chunksize)
          for (j in 1:length(spl)) {
            ind <- spl[[j]]
            Ej <- Matrix::sparseMatrix(i = ind, j = seq_along(ind), x = 1, dims = c(nv, length(ind)))
            v <- Matrix::solve(object$cholprecondHessian, Ej)
            se[ind] <- object$diagHessian[ind] * sqrt(v[cbind(ind, seq_along(ind))])
            if (progress) setTxtProgressBar(pb, j)
          }
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
            mat <- .solve_pchol(object$cholprecondHessian, z)
            mat <- object$beta + object$diagHessian * mat
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
          J <- object$diagHessian * matrix(J, ncol = object$np)
          if (loop) {
          # for (j in 1:object$n) {
          #   ind <- ind0 + j
          #   ej <- Matrix::sparseVector(i = ind, x = J[j, ], length = nv)
          #   temp <- Matrix::solve(object$cholprecondHessian, ej)[ind]
          #   sek[j] <- sqrt(sum(ej[ind] * temp))
          #   if (progress) setTxtProgressBar(pb, j)
          # }
            spl <- split(1:object$n, c(0:(object$n - 1)) %/% chunksize)
            for (j in 1:length(spl)) {
              ind <- spl[[j]]
              Jj <- as.vector(t(J[ind, ]))
              indr <- as.integer(outer(ind0, ind, FUN = '+'))
              indc <- rep(seq_along(ind), each = object$np)
              Ej <- Matrix::sparseMatrix(i = indr, j = indc, x = Jj, dims = c(nv, length(ind)))
              temp <- Matrix::solve(object$cholprecondHessian, Ej)[indr, ]
              sek[ind] <- sqrt(Matrix::colSums(Ej[indr, ] * temp))
              if (progress) setTxtProgressBar(pb, j)
            }
          } else {
            indr <- as.integer(outer(ind0, 1:object$n, FUN = '+'))
            indc <- rep(1:object$n, each = object$np)
            EJ <- Matrix::sparseMatrix(i = indr, j = indc, x = as.vector(t(J)), dims = c(nv, object$n))
            sek <- sqrt(Matrix::colSums(EJ * Matrix::solve(object$cholprecondHessian, EJ)))
          }
        }
        se[[k]] <- matrix(sek, object$nx, object$ny)
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

.fitted_values <- function(pars, likdata) {
  pl <- split(pars, likdata$psplit)
  t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
}
