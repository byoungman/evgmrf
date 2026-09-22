.args0 <- list(delta = .1, mult = 1, C = 1, tau = NULL, nper = 1)

.checks <- function(model, order) {
  if (!.check_multiple(order, model))
    stop('length(order) not multiple of length(model) or vice-versa.')
}

.model2formula <- function(model) {
  if (is.character(model)) {
    if (model %in% c('icar', 'bym2', 'bym3')) {
      return( ~ -1)
    } else {
      return( ~ 1)
    }
  } else {
    if (is.na(model)) {
      return( ~ 1)
    } else {
      stop('Unrecognized model')
    }
  }
}

.list2mat <- function(list) {
  rows <- lapply(list, seq_along)
  cols <- rep(seq_along(rows), sapply(rows, length))
  rows <- unlist(rows)
  n <- length(list)
  m <- max(rows)
  out <- matrix(NA, m, n)
  out[cbind(rows, cols)] <- unlist(list)
  out
}

.list2array <- function(list) {
  rows <- lapply(list, function(x) as.matrix(expand.grid(seq_len(nrow(x)), seq_len(ncol(x)))))
  for (i in seq_along(rows)) rows[[i]] <- cbind(rows[[i]], i)
  ind <- do.call(rbind, rows)
  out <- array(NA, apply(ind, 2, max))
  out[ind] <- unlist(list)
  out
}

.pend012 <- function(pars, fn, lst, deriv = 0, eps = 1e-4) {
  out <- list()
  lst$x <- pars
  f0 <- do.call(fn, lst)
  out[[1]] <- sum(f0)
  if (deriv == 0)
    return(out[[1]])
  ph <- pars + eps
  pl <- pars - eps
  lst$x <- ph
  fh <- do.call(fn, lst)
  lst$x <- pl
  fl <- do.call(fn, lst)
  out[[2]] <- .5 * (fh - fl) / eps
  out[[3]] <- (fh + fl - 2 * f0) / (eps^2)
  out
}

.d0_Q <- function(pars, likdata, likfns, Q, hyper, diag = FALSE) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  out0 <- likdata$mult * likfns$d0(as.matrix(pm), likdata)
  if (!is.null(likdata$bymfns)) {
    for (i in seq_along(pl)) {
      if (!is.null(likdata$bymfns[[i]])) {
        xi <- pl[[i]][likdata$id_bym2[[i]]]
        parsi <- hyper[[i]][names(formals(likdata$bymfns[[i]]))[-1]]
        out0 <- out0 + .pend012(xi, likdata$bymfns[[i]], parsi)
      }
    }
  }
  out <- out0 + .5 * crossprod(pars, Q %*% pars)[1, 1]
  if (!is.finite(out))
    out <- 1e20
  attr(out, 'unpenalized') <- out0
  out
}

.d12_Q <- function(pars, likdata, likfns, Q, hyper) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  gH <- likfns$d12(pm, likdata)
  out <- list(g = likdata$mult * as.vector(gH[[1]] %*% likdata$X) + as.vector(Q %*% pars))
  n <- length(likdata$z)
  p <- nrow(pm)
  r1 <- n * rep(0:(p - 1), p:1)
  c1 <- n * unlist(sapply(1:p, function(i) i:p - 1))
  n2 <- rep(1:n, each = sum(1:p))
  r2 <- r1 + n2
  c2 <- c1 + n2
  H <- gH$H
  H <- Matrix::sparseMatrix(r2, c2, x = as.vector(t(H)), symmetric = TRUE)
  H <- likdata$mult * crossprod(likdata$X, H %*% likdata$X)
  out$H <- H + Q
  gl <- Hl <- lapply(pl, function(x) 0 * x)
  if (!is.null(likdata$bymfns)) {
    for (i in seq_along(pl)) {
      if (!is.null(likdata$bymfns[[i]])) {
        parsi <- hyper[[i]][names(formals(likdata$bymfns[[i]]))[-1]]
        temp <- .pend012(pl[[i]][likdata$id_bym2[[i]]], fn = likdata$bymfns[[i]], parsi, deriv = 2)
        gl[[i]][likdata$id_bym2[[i]]] <- temp[[2]]
        Hl[[i]][likdata$id_bym2[[i]]] <- temp[[3]]
      }
    }
    out$g <- out$g + unlist(gl)
    H0 <- Matrix::Diagonal(n = length(out$g), x = unlist(Hl))
    out$H <- out$H + H0
  }
  out
}

.perturb <- function(A, b = NULL, tol = 1e-1, mult = 1e2) {
  d0 <- Matrix::diag(A)
  test <- suppressWarnings(try(.cholAb(A, b), silent = TRUE))
  while(inherits(test, "try-error")) {
    Matrix::diag(A) <- d0 + tol
    test <- suppressWarnings(try(.cholAb(A, b), silent = TRUE))
    tol <- mult * tol
    if (!is.finite(tol))
      stop("Can't perturb Hessian to be positive definite.")
  }
  attr(A, 'chol') <- test
  A
}

.perturb_super <- function(A, b = NULL, tol = 1e-1, mult = 1e2, super = TRUE) {
  A0 <- A
  d0 <- Matrix::diag(A)
  test <- suppressWarnings(try(Matrix::Cholesky(A, LDL = FALSE, super = super), silent = TRUE))
  while(inherits(test, "try-error")) {
    Matrix::diag(A) <- d0 + tol
    test <- suppressWarnings(try(Matrix::Cholesky(A, LDL = FALSE, super = super), silent = TRUE))
    tol <- mult * tol
    if (!is.finite(tol))
      stop("Can't perturb Hessian to be positive definite.")
  }
  attr(A, 'chol') <- test#.chol_logdet_solve(A, b)
  A
}

.Cholesky0 <- function(rho, Qd, ridge = 1) {
  Q <- .mQ(rho, Qd)
  Q <- as(Q, 'symmetricMatrix')
  n <- nrow(Q)
  Q <- Q + Matrix::Diagonal(n, rep(ridge, n))
  Matrix::Cholesky(Q, LDL = FALSE, super = TRUE)
}

.search_Q0 <- function(pars, likdata, likfns, Q, hyper) {
  gH <- .d12_Q(pars, likdata, likfns, Q, hyper)
  H <- gH$H
  D <- Matrix::Diagonal(nrow(H), 1 / sqrt(pmax(Matrix::diag(H), 1e-8)))
  H <- D %*% H %*% D
  b <- as.vector(D %*% gH$g)
  H <- .perturb(H, b, likdata$opts.perturb$tol, likdata$opts.perturb$mult)
  cholH <- attr(H, 'chol')
  stp <- D %*% cholH$z
  ldet <- cholH$logdet_A
  ldet <- ldet  - 2 * sum(log(Matrix::diag(D)))
  iD <- Matrix::Diagonal(nrow(H), sqrt(Matrix::diag(H)))
  H <- iD %*% H %*% iD
  attr(H, 'ldet') <- ldet
  attr(gH$g, 'ldet') <- ldet
  attr(stp, 'gradient') <- gH$g
  stp
}

.search_Q <- function(pars, likdata, likfns, Q, hyper, kept = NULL, diag = FALSE) {
  gH <- .d12_Q(pars, likdata, likfns, Q, hyper)
  H <- gH$H
  d <- pmax(Matrix::diag(H), 1e-8)
  if (!diag) {
    D <- Matrix::Diagonal(nrow(H), 1 / sqrt(d))
    H <- D %*% H %*% D
    b <- as.vector(D %*% gH$g)
    if (likdata$control$inner_optim == 'Cholesky') {
      H <- .perturb_super(H, likdata$chol0, likdata$control$perturb.tol, likdata$control$perturb.mult, likdata$control$super)
      cholH <- attr(H, 'chol')
      stp <- D %*% Matrix::solve(cholH, b)
      ldet <- as.vector(Matrix::determinant(cholH, sqrt = FALSE)$modulus)
    } else {
      H <- .perturb(H, b, likdata$control$perturb.tol, likdata$control$perturb.mult)
      cholH <- attr(H, 'chol')
      stp <- D %*% cholH$z
      ldet <- cholH$logdet_A
    }
    ldet <- ldet  - 2 * sum(log(Matrix::diag(D)))
  } else {
    stp <- gH$g / d
    ldet <- sum(log(Matrix::diag(H)))
  }
  attr(H, 'ldet') <- ldet
  attr(gH$g, 'ldet') <- ldet
  if (any(!is.finite(gH$g)))
    stop('Non-finite gradient')
  attr(stp, 'gradient') <- gH$g
  attr(stp, 'H0') <- gH$H
  attr(stp, 'precondHessian') <- H
  attr(stp, 'diagHessian') <- D
  if (!diag) {
    attr(stp, 'cholprecondHessian') <- cholH
    iD <- Matrix::Diagonal(nrow(H), sqrt(d))
    H <- iD %*% H %*% iD
    attr(stp, 'idiagHessian') <- iD
  }
  attr(stp, 'Hessian') <- H
  stp
}

## Shared functions

.split2 <- function(x, id) {
  out <- list()
  uid <- unique(id)
  for (i in uid) {
    out[[uid[i]]] <- x[id == uid[i]]
  }
  out
}

.check_multiple <- function(x, y) {
  ratio <- length(x) / length(y)
  (ratio == round(ratio)) | (1 / ratio == round(1 / ratio))
}
