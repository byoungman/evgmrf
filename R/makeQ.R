.makeQ_data <- function(nx, ny, model, order, alpha, n_null, W = NULL, bymfns) {
  if (is.null(W)) {
    n <- nx * ny
    given_W <- FALSE
  } else {
    n <- nrow(W)
    given_W <- TRUE
  }
  id <- seq_along(n_null)
  mods <- ords <- alphs <- rep(NA, length(id))
  mods[id] <- model
  ords[id] <- order
  alphs[id] <- alpha
  target <- list(icar = 1, 
                 car = c(1, -4), 
                 bym = c(1, -4), 
                 bym2 = c(1, -4),
                 bym3 = c(1),
                 bym4 = c(1))[mods]
  app <- lapply(ords, function(x) 0 * x[-1])
  target <- lapply(seq_along(target), function(i) c(target[[i]], app[[i]]))
  reps <- sapply(target, length)
  if (!is.null(bymfns)) {
    for (i in 1:length(bymfns)) {
      fnsi <- bymfns[[i]]
      if (!is.null(fnsi)) {
        nargs <- length(unlist(formals(fnsi))) - 1
        reps[i] <- reps[i] + nargs
      }
    }
  }
  hyper_spl <- rep(id, reps)
  Ql <- list()
  if (!given_W) {
    for (i in 1:max(unlist(order))) {
      Ql[[i]] <- .makeQ_order(nx, ny, i)
    }
  }
  target <- target[reps > 0]
  list(n = n, n_null = n_null, W = W,
       ord = ords, mod = mods, alph = alphs, spl = hyper_spl, np = length(id), 
       target = target, Ql = Ql, nx = nx, ny = ny)
}

.mQ <- function(rho, Qd, alpha.tol = 1e-6) {
  mods <- Qd$mod
  ords <- Qd$ord
  alphs <- Qd$alph
  hyper_spl <- lapply(1:Qd$np, function(.) numeric(0))
  hyper_spl[!is.na(Qd$mod)] <- .split2(rho, Qd$spl)
  # hyper_spl <- .split2(rho, Qd$spl)
  Ql <- list()
  for (i in 1:Qd$np) {
    Ql[[i]] <- .makeQ_any(hyper_spl[[i]], Qd, mods[i], ords[[i]], alphs[[i]], Qd$n_null[i], Qd$R[[i]], alpha.tol)
  }
  Q <- Matrix::.bdiag(Ql)
  attr(Q, 'logdet') <- sum(sapply(Ql, attr, 'logdet'))
  attr(Q, 'test') <- unique(unlist(sapply(Ql, attr, 'test')))
  attr(Q, 'pars') <- lapply(Ql, attr, "pars")
  attr(Q, 'splpars') <- hyper_spl
  Q
}

.makeQ_any <- function(pars, Qd, model, order, alpha, n_null, R, rho.tol = 1e-6) {
  logdet <- 0
  alpha.id <- names(pars) == ''
  if (!is.na(model)) {
    lambda <- exp(pars['lambda'])
    if (model == 'car') {
      rho.tol <- rho.tol + (1 - rho.tol) * pnorm(pars['rho'])
    } 
    rho <- 1 - rho.tol
    if (any(order > 1) & is.na(alpha)) {
      alpha <- pnorm(pars[alpha.id])
    } else {
      if (is.na(alpha))
        alpha <- 1
    }
    Q <- .makeQ(Qd, rho, alpha, order)
    if (model != 'bym')
      Q <- Q + Matrix::Diagonal(n = Qd$n, x = 1 - rho)
    if (model == 'icar') {
      logdet <- (Qd$n - 1) * pars['lambda']
    } else {
      if (model == 'car') {
        logdet <- .ldchol(Q)#Matrix::determinant(Q)$modulus
        logdet <- logdet + Qd$n * pars['lambda']
      }
    }
    Q <- lambda * Q
    if (model == 'bym') {
      Q <- Q + Matrix::Diagonal(n = Qd$n, x = rho.tol + exp(pars['rho']))
      logdet <- .ldchol(Q)
    }
    if (model == 'bym2') {
      lkappa <- pars['lambda']
      kappa <- exp(lkappa)
      lrho <- pnorm(pars['rho'], log.p = TRUE)
      l1mrho <- pnorm(pars['rho'], log.p = TRUE, lower.tail = FALSE)
      rhoc <- pnorm(pars['rho'], lower.tail = FALSE)
      logdet <- (Qd$n - 1) * (lkappa - lrho) + Qd$n * (lkappa - l1mrho)
      Q <- Q * exp(-lrho)
      Q <- list(Q, Matrix::Diagonal(n = Qd$n, x = exp(lkappa - l1mrho)))
      Q <- Matrix::bdiag(Q)
    }
    if (model %in% paste('bym', 3:4, sep = '')) {
      p1 <- lambda
      logdet <- (Qd$n - 1) * log(p1)
      Q <- list(Q, Matrix::Diagonal(n = Qd$n, x = 0))
      Q <- Matrix::bdiag(Q)
    }
  } else {
    Q <- Matrix::Diagonal(0)
  }
  if (n_null > 0) {
    Q <- Matrix::.bdiag(list(Q, Matrix::Diagonal(n_null, numeric(n_null))))
    Q <- as(Q, 'CsparseMatrix')
  }
  attr(Q, 'logdet') <- logdet
  if (!is.na(model)) {
    if (model == 'bym4')
      attr(Q, 'pars') <- pars
  }
  Q
}

.model2n <- function(model) {
  if (is.na(model))
    return(0)
  if (model == 'icar') {
    return(1)
  } else {
    return(2)
  }
}

.make_D1 <- function(n) {
  if (n == 1) {
    out <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(1, 1))
  } else {
    i <- 1:(n - 1)
    j <- c(1:(n - 1), 2:n)
    x <- c(rep(-1, n - 1), rep(1, n - 1))
    out <- Matrix::sparseMatrix(i = rep(i, 2), j = j, x = x, dims = c(n - 1, n))
  }
  out
}

.make_D2 <- function(n) {
  i <- c(1:(n-2), 1:(n-2), 1:(n-2))
  j <- c(1:(n-2), 2:(n-1), 3:n)
  x <- c(rep(1, n-2), rep(-2, n-2), rep(1, n-2))
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n-2, n))
}

.make_D3 <- function(n) {
  i <- rep(1:(n - 3), each = 4)
  j <- c(1:(n - 3), 2:(n - 2), 3:(n - 1), 4:n)
  x <- rep(c(1, -3, 3, -1), n - 3)
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n - 3, n))
}

.makeQ_order <- function(nx, ny, order) {
  if (order == 1) {
    Dx <- .make_D1(nx)
    Dy <- .make_D1(ny)
  }
  if (order == 2) {
    Dx <- .make_D2(nx)
    Dy <- .make_D2(ny)
  }
  if (order == 3) {
    Dx <- .make_D3(nx)
    Dy <- .make_D3(ny)
  }
  Qx <- Matrix::kronecker(Matrix::Diagonal(ny), Matrix::crossprod(Dx))
  Qy <- Matrix::kronecker(Matrix::crossprod(Dy), Matrix::Diagonal(nx))
  as(Qx + Qy, "generalMatrix")
}

.makeQ <- function(Qd, rho, alpha, order) {
  if (!is.null(Qd$W)) {
    D <- Matrix::Diagonal(x = Matrix::rowSums(Qd$W))
    Q <- D - rho * Qd$W
  } else {
    Q <- 0 * Qd$Ql[[1]]
    alpha <- c(1, alpha)
    # if (any(order > 1)) {
      for (i in seq_along(order))
        Q <- Q + alpha[i] * Qd$Ql[[order[i]]]
  # }
  }
  return(Q)
}
