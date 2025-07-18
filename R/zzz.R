.inits_model <- function(model, order = 1, alpha = NA) {
  if (is.na(model)) {
    out <- numeric(0)
  } else {
    if (model == 'icar') {
      out <- c(lambda = 1)
    } else {
      if (model == 'car') {
        out <- c(lambda = 1, rho = -4)
      } else {
        if (model == 'bym') {
          out <- c(lambda = 1, rho = -4)
        }
      }
    }
  }
  if (order > 1 & is.na(alpha))
    out <- c(out, alpha = 0)
  out        
}

.makeQ_data <- function(nx, ny, model, order, alpha, n_null, W = NULL) {
  if (is.null(W)) {
    n <- nx * ny
    W <- spatstat.sparse::gridadjacencymatrix(c(nx, ny), diagonal = FALSE)
  } else {
    n <- nrow(W)
  }
  # d <- Matrix::colSums(W)
  # ldetD <- sum(log(d))
  # D <- Matrix::Diagonal(n, d)
  # B <- Matrix::Diagonal(n, 1 / Matrix::colSums(W)) %*% W
  id <- seq_along(n_null)
  mods <- ords <- alphs <- rep(NA, length(id))
  mods[id] <- model
  ords[id] <- order
  alphs[id] <- alpha
  alphs[ords == 1] <- NA
  target <- list(icar = 1, car = c(1, -4), bym = c(1, -4))[mods]
  if (any(ords > 1))
    target[ords > 1] <- lapply(target[ords > 1], function(x) c(x, 0))
  reps <- sapply(mods, .model2n) + as.integer(ords > 1) - as.integer(is.finite(alphs))
  hyper_spl <- rep(id, reps)
  Ql <- list()
  # for (i in 1:3) {
  #   if (i %in% order)
  #     Ql[[i]] <- .make_Q(nx, ny, 1)
  # }
  # list(D = D, B = B, W = W, n = n, ldetD = ldetD, n_null = n_null, 
  #      ord = ords, mod = mods, spl = rho_spl, np = length(id), 
  #      target = target, Ql = Ql, nx = nx, ny = ny)
  # W <- .Wlist(W, max(ords))
  list(n = n, n_null = n_null, W = W,
       ord = ords, mod = mods, alph = alphs, spl = hyper_spl, np = length(id), 
       target = target, Ql = Ql, nx = nx, ny = ny)
}

.Wlist <- function(W, pow) {
  out <- list(W)
  if (pow > 1) {
    for (i in 2:pow) {
      out[[i]] <- out[[i - 1]] %*% W
    }
  }
  out
}

.matpow <- function(W0, pow) {
  W <- W0
  if (pow > 1) {
    for (i in 2:pow) {
      W <- W %*% W0
    }
  }
  W
}

.makeQ <- function(Wlist, rho, alpha) {
  alpha_seq <- alpha^(seq_along(Wlist) - 1)
  Wlist <- mapply('*', Wlist, alpha_seq)
  W <- Reduce('+', Wlist)
  D <- Matrix::Diagonal(x = Matrix::rowSums(W))
  D - rho * W
}

.makeQ <- function(W0, rho, alpha, order) {
  W <- W0
  if (order > 1) {
    for (i in 2:order)
      W <- W + (alpha^(i - 1)) * .matpow(W0, i)
  }
  # alpha_seq <- alpha^(seq_along(Wlist) - 1)
  # Wlist <- mapply('*', Wlist, alpha_seq)
  # W <- Reduce('+', Wlist)
  D <- Matrix::Diagonal(x = Matrix::rowSums(W))
  D - rho * W
}

.makeQ_any <- function(pars, Qd, model, order, alpha, n_null, R, rho.tol = 1e-6) {
  if (!is.na(model)) {
    lambda <- exp(pars['lambda'])
    if (model == 'car') {
      rho.tol <- rho.tol + (1 - rho.tol) * pnorm(pars['rho'])
    } 
    rho <- 1 - rho.tol
    if (order > 1 & is.na(alpha)) {
      alpha <- pnorm(pars['alpha'])
    } else {
      if (is.na(alpha))
        alpha <- 1
    }
    Q <- .makeQ(Qd$W, rho, alpha, order)
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
  } else {
    Q <- Matrix::Diagonal(0)
  }
  if (n_null > 0) {
    Q <- Matrix::.bdiag(list(Q, Matrix::Diagonal(n_null, numeric(n_null))))
    Q <- as(Q, 'CsparseMatrix')
  }
  attr(Q, 'logdet') <- logdet
  Q
}

.mQ <- function(rho, Qd, alpha.tol = 1e-6) {
  mods <- Qd$mod
  ords <- Qd$ord
  alphs <- Qd$alph
  hyper_spl <- .split2(rho, Qd$spl)
  Ql <- list()
  for (i in 1:Qd$np) {
    Ql[[i]] <- .makeQ_any(hyper_spl[[i]], Qd, mods[i], ords[i], alphs[[i]], Qd$n_null[i], Qd$R[[i]], alpha.tol)
  }
  Q <- Matrix::.bdiag(Ql)
  attr(Q, 'logdet') <- sum(sapply(Ql, attr, 'logdet'))
  Q
}
