.inits_model <- function(model, order = 1, alpha = NA, val) {
  if (is.na(model)) {
    out <- numeric(0)
  } else {
    if (model == 'icar') {
      out <- c(lambda = 1)
    } else {
      if (model == 'car') {
        out <- c(lambda = 1, rho = -4)
      } else {
        if (substr(model, 1, 3) == 'bym') {
          if (model == 'bym3') {
            out <- c(lambda = val)
          } else {
            if (model == 'bym4') {
              out <- c(lambda = 1, s = log(.1))
            } else {
              out <- c(lambda = val, rho = -1.6)
              out <- c(lambda = -3, rho = val)
            }
          }
        }
      }
    }
  }
  if (order > 1 & is.na(alpha))
    out <- c(out, -(2:order - 2))
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
  target <- list(icar = 1, 
                 car = c(1, -4), 
                 bym = c(1, -4), 
                 bym2 = c(1, -4),
                 bym3 = c(1),
                 bym4 = c(1, -4))[mods]
  app <- lapply(ords, function(x) rep(0, x - 1))
  target <- lapply(seq_along(target), function(i) c(target[[i]], app[[i]]))
  # if (any(ords > 1))
  #   target[ords > 1] <- lapply(target[ords > 1], function(x) c(x, rep(0, ords[ords > 1] - 1)))
  # reps <- sapply(mods, .model2n) + as.integer(ords > 1) - as.integer(is.finite(alphs))
  reps <- sapply(target, length)
  hyper_spl <- rep(id, reps)
  Ql <- list()
  for (i in 1:max(order)) {
      Ql[[i]] <- .makeQ_order(nx, ny, i)
  }
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

# .makeQ <- function(Wlist, rho, alpha) {
#   alpha_seq <- alpha^(seq_along(Wlist) - 1)
#   Wlist <- mapply('*', Wlist, alpha_seq)
#   W <- Reduce('+', Wlist)
#   D <- Matrix::Diagonal(x = Matrix::rowSums(W))
#   D - rho * W
# }
# 
# .makeQ <- function(W0, rho, alpha, order) {
#   W <- W0
#   if (order > 1) {
#     for (i in 2:order)
#       W <- W + (alpha^(i - 1)) * .matpow(W0, i)
#   }
#   # alpha_seq <- alpha^(seq_along(Wlist) - 1)
#   # Wlist <- mapply('*', Wlist, alpha_seq)
#   # W <- Reduce('+', Wlist)
#   D <- Matrix::Diagonal(x = Matrix::rowSums(W))
#   D - rho * W
# }

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
  attr(Q, 'test') <- unique(unlist(sapply(Ql, attr, 'test')))
  attr(Q, 'pars') <- lapply(Ql, attr, "pars")
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
    if (order > 1 & is.na(alpha)) {
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
      tau <- lambda
      rho <- pnorm(pars['rho'])
      p1 <- tau / rho
      lp1 <- pars['lambda'] - log(rho)
      p2 <- tau / pnorm(pars['rho'], lower.tail = FALSE)
      lp2 <- pars['lambda'] - log(1 - rho)
      # logdet <- (Qd$n - 1) * lp1 + Qd$n * lp2
      p1 <- lambda
      p2 <- exp(pars['rho'])
      logdet <- (Qd$n - 1) * log(p1) + 0 * Qd$n * log(p2)
      # Q <- Q / rho
      Q <- list(Q, Matrix::Diagonal(n = Qd$n, x = 0 * p2))
      Q <- Matrix::bdiag(Q)
      attr(Q, 'test') <- p2
    }
    if (model %in% paste('bym', 3:4, sep = '')) {
      tau <- lambda
      p1 <- lambda
      logdet <- (Qd$n - 1) * log(p1)
      Q <- list(Q, Matrix::Diagonal(n = Qd$n, x = 0))
      Q <- Matrix::bdiag(Q)
      attr(Q, 'test') <- 0
    }
  } else {
    Q <- Matrix::Diagonal(0)
  }
  if (n_null > 0) {
    Q <- Matrix::.bdiag(list(Q, Matrix::Diagonal(n_null, numeric(n_null))))
    Q <- as(Q, 'CsparseMatrix')
  }
  attr(Q, 'logdet') <- logdet
  if (model == 'bym2')
    attr(Q, 'test') <- p2
  if (model == 'bym4')
    attr(Q, 'pars') <- pars
  Q
}

.makeQ <- function(Qd, rho, alpha, order) {
  W0 <- W <- Qd$W
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
  Q <- Qd$Ql[[1]]
#  Ql <- mapply('*', Qd$Ql, c(1, alpha))
#  return(Reduce('+', Ql))
#  browser()
  if (order > 1) {
    for (i in 2:order)
      Q <- Q + alpha[i - 1] * Qd$Ql[[i]]
  }
  return(Q)
  # alpha_seq <- alpha^(seq_along(Wlist) - 1)
  # Wlist <- mapply('*', Wlist, alpha_seq)
  # W <- Reduce('+', Wlist)
  D <- Matrix::Diagonal(x = Matrix::rowSums(W))
  D - rho * W
}
