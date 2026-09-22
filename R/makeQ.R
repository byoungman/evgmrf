.inits_model <- function(model, order = 1, val, bymfns, hyper) {
  if (any(order > 1) && model != 'icar')
    stop('order > 1 currently only possible for ICAR model.')
  if (is.null(hyper$kappa) || hyper$kappa < 0)
    hyper$kappa = 2.7
  if (model == 'car') {
    if (is.null(hyper$rho) || hyper$rho < 0)
      hyper$rho <- .9
  }
  if (model == 'bym') {
    if (is.null(hyper$epsilon) || hyper$epsilon < 0) {
      hyper$epsilon <- .1 * hyper$kappa
      hyper$kappa <- .9 * hyper$kappa
    }
  }
  if (model == 'bym2') {
    if (is.null(hyper$rho) || hyper$rho < 0) {
      hyper$rho <- .9
    }
  }
  if (length(order) > 1 && any(order > 1))
    hyper$nu <- .5^(order[order > 1] - 1)
  if (!missing(val))
    hyper$kappa <- val
  hyper
}

.null_basis_1d <- function(m, q) {
  x <- seq_len(m)
  V <- outer(x, 0:(q - 1), `^`)
  qr.Q(qr(V))
}

.null_basis_Qn <- function(nx, ny, q_min) {
  Ux <- .null_basis_1d(nx, q_min)
  Uy <- .null_basis_1d(ny, q_min)
  as.matrix(kronecker(Uy, Ux))
}

.sparseinv_sQ <- function(Q, A, eps = NULL) {
  
  n <- nrow(Q)
  Q <- as(Q, "CsparseMatrix")
  Q <- Matrix::forceSymmetric(Q)
  
  A <- as(A, "CsparseMatrix")
  k <- nrow(A)
  
  # NOTE: do NOT perturb with c * crossprod(A) -- for a dense constraint
  # matrix A (e.g. the all-ones sum-to-zero direction) that outer product is
  # a fully dense n x n matrix, so Q + c*crossprod(A) destroys Q's sparsity
  # entirely and CHOLMOD runs out of memory trying to factorize it at large
  # n. A diagonal ridge preserves Q's sparsity pattern exactly; the
  # constraint correction is still applied exactly afterward via Woodbury.
  if (is.null(eps)) {
    d <- Matrix::diag(Q)
    eps <- 1e-8 * mean(d[d > 0])
  }
  Qeps <- as(Q + Matrix::Diagonal(n, eps), "CsparseMatrix")
  
  cp   <- sparseinv::cholPermute(Qeps)
  Sinv <- sparseinv::Takahashi_Davis(Qeps, cholQp = cp$Qpermchol, P = cp$P)
  d_star <- Matrix::diag(Sinv)
  
  ch <- Matrix::Cholesky(Qeps, perm = TRUE, LDL = FALSE)
  W  <- as.matrix(Matrix::solve(ch, Matrix::t(as.matrix(A)), system = "A"))   # n x k
  V  <- as.matrix(A %*% W)                                    # k x k
  
  # NOTE: colSums (not a bare transpose) -- the previous version only gave
  # the right answer for k = 1 (order-1 penalties); for k > 1 (order >= 2,
  # where q_min^2 >= 4) the transpose silently returned the wrong sQ.
  Y <- forwardsolve(t(chol(V)), t(W))        # k x n
  correction <- Matrix::colSums(Y^2)         # length n
  d_constr <- d_star - correction
  
  exp(mean(log(d_constr)))
}

.makeQ_data <- function(nx, ny, model, order, n_null, W = NULL, bymfns, hyper) {
  if (is.null(W)) {
    n <- nx * ny
    given_W <- FALSE
  } else {
    n <- try(nrow(W), silent = TRUE)
    n <- nrow(W)
    given_W <- TRUE
    if (!(inherits(W, "dgCMatrix") || inherits(W, "dsCMatrix")))
      stop('W must inherit class "dgCMatrix" or "dsCMatrix"')
    if (!all(W@x %in% c(0, 1)))
      stop('All elements of W must be in {0, 1}.')
    if (!Matrix::isSymmetric(W))
      stop('W must be symmetric')
    if (!all(Matrix::diag(W) == 0))
      stop('All diagonal elements of W must be zero.')
  }
  n_pars <- length(n_null)
  id <- seq_along(n_null)
  mods <- ords <- nus <- rep(NA, length(id))
  mods[id] <- model
  ords[id] <- order
  target <- list(icar = 1, 
                 car = c(1, -4), 
                 bym = c(1, -4), 
                 bym2 = c(1, -4),
                 bym3 = c(1))[mods]
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
  for (i in seq_along(hyper)) {
    if (is.null(hyper[[i]]$kappa))
      hyper[[i]]$kappa <- -1
    if (length(order[[i]]) > 1) {
      if (is.null(hyper[[i]]$nu)) {
        nu <- numeric(0)
        if (2 %in% order[[i]])
          nu <- c(nu, -1)
        if (3 %in% order[[i]])
          nu <- c(nu, -1)
        hyper[[i]]$nu <- nu
      }
    }
    if (mods[[i]] %in% c('car', 'bym2')) {
      if (is.null(hyper[[i]]$rho))
        hyper[[i]]$rho <- -1
    }
    if (mods[[i]] %in% 'bym') {
      if (is.null(hyper[[i]]$epsilon))
        hyper[[i]]$epsilon <- -1
    }
    if (mods[[i]] %in% 'bym3') {
      args <- names(formals(bymfns[[i]]))[-1]
      for (j in seq_along(args)) {
        if (is.null(hyper[[i]][[args[j]]]))
          hyper[[i]][[args[j]]] <- -1
      }
    }
  }
  hyper[is.na(model)] <- NULL
  hyper_swap <- lapply(seq_along(hyper), 
                       function(i) if (any(unlist(hyper[[i]]) == -1)) 
                         cbind(i, which(unlist(hyper[[i]]) == -1)))
  hyper_null <- lapply(hyper_swap, is.null)
  if (all(unlist(hyper_null))) {
    hyper_swap <- cbind(numeric(0), numeric(0))
  } else {
    hyper_swap <- do.call(rbind, hyper_swap)
  }
  Ql <- list()
  if (!given_W) {
    for (i in 1:max(unlist(order))) {
      Ql[[i]] <- .makeQ_order(nx, ny, i)
    }
  }
  target <- target[reps > 0]
  need2scale <- mods %in% c('bym2', 'bym3')
  lsQs <- numeric(length(order))
  if (any(need2scale)) {
    order2scale <- sapply(ords[need2scale], min)
    A2calculate <- unique(order2scale)
    A_null <- lapply(A2calculate, .null_basis_Qn, nx = nx, ny = ny)
    sinv <- numeric(3)
    sinv[A2calculate] <- sapply(seq_along(A2calculate), function(i) .sparseinv_sQ(Ql[[i]], Matrix::t(A_null[[i]])))
    lsQs[order2scale] <- log(sinv[order2scale])
  }
  list(n = n, n_null = n_null, W = W,
       ord = ords, mod = mods, nu = nus, spl = NULL, np = length(id), 
       target = target, Ql = Ql, nx = nx, ny = ny, hyper = hyper,
       hyper_swap = hyper_swap, lsQs = lsQs)
}

.hyper2pars <- function(hyper, swap) {
  out <- numeric(nrow(swap))
  if (nrow(swap) > 0) {
    for (i in 1:nrow(swap)) {
      ind1 <- swap[i, 1]
      ind2 <- rownames(swap)[i]
      out[i] <- unlist(hyper[[ind1]])[ind2]
      if (ind2 %in% c('kappa', 'epsilon')) {
        out[i] <- log(out[i])
      } else {
        if (substr(ind2, 1, 5) %in% c('rho', 'nu'))
          out[i] <- qnorm(out[i])
      }
    }
  }
  out
}

.pars2hyper <- function(pars, hyper, swap) {
  rownames(swap)[substr(rownames(swap), 1, 2) == 'nu'] <- 'nu'
  names(pars) <- rownames(swap)
  pars_list <- list()
  i_swap <- unique(swap[, 1])
  for (i in i_swap) {
    pars_list[[i]] <- pars[swap[, 1] == i]
    temp <- split(pars_list[[i]], names(pars_list[[i]]))
    for (j in 1:length(temp)) {
      if (names(temp)[j] %in% c('kappa', 'epsilon'))
        temp[[j]] <- exp(temp[[j]])
      if (names(temp)[j] %in% c('rho', 'nu'))
        temp[[j]] <- pnorm(temp[[j]])
    }
    hyper[[i]][names(temp)] <- temp
  }
  hyper
}

.mQ <- function(hyper, Qd, alpha.tol = 1e-6) {
  mods <- Qd$mod
  ords <- Qd$ord
  nus <- Qd$nu
  lsQs <- Qd$lsQs
  Ql <- list()
  for (i in 1:Qd$np) {
    Ql[[i]] <- .makeQ_any(hyper[[i]], Qd, mods[i], ords[[i]], lsQs[i], Qd$n_null[i], Qd$R[[i]], alpha.tol)
  }
  Q <- Matrix::.bdiag(Ql)
  attr(Q, 'logdet') <- sum(sapply(Ql, attr, 'logdet'))
  attr(Q, 'test') <- unique(unlist(sapply(Ql, attr, 'test')))
  attr(Q, 'pars') <- lapply(Ql, attr, "pars")
  Q
}

.makeQ_any <- function(pars, Qd, model, order, lsQ, n_null, R, rho.tol = 1e-6) {
  logdet <- 0
  nu <- c(1, NA, NA)
  if (!is.na(model)) {
    kappa <- pars$kappa
    lkappa <- log(kappa)
    if (model == 'car') {
      rho.tol <- rho.tol + (1 - rho.tol) * pars$rho
      Q_rank <- Qd$n
    } else {
      Q_rank <- Qd$n - max(order)
    }
    rho <- 1 - rho.tol
    if (any(order > 1) && length(order) > 1) {
      nu[order[order > 1]] <- pars$nu
    }
    Q <- .makeQ(Qd, rho, nu, order)
    if (model %in% c('icar', 'car'))
      Q <- Q + Matrix::Diagonal(n = Qd$n, x = 1 - rho)
    if (model == 'icar') {
      logdet <- Q_rank * lkappa
    } 
    if (model == 'car') {
      logdet <- .ldchol(Q)
      logdet <- logdet + Q_rank * lkappa
    }
    Q <- kappa * Q
    if (model == 'bym') {
      Q <- Q + Matrix::Diagonal(n = Qd$n, x = pars$epsilon)
      logdet <- .ldchol(Q)
    }
    if (model == 'bym2') {
      rho <- rho.tol + (1 - rho.tol) * pars$rho
      lrho <- log(rho)
      l1mrho <- log(1 - rho)
      logdet <- Q_rank * (lkappa - lrho + lsQ) + Qd$n * (lkappa - l1mrho)
      Q <- Q * exp(-lrho)
      Q <- list(Q, Matrix::Diagonal(n = Qd$n, x = exp(lkappa - l1mrho)))
      Q <- Matrix::bdiag(Q)
    }
    if (model == 'bym3') {
      logdet <- Q_rank * lkappa
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
  Q
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
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n - 2, n))
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

.makeQ <- function(Qd, rho, nu, order) {
  if (!is.null(Qd$W)) {
    D <- Matrix::Diagonal(x = pmax(Matrix::rowSums(Qd$W), 1))
    Q <- D - rho * Qd$W
  } else {
    Q <- 0 * Qd$Ql[[1]]
    for (i in seq_along(order))
      Q <- Q + nu[i] * Qd$Ql[[order[i]]]
  }
  return(Q)
}