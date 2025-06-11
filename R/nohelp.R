.checks <- function(model, order) {
  if (!.check_multiple(order, model))
    stop('length(order) not multiple of length(model) or vice-versa.')
}

.inits_model <- function(model) {
  if (is.na(model)) {
    out <- numeric(0)
  } else {
    if (model == 'icar') {
      out <- 1
    } else {
      if (model == 'car') {
        out <- c(1, -4)
      } else {
        if (model == 'bym') {
          out <- c(1, -4)
        }
      }
    }
  }
  out        
}

.control.evgmrf <- function() {
  list(eps = 5e-3, it0 = 20, step_size = 0.2, reml_eps = 5e-3, reml_direction = 'ad-hoc',
       reml_steptol = 1e-2, reml_stepmax = 1, reml_itlim = 1e2, inner_optim = 'chol',
       alpha.tol = 1e-6, grad_mult = 0, par_mult = 1, super = TRUE, update = FALSE, 
       openmp = TRUE, threads = 0, perturb.tol = 1e-2, perturb.mult = 5, 
       perturb.method = 'chol', perturb.tol.eigen = 1e-3, perturb.mult.eigen = 10)
}

## REML functions

.d0_Q <- function(pars, likdata, likfns, Q) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  out <- likdata$mult * likfns$d0(as.matrix(pm), likdata)
  out <- out + .5 * crossprod(pars, Q %*% pars)[1, 1]
  if (!is.finite(out))
    out <- 1e20
  out
}

.d1_Q <- function(pars, likdata, likfns, Q) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  out <- likfns$d12(pm, likdata)[[1]]
  out <- likdata$mult * as.vector(out %*% likdata$X)
  out <- out + as.vector(Q %*% pars)
  print(c(min(out), quantile(out, c(.01, seq(.1, .9, by = .1), .99)), mean(out), mean(abs(out)), median(abs(out)), sqrt(sum(out^2))))
  out
}

# blockChol3 <- function(x) {
#   l <- 0 * x
#   l[, 1] <- sqrt(x[, 1])
#   l[, 2] <- x[, 2] / l[, 1]
#   l[, 3] <- x[, 3] / l[, 1]
#   l[, 4] <- sqrt(x[, 4] - l[, 2] * l[, 2])
#   l[, 5] <- (x[, 5] - l[, 3] * l[, 2]) / l[, 4]
#   l[, 6] <- sqrt(x[, 6] - l[, 3] * l[, 3] - l[, 5] * l[, 5])
#   l
# }

.d2_Q <- function(pars, likdata, likfns, Q) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  out <- likfns$d12(as.matrix(pm), likdata)
  H0 <- out$H
  n <- length(likdata$z)
  p <- nrow(pm)
  r1 <- n * rep(0:(p - 1), p:1)
  c1 <- n * unlist(sapply(1:p, function(i) i:p - 1))
  n2 <- rep(1:n, each = sum(1:p))
  r2 <- r1 + n2
  c2 <- c1 + n2
  out <- Matrix::sparseMatrix(r2, c2, x = as.vector(t(H0)), symmetric = TRUE)
  out <- likdata$mult * crossprod(likdata$X, out %*% likdata$X)
  out <- out + Q
  out
}

.d2_Q_diag <- function(pars, likdata, likfns, Q) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  out <- likfns$d12(as.matrix(pm), likdata)
  H0 <- out$H
  n <- length(likdata$z)
  p <- nrow(pm)
  r1 <- n * rep(0:(p - 1), p:1)
  c1 <- n * unlist(sapply(1:p, function(i) i:p - 1))
  n2 <- rep(1:n, each = sum(1:p))
  r2 <- r1 + n2
  c2 <- c1 + n2
  out <- Matrix::sparseMatrix(r2, c2, x = as.vector(t(H0)), symmetric = TRUE)
  out <- likdata$mult * crossprod(likdata$X, out %*% likdata$X)
  out <- out + Q
  Matrix::diag(out)
}


.d12_Q <- function(pars, likdata, likfns, Q) {
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
  out
}

# # A should be a symmetric sparse matrix, e.g., class "dgCMatrix"
# smallest_eigen <- function(A, which = "SA") {
#   RSpectra::eigs_sym(A, k = 1, which = which)$values
# }
# 
# eig1_fn <- function(A) {
#   # RSpectra::eigs_sym(A, 1, which = 'SA', opts = list(tol = 1e-10, ncv = 5))$values
#   # RSpectra::eigs_sym(H, 1, which = 'SA', opts = list(tol = 1e-10, ncv = 5, retvec = FALSE, maxitr = 1e4))$values
#   RSpectra::eigs_sym(A, 1, which = 'SA', opts = list(tol = 1e-6, ncv = 20, retvec = FALSE, maxitr = 1e4))$values
# }
# 
# eign_fn <- function(A) {
#   # RSpectra::eigs_sym(A, 1, which = 'SA', opts = list(tol = 1e-10, ncv = 5))$values
#   # RSpectra::eigs_sym(H, 1, which = 'SA', opts = list(tol = 1e-10, ncv = 5, retvec = FALSE, maxitr = 1e4))$values
#   RSpectra::eigs_sym(A, 1, which = 'LA', opts = list(tol = 1e-2, ncv = 20, retvec = FALSE, maxitr = 1e4))$values
# }

.perturb_eigen <- function(A, b, tol = 1e-6, mult = 1) {
  d0 <- Matrix::diag(A)
  ev1 <- .eig1(A)
  if (ev1 - tol < 0) {
    cond <- TRUE
    while (cond) {
      eps <- tol + abs(ev1)
      A <- A + Matrix::Diagonal(n = nrow(A), x = eps)
      test <- suppressWarnings(try(.cholAb(A, b), silent = TRUE))
      cond <- inherits(test, "try-error")
      tol <- mult * tol
    }
    attr(A, 'chol') <- test
  } else {
    attr(A, 'chol') <- .cholAb(A, b)
  }
  A
}

# perturb4 <- function(A, tol = 1e-2) {
#   d0 <- Matrix::diag(A)
#   eig1 <- eig1_fn(A)
#   if (eig1 - tol < 0) {
#     eps <- tol + abs(eig1)
#     Matrix::diag(A) <- d0 + eps
#   }
#   A
# }

.perturb <- function(A, b = NULL, tol = 1e-1, mult = 1e2) {
  d0 <- Matrix::diag(A)
  test <- suppressWarnings(try(.cholAb(A, b), silent = TRUE))
  # print(paste('tol =', 0))
  while(inherits(test, "try-error")) {
    # print(tol)
    Matrix::diag(A) <- d0 + tol
    test <- suppressWarnings(try(.cholAb(A, b), silent = TRUE))
    tol <- mult * tol
    if (!is.finite(tol))
      stop("Can't perturb Hessian to be positive definite.")
  }
  attr(A, 'chol') <- test#.chol_logdet_solve(A, b)
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

.perturb_update <- function(A, chol0, tol = 1e-1, mult = 1e2) {
  d0 <- Matrix::diag(A)
  test <- suppressWarnings(try(Matrix::update(chol0, A), silent = TRUE))
  while(inherits(test, "try-error")) {
    # D <- Matrix::Diagonal(d0 + tol, n = length(d0))
    # test <- suppressWarnings(Matrix::updown('+', D, chol0))
    Matrix::diag(A) <- d0 + tol
    test <- suppressWarnings(try(Matrix::update(chol0, A), silent = TRUE))
    tol <- mult * tol
    if (!is.finite(tol))
      stop("Can't perturb Hessian to be positive definite.")
  }
  attr(A, 'chol') <- test
  A
}

# perturb2 <- function(A, b = NULL, L0) {
#   d0 <- Matrix::diag(A)
#   eps <- 1e-8
#   L <- suppressWarnings(try(Matrix::update(L0, A), silent = TRUE))
#   while(inherits(L, "try-error")) {
#     print(eps)
#     Matrix::diag(A) <- d0 + eps
#     L <- suppressWarnings(try(Matrix::update(L0, A), silent = TRUE))
#     eps <- 1e2 * eps
#     if (!is.finite(eps))
#       stop("Can't perturb Hessian to be positive definite.")
#   }
#   attr(A, 'chol') <- .chol_logdet_solve(A, b)
#   A
# }

.Cholesky0 <- function(pars, likdata, likfns, Q, super = FALSE) {
  gH <- .d12_Q(pars, likdata, likfns, Q)
  H <- gH$H
  D <- Matrix::Diagonal(nrow(H), 1 / sqrt(pmax(Matrix::diag(H), 1e-8)))
  H <- D %*% H %*% D
  b <- as.vector(D %*% gH$g)
  d0 <- Matrix::diag(H)
  eps <- 1e-8
  L <- suppressWarnings(try(Matrix::Cholesky(H, LDL = FALSE, super = super), silent = TRUE))
  while(inherits(L, "try-error")) {
    print(eps)
    Matrix::diag(H) <- d0 + eps
    L <- suppressWarnings(try(Matrix::Cholesky(H, LDL = FALSE, super = super), silent = TRUE))
    eps <- 1e2 * eps
  }
  L
}

.Cholesky0 <- function(rho, Qd, ridge = 1) {
  Q <- .mQ(rho, Qd)
  Q <- as(Q, 'symmetricMatrix')
  n <- nrow(Q)
  Q <- Q + Matrix::Diagonal(n, rep(ridge, n))
  Matrix::Cholesky(Q, LDL = FALSE, super = TRUE)
}

.search_Q0 <- function(pars, likdata, likfns, Q) {
  gH <- .d12_Q(pars, likdata, likfns, Q)
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

.search_Q <- function(pars, likdata, likfns, Q) {
  gH <- .d12_Q(pars, likdata, likfns, Q)
  H <- gH$H
  # if (.ld$precondition) {
  d <- pmax(Matrix::diag(H), 1e-8)
    D <- Matrix::Diagonal(nrow(H), 1 / sqrt(d))
    H <- D %*% H %*% D
    b <- as.vector(D %*% gH$g)
    if (likdata$control$perturb.method == 'eigen') {
      H <- .perturb_eigen(H, b, likdata$control$perturb.tol.eigen, likdata$control$perturb.mult.eigen)
      cholH <- attr(H, 'chol')
      stp <- D %*% cholH$z
      ldet <- cholH$logdet_A
    } else {
      if (likdata$control$inner_optim == 'Cholesky') {
        if (likdata$control$update) {
          H <- .perturb_update(H, likdata$chol0, likdata$control$perturb.tol, likdata$control$perturb.mult)
        } else {
          H <- .perturb_super(H, likdata$chol0, likdata$control$perturb.tol, likdata$control$perturb.mult, likdata$control$super)
        }
        cholH <- attr(H, 'chol')
        stp <- D %*% Matrix::solve(cholH, b)
        ldet <- as.vector(Matrix::determinant(cholH, sqrt = FALSE)$modulus)
      } else {
      H <- .perturb(H, b, likdata$control$perturb.tol, likdata$control$perturb.mult)
      cholH <- attr(H, 'chol')
      stp <- D %*% cholH$z
      ldet <- cholH$logdet_A
    }
  }
  ldet <- ldet  - 2 * sum(log(Matrix::diag(D)))
  attr(H, 'ldet') <- ldet
  attr(gH$g, 'ldet') <- ldet
  attr(stp, 'gradient') <- gH$g
  attr(stp, 'H0') <- gH$H
  attr(stp, 'cholprecondHessian') <- cholH
  attr(stp, 'precondHessian') <- H
  iD <- Matrix::Diagonal(nrow(H), sqrt(d))
  H <- iD %*% H %*% iD
  attr(stp, 'diagHessian') <- D
  attr(stp, 'idiagHessian') <- iD
  attr(stp, 'Hessian') <- H
  stp
}

.inits <- function(pars, likdata, likfns, Qd, makeQ) {
  beta <- attr(pars, 'beta')
  Q <- makeQ(pars, Qd)
  fit <- .newton(beta, .d0_Q, .search_Q, likdata = likdata, likfns = likfns, Q = Q, stepmax = 3)
  out <- fit$par
  attr(out, 'gradconv') <- fit$gradconv
  out
}
# 
# .reml0_test <- function(i) {
#   likdata2 <- likdata
#   likdata2$z <- likdata2$z[1:i]
#   likdata2$u <- likdata2$u[1:i]
#   beta2 <- matrix(beta, ncol = 3)[1:5, ]
#   .d0_Q(

.beta0 <- function(pars, Qd, likdata, likfns, makeQ, it0) {
  beta <- attr(pars, 'beta')
  Q <- makeQ(pars, Qd)
  .newton(beta, .d0_Q, .search_Q0, likdata = likdata, likfns = likfns, Q = Q, itlim = it0)
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


# .quick_gev <- function(y) {
#   psi0 <- sqrt(6 * var(unlist(y), na.rm = TRUE)) / pi
#   mu0 <- mean(unlist(y), na.rm = TRUE) - 0.57722 * psi0
#   inits <- c(mu0, log(psi0), .1)
#   nlminb(inits, gev0, gev1, gev2, yv = y)$par
# }

