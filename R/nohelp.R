.args0 <- list(delta = .1, mult = 1, C = 1, tau = NULL, nper = 1)

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
        if (model == 'bym4') {
          out <- 1
        } else {
          out <- c(1, -4)
        }
      }
    }
  }
  out        
}

.control.evgmrf <- function() {
  out <- list(eps = 5e-3, it0 = 20, step_size = 0.2, reml_eps = 5e-3, reml_direction = 'ad-hoc',
       reml_steptol = 5e-3, reml_stepmax = 1, reml_itlim = 1e2, inner_optim = 'chol',
       alpha.tol = 1e-6, grad_mult = 0, par_mult = .1, update = FALSE, 
       openmp = FALSE, threads = 0, perturb.tol = 1e-2, perturb.mult = 5, 
       perturb.method = 'chol', perturb.tol.eigen = 1e-3, perturb.mult.eigen = 10,
       cv_eps = 5e-2, cv_gradtol = .05, cv_steptol = 1e-4, super = FALSE,
       sandwich = FALSE)
  if (out$inner_optim == 'Cholesky' | out$perturb.method == 'Cholesky')
    out$super <- TRUE
  out
}

## REML functions

.dst0 <- function(x, omega = 1, alpha = 0, nu = 10, log = FALSE) {
  delta <- alpha / sqrt(1 + alpha^2)
  b_nu <- sqrt(nu / pi) * gamma((nu - 1) / 2) / gamma(nu / 2)
  xi <- -omega * b_nu * delta  # shift so mean = 0
  sn::dst(x, xi = xi, omega = omega, alpha = alpha, nu = nu, log = log)
}

.dmixexp0 <- function(x, p = 0.5, rate = 1, log = FALSE) {
  if (p < 0 || p > 1) stop("p must be between 0 and 1")
  if (rate <= 0) stop("rate must be positive")
  
  dens <- numeric(length(x))
  dens[x == 0] <- 1 - p
  dens[x > 0] <- p * rate * exp(-rate * x[x > 0])
  dens[x < 0] <- 0
  
  if (log) dens <- log(dens)
  
  return(dens)
}

.dmixexpnorm <- function(x, p = 0.5, rate = 1, sd = 1, log = FALSE) {
  if (p < 0 || p > 1) stop("p must be between 0 and 1")
  if (rate <= 0) stop("rate must be positive")
  if (sd <= 0) stop("sd must be positive")
  
  # Exponential pdf (only for x >= 0)
  f_exp <- ifelse(x >= 0, rate * exp(-rate * x), 0)
  
  # Normal pdf (symmetric, defined everywhere)
  f_norm <- dnorm(x, mean = 0, sd = sd)
  
  # Mixture
  dens <- p * f_exp + (1 - p) * f_norm
  
  if (log) dens <- log(dens)
  return(dens)
}

.temp <- function(x, s = .1) {
  ifelse (x <= 2 * s, 
          dnorm(x, mean = 0, sd = s), 
          (1 - pnorm(2)) * dunif(x, 2 * s, 4))
}

.temp <- function(x) exp(-1/x) / (x * x)

.dtest <- function(x, p = 0.95, sd = 1, s = 1, log = FALSE) {
  # dens <- p * dnorm(x, mean = 0, sd = s) + (1 - p) * dbeta(x / 4, 1.1, 1.1)
  # dens <- p * dnorm(x, mean = 0, sd = s) + (1 - p) * dnorm(x, 2, 2)
  # dens <- p * dnorm(x, mean = 0, sd = s) + (1 - p) * dnorm(x, 0, 20)
  # dens <- .temp(x, s)
  # dens <- p * dnorm(x, mean = 0, sd = .1) + (1 - p) * dnorm(x, mean = 0, sd = .15)
  dens <- dcauchy(x, location = 0, scale = .1)
  # dens <- dnorm(x, mean = 0, sd = .1)
  if (log) dens <- log(dens)
  return(dens)
}

.d2unif <- function(x, p = 0.95, log = FALSE) {
  dens <- p * dunif(x, -0.1, 0.1) + (1 - p) * dexp(x)
  if (log) dens <- log(dens)
  return(dens)
}

dstudent <- function(x, df, mean = 0, scale = 1, log = FALSE) {
  z <- (x - mean) / scale
  if (log)
    dt(z, df, log = TRUE) - log(scale)
  else
    dt(z, df) / scale
}

dsech <- function(x, location = 0, scale = 1, log = FALSE) {
  if (scale <= 0) stop("scale must be positive")
  z <- (x - location) / scale
  logdens <- -log(2 * scale) - log(cosh(pi * z / 2))
  if (log) logdens else exp(logdens)
}

# Exponentially Modified Gaussian (EMG) PDF
dEMG <- function(y, mu = 0, sigma = 1, lambda = 1, log = FALSE) {
  u <- (mu + lambda * sigma^2 - y) / (sqrt(2) * sigma)
  
  # Compute log density for numerical stability
  log_pdf <- log(lambda) - log(2) +
    .5 * lambda * (2 * mu + lambda * sigma^2 - 2 * y) + pnorm(-u * sqrt(2), log = TRUE)
    # log(pnorm(-u * sqrt(2)))  # since erfc(x) = 2 * pnorm(-x * sqrt(2))

  if (log) return(log_pdf)
  else return(exp(log_pdf))
}

# Function to calculate the PDF of the Normal Inverse Gaussian distribution
# Parameters:
# x: Vector of quantiles.
# alpha: Tail heaviness parameter.
# beta: Skewness parameter.
# delta: Scale parameter.
# mu: Location parameter.
# log: Logical; if TRUE, log-density is returned.

dnig <- function(x, alpha, beta, delta, mu, log = FALSE) {
  
  # Ensure constraints are met (alpha > 0, delta > 0, |beta| < alpha)
  if (alpha <= 0 || delta <= 0 || any(abs(beta) >= alpha)) {
    stop("Parameters must satisfy: alpha > 0, delta > 0, and |beta| < alpha.")
  }
  
  # Pre-calculate common terms
  gamma <- sqrt(alpha^2 - beta^2)
  z <- x - mu
  
  # Calculate K_1 (Modified Bessel function of the third kind, order 1)
  # R uses 'besselK(x, nu)' where nu is the order.
  bessel_term <- besselK(alpha * sqrt(delta^2 + z^2), 1)
  
  # Calculate the NIG density f(x)
  # The formula is: f(x) = (alpha * delta / pi) * K_1(alpha * sqrt(delta^2 + z^2)) * #                     * exp(delta * gamma + beta * z) / sqrt(delta^2 + z^2)
  
  # Term 1: The constant part
  constant_term <- (alpha * delta / pi) * exp(delta * gamma)
  
  # Term 2: The variable part
  variable_term <- bessel_term * exp(beta * z) / sqrt(delta^2 + z^2)
  
  # Full density
  density <- constant_term * variable_term
  
  # Handle log argument
  if (log) {
    return(log(density))
  } else {
    return(density)
  }
}

.nldfrech <- function(x, s = 1, lambda = 1) {
  # x <- abs(x)
  # return(-dexp(x, rate = 1 / s, log = TRUE))
  # return(-dnorm(x, 0, s, log = TRUE))
  # return(-dcauchy(x, 0, s, log = TRUE))
  # return(-sn::dst(x,  xi = 0, omega = 1, alpha = 5, nu = 5, log = TRUE))
  # return(-.dst0(x, omega = s, alpha = 5, nu = 5, log = TRUE))
  # return(-.dmixexpnorm(x, p = 0.5, rate = s, sd = .01, log = TRUE))
  alpha <- 100
  beta <- sqrt(alpha^2 - 1)
  mu <- -s#- s * beta / sqrt(alpha * alpha - beta * beta)
  return(-log(dnig(x, alpha, beta, s, mu)))
  return(-log(.9 * dnorm(x, 0, s) * .1 * dgamma(x, shape = 2)))
  return(-dEMG(x, 0, s, lambda, log = TRUE))
  return(-dsech(x, 0, .5, log = TRUE))
  return(-dstudent(x, 1, 0, .5, log = TRUE))
  return(-dlogis(x, 0, .5, log = TRUE))
  return(-dnorm(x, 0, .5, log = TRUE))
  return(-dcauchy(x, scale = .1, log = TRUE))
  return(-.dtest(x, s = s, log = TRUE))
  # return(-.d2unif(x, log = TRUE))
  out <- log(s) - log(alpha)
  x <- (x - m) / s
  out + (1 + alpha) * log(x) + x^(-alpha)
}

.pend012 <- function(pars, s = .1, lambda = 1, deriv = 0, eps = 1e-4) {
  out <- list()
  f0 <- .nldfrech(pars, s, lambda)
  out[[1]] <- sum(f0)
  if (deriv == 0)
    return(out[[1]])
  ph <- pars + eps
  pl <- pars - eps
  fh <- .nldfrech(ph, s, lambda)
  fl <- .nldfrech(pl, s, lambda)
  out[[2]] <- .5 * (fh - fl) / eps
  out[[3]] <- (fh + fl - 2 * f0) / (eps^2)
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

.d0_Q <- function(pars, likdata, likfns, Q) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  # browser()
  # n_test <- ncol(pm)
  # pm2 <- pm[, 1:n_test]
  # likdata2 <- likdata
  # likdata2$z <- likdata2$z[1:n_test]
  # g1 <- numDeriv::grad(function(x) likfns$d0(matrix(x, 3), likdata2), pm2)
  # g2 <- likfns$d12(pm2, likdata2)[[1]]
  # H1 <- numDeriv::hessian(function(x) likfns$d0(matrix(x, 3), likdata2), pm2)
  # H2 <- likfns$d12(pm2, likdata2)[[2]]
  # g1 <- likfns$d1(as.matrix(pm2), likdata2)
  # numDeriv::hessian(function(x) likfns$d0(matrix(x, 3), likdata2), pm2)[1:3, 1:3]
  # matrix(likfns$d12(as.matrix(pm2), likdata2)[[2]][1, c(1, 2, 3, 2, 4, 5, 3, 5, 6)], 3, 3)
  out <- likdata$mult * likfns$d0(as.matrix(pm), likdata)
  # maybe reinstate this with model = bym4 identifier
  if (!is.null(likdata$bymfns)) {
  # if (any(unlist(likdata$id_bym2))) {
    # temp <- exp(unlist(attr(Q, 'pars')))
    for (i in seq_along(pl)) {
      if (!is.null(likdata$bymfns[[i]])) {
        xi <- pl[[i]][likdata$id_bym2[[i]]]
        parsi <- attr(Q, 'splpars')[[i]][-1]
        if (length(parsi) > 0) {
          parsi <- as.list(parsi)
        }
        out <- out + .pend012(xi, likdata$bymfns[[i]], parsi)
      }
    }
    # out <- out + .pend012(pars[unlist(likdata$id_bym2)], s = temp[2])
  # }
  }
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
  if (any(!is.finite(pars)))
    browser()
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
  # if (likdata$bym4) {
  # # if (any(unlist(likdata$id_bym2))) {
  #   temp <- exp(unlist(attr(Q, 'pars')))
  #   browser()
  #   temp <- .pend012(pars[unlist(likdata$id_bym2)], s = temp[2], deriv = 2)
  #   # temp <- .pend012(pars[unlist(likdata$id_bym2)], s = attr(Q, 'test'), deriv = 2)
  #   id <- unlist(likdata$id_bym2)
  #   g0 <- H0 <- numeric(length(out$g))
  #   g0[id] <- temp[[2]]
  #   out$g <- out$g + g0
  #   H0[id] <- temp[[3]]
  #   H0 <- Matrix::Diagonal(n = length(H0), x = H0)
  #   out$H <- out$H + H0
  # }
  if (!is.null(likdata$bymfns)) {
    # if (any(unlist(likdata$id_bym2))) {
    # temp <- exp(unlist(attr(Q, 'pars')))
    gl <- Hl <- list()
    for (i in seq_along(pl)) {
      if (!is.null(likdata$bymfns[[i]])) {
        parsi <- attr(Q, 'splpars')[[i]][-1]
        if (length(parsi) > 0) {
          parsi <- as.list(parsi)
        }
        temp <- .pend012(pl[[i]][likdata$id_bym2[[i]]], fn = likdata$bymfns[[i]], parsi, deriv = 2)
        gl[[i]] <- Hl[[i]] <- numeric(length(pl[[i]]))
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
  if (any(!is.finite(gH$g)))
    stop('Non-finite gradient')
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

