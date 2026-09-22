.reml0 <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps, direction, kept = NULL) {
  out <- try(.reml0_try(pars, likdata, likfns, Qd, alpha, makeQ, eps, direction, kept = NULL), silent = TRUE)
  if (inherits(out, 'try-error')) {
    out <- 1e20
  }
  out
}

.reml0_try <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps, direction, kept = NULL) {
  beta <- attr(pars, 'beta')
  hyper <- .pars2hyper(pars, Qd$hyper, Qd$hyper_swap)
  Q <- makeQ(hyper, Qd, likdata$control$alpha.tol)
  if (attr(pars, 'first')) 
    likdata$control$inner$itlim <- 1e3
  fit <- .newton_step_inner(beta, .d0_Q, .search_Q, likdata = likdata, likfns = likfns, 
                                    Q = Q, hyper = hyper, control = likdata$control$inner)
  beta1 <- fit$par
  out <- fit$objective
  attr(out, 'unpenalized') <- attr(fit$objective, 'unpenalized')
  attr(out, 'penalized') <- as.numeric(fit$objective)
  out <- out + (1 / likdata$gamma) * .5 * attr(fit$gradient, 'ldet')
  out <- out - (1 / likdata$gamma) * .5 * attr(Q, 'logdet')
  # out <- out + likdata$control$par_mult * sum(unlist(mapply('-', split(pars, Qd$spl), Qd$target))^2)
  attr(out, 'beta') <- fit$par
  attr(out, 'gradient') <- fit$gradient
  attr(out, 'iterations') <- fit$iterations
  attr(out, 'precondHessian') <- fit$precondHessian
  attr(out, 'Hessian') <- fit$Hessian
  attr(out, 'H0') <- fit$H0
  attr(out, 'cholprecondHessian') <- fit$cholprecondHessian
  attr(out, 'diagHessian') <- fit$diagHessian
  attr(out, 'idiagHessian') <- fit$idiagHessian
  attr(out, 'betal') <- NULL
  out
}

.reml1 <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps = 5e-3, direction = 'backward') {
  f0 <- .reml0(pars, likdata, likfns, Qd, alpha, makeQ, eps)
  b0 <- attr(f0, 'beta')
  np <- length(pars)
  if (length(eps) == 1)
    eps <- rep(eps, np)
  if (length(direction) == 1) {
    if (direction == 'ad-hoc') {
      direction <- ifelse(pars > 1, 'backward', 'forward')
    } else {
      direction <- rep(direction, np)
    }
  }
  g <- numeric(np)
  for (i in 1:np) {
    if (direction[i] %in% c('forward', 'central')) {
      pe <- replace(pars, i, pars[i] + eps[i])
      attr(pe, "beta") <- b0
      fu <- .reml0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'forward') {
        g[i] <- as.vector(fu - f0) / eps[i]
      }
    }
    if (direction[i] %in% c('backward', 'central')) {
      pe <- replace(pars, i, pars[i] - eps[i])
      attr(pe, "beta") <- b0
      fl <- .reml0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'backward') {
        g[i] <- as.vector(f0 - fl) / eps[i]
      } else {
        g[i] <- .5 * as.vector(fu - fl) / eps[i]
      }
    }
  }
  g
}

.reml_step <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps = 5e-3, direction = 'backward', kept = NULL) {
  f0 <- .reml0(pars, likdata, likfns, Qd, alpha, makeQ, eps)
  b0 <- attr(f0, 'beta')
  np <- length(pars)
  if (length(eps) == 1)
    eps <- rep(eps, np)
  if (length(direction) == 1) {
    if (direction == 'ad-hoc') {
      direction <- ifelse(pars > 1, 'backward', 'forward')
    } else {
      direction <- rep(direction, np)
    }
  }
  g <- h <- numeric(np)
  for (i in 1:np) {
    if (direction[i] %in% c('forward', 'central')) {
      pe <- replace(pars, i, pars[i] + eps[i])
      attr(pe, "beta") <- b0
      fu <- .reml0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'forward') {
        pe <- replace(pars, i, pars[i] + 2 * eps[i])
        attr(pe, "beta") <- attr(fu, 'beta')
        fuu <- .reml0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
        g[i] <- as.vector(fu - f0) / eps[i]
        h[i] <- as.vector(fuu - 2 * fu + f0) / eps[i] / eps[i]
      }
    }
    if (direction[i] %in% c('backward', 'central')) {
      pe <- replace(pars, i, pars[i] - eps[i])
      attr(pe, "beta") <- b0
      fl <- .reml0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'backward') {
        pe <- replace(pars, i, pars[i] - 2 * eps[i])
        attr(pe, "beta") <- attr(fl, 'beta')
        fll <- .reml0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
        g[i] <- as.vector(f0 - fl) / eps[i]
        h[i] <- as.vector(f0 - 2 * fl + fll) / eps[i] / eps[i]
      } else {
        g[i] <- .5 * as.vector(fu - fl) / eps[i]
        h[i] <- .5 * as.vector(fu - 2 * f0 + fl) / eps[i] / eps[i]
      }
    }
  }
  out <- as.matrix(g / pmax(abs(h), 1))
  attr(out, 'gradient') <- g
  out
}
