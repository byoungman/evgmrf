.reml0 <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps, direction) {
  beta <- attr(pars, 'beta')
  Q <- makeQ(pars, Qd, likdata$control$alpha.tol)
  if (is.null(attr(pars, 'first'))) {
    fit <- .newton(beta, .d0_Q, .search_Q, likdata = likdata, likfns = likfns, Q = Q, stepmax = 3)
  } else {
    fit <- .newton(beta, .d0_Q, .search_Q, likdata = likdata, likfns = likfns, Q = Q, stepmax = 3, itlim = 1e3)
  }
  beta1 <- fit$par
  out <- fit$objective
  out <- out + .5 * attr(fit$gradient, 'ldet')
  out <- out - .5 * attr(Q, 'logdet')
  out <- out + likdata$control$par_mult * sum((pars - Qd$target)^2)
  out <- out + likdata$control$grad_mult * as.numeric(any(abs(fit$gradient) > 1))
  attr(out, 'beta') <- fit$par
  attr(out, 'gradient') <- fit$gradient
  attr(out, 'iterations') <- fit$iterations
  attr(out, 'precondHessian') <- fit$precondHessian
  attr(out, 'Hessian') <- fit$Hessian
  attr(out, 'H0') <- fit$H0
  attr(out, 'cholprecondHessian') <- fit$cholprecondHessian
  attr(out, 'diagHessian') <- fit$diagHessian
  attr(out, 'idiagHessian') <- fit$idiagHessian
  out
}

# .reml0_bfgs <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps, direction) {
#   beta <- attr(pars, 'beta')
#   Q <- makeQ(pars, Qd)
#   if (is.null(attr(pars, 'first'))) {
#     fit <- newton(beta, .d0_Q, .search_Q, likdata = likdata, likfns = likfns, Q = Q, stepmax = 3)
#   } else {
#     fit <- newton(beta, .d0_Q, .search_Q, likdata = likdata, likfns = likfns, Q = Q, stepmax = 3, itlim = 1e3)
#   }
#   beta1 <- fit$par
#   out <- fit$objective
#   out <- out + .5 * attr(fit$gradient, 'ldet')
#   out <- out - .5 * attr(Q, 'logdet')
#   attr(out, 'beta') <- fit$par
#   attr(out, 'gradient') <- fit$gradient
#   attr(out, 'iterations') <- fit$iterations
#   attr(out, 'precondHessian') <- fit$precondHessian
#   attr(out, 'Hessian') <- fit$Hessian
#   attr(out, 'cholprecondHessian') <- fit$cholprecondHessian
#   attr(out, 'diagHessian') <- fit$diagHessian
#   attr(out, 'idiagHessian') <- fit$idiagHessian
#   out
# }

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

.reml_step <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps = 5e-3, direction = 'backward') {
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
