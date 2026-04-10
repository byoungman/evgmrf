.d0_cv <- function (pars, likdata, likfns) {
  here <- sapply(likdata$z, function(x) any(is.finite(x)))
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% 
                                                        pl[[i]])))
  likdata$z <- likdata$z[here]
  likdata$w <- likdata$w[here]
  out <- likdata$mult * likfns$d0(as.matrix(pm[, here, drop = FALSE]), likdata)
  if (!is.finite(out)) 
    out <- 1e+20
  out
}

.cv0 <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps, direction) {
  beta <- attr(pars, 'beta')
  Q <- makeQ(pars, Qd, likdata[[1]]$control$alpha.tol)
  out <- c()
  betas <- list()
  for (i in seq_along(likdata[-1])) {
    if (is.null(attr(pars, 'first'))) {
      fit <- try(.newton(attr(pars, 'betal')[[i]], .d0_Q, .search_Q, likdata = likdata[[i]], likfns = likfns, Q = Q, stepmax = 3), silent = FALSE)
    } else {
      if (is.null(attr(pars, 'betal'))) {
        fit <- try(.newton(beta, .d0_Q, .search_Q, likdata = likdata[[i]], likfns = likfns, Q = Q, stepmax = 3, itlim = 1e3), silent = FALSE)
      } else {
        fit <- try(.newton(attr(pars, 'betal')[[i]], .d0_Q, .search_Q, likdata = likdata[[i]], likfns = likfns, Q = Q, stepmax = 3), silent = FALSE)
      }
    }
    if (!inherits(fit, 'try-error')) {
      if (any(!is.finite(fit$gradient))) {
        out <- c(out, 1e20)
      } else {
        # out <- c(out, .d0_Q(fit$par, likdata$valid[[i]], likfns, 0 * Q))
        out <- c(out, .d0_cv(fit$par, likdata$valid[[i]], likfns))      }
      betas[[i]] <- fit$par
    } else {
      betas[[i]] <- beta
    }
  }
  if (length(out) == 0)
    out <- 1e20
  out <- mean(out)
  # out <- out + .5 * attr(fit$gradient, 'ldet')
  # out <- out - .5 * attr(Q, 'logdet')
  out <- out + .1 * sum(pars^2)
  # out <- out + likdata$control$par_mult * sum((pars - unlist(Qd$target))^2)
  # out <- out + likdata$control$grad_mult * as.numeric(any(abs(fit$gradient) > 1))
  attr(out, 'beta') <- fit$par
  attr(out, 'betal') <- betas
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

.cv1 <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps = 5e-2, direction = 'backward') {
  f0 <- .cv0(pars, likdata, likfns, Qd, alpha, makeQ, eps)
  b0 <- attr(f0, 'beta')
  bl <- attr(f0, 'betal')
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
      attr(pe, "betal") <- bl
      attr(pe, "names") <- attr(pars, "names")
      fu <- .cv0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'forward') {
        g[i] <- as.vector(fu - f0) / eps[i]
      }
    }
    if (direction[i] %in% c('backward', 'central')) {
      pe <- replace(pars, i, pars[i] - eps[i])
      attr(pe, "beta") <- b0
      attr(pe, "betal") <- bl
      attr(pe, "names") <- attr(pars, "names")
      fl <- .cv0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'backward') {
        g[i] <- as.vector(f0 - fl) / eps[i]
      } else {
        g[i] <- .5 * as.vector(fu - fl) / eps[i]
      }
    }
  }
  g
}

.cv_step <- function(pars, likdata, likfns, Qd, alpha, makeQ, eps = 5e-2, direction = 'backward') {
  f0 <- .cv0(pars, likdata, likfns, Qd, alpha, makeQ, eps)
  b0 <- attr(f0, 'beta')
  bl <- attr(f0, 'betal')
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
      attr(pe, "betal") <- bl
      fu <- .cv0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'forward') {
        pe <- replace(pars, i, pars[i] + 2 * eps[i])
        attr(pe, "beta") <- attr(fu, 'beta')
        fuu <- .cv0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
        g[i] <- as.vector(fu - f0) / eps[i]
        h[i] <- as.vector(fuu - 2 * fu + f0) / eps[i] / eps[i]
      }
    }
    if (direction[i] %in% c('backward', 'central')) {
      pe <- replace(pars, i, pars[i] - eps[i])
      attr(pe, "beta") <- b0
      attr(pe, "betal") <- bl
      fl <- .cv0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
      if (direction[i] == 'backward') {
        pe <- replace(pars, i, pars[i] - 2 * eps[i])
        attr(pe, "beta") <- attr(fl, 'beta')
        attr(pe, "betal") <- attr(fl, 'betal')
        fll <- .cv0(pe, likdata, likfns, Qd, alpha, makeQ, eps)
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
