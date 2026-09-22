## r-largest functions

.quick_rlarge <- function(y, delta = 0, hessian = FALSE, derivs = 2) {
  psi0 <- sqrt(6 * var(y[, 1], na.rm = TRUE)) / pi
  mu0 <- mean(y[, 1], na.rm = TRUE) - 0.57722 * psi0
  inits <- c(mu0, log(psi0), -log(.36))
  if (derivs == 2)
    out <- nlminb(inits, .rlarged0, .rlarged1, .rlarged2, yv = y)$par
  if (derivs == 1)
    out <- nlminb(inits, .rlarged0, .rlarged1, yv = y)$par
  if (derivs == 0)
    out <- nlminb(inits, .rlarged0, yv = y)$par
  if (hessian)
    attr(out, 'hessian') <- .rlarged2(out, y)
  out
}

.d0_rlarge <- function(pars_mat, likdata) {
  if (likdata$openmp) {
    out <- .rlargegmrfld0_omp(pars_mat, likdata$z, likdata$w, likdata$threads)
  } else {
    out <- .rlargegmrfld0(pars_mat, likdata$z, likdata$w)
  }
  out
}

.d12_rlarge <- function(pars_mat, likdata) {
  if (likdata$openmp) {
    gH <- .rlargegmrfld12_omp(pars_mat, likdata$z, likdata$w, likdata$threads)
  } else {
    gH <- .rlargegmrfld12(pars_mat, likdata$z, likdata$w)
  }
  list(g = as.vector(gH[, 1:3]), H = gH[, -c(1:3)]) 
}

.rlarge_fns <- list(d0 = .d0_rlarge, d12 = .d12_rlarge)
.rlarge_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
.rlarge_fns$names <- list(link = c('location', 'logscale', 'transshape'),
                    response = c('location', 'scale', 'shape'))