## Point process functions

.quick_tpp <- function(y, m, u, delta) {
  if (is.null(delta))
    delta <- 1
  mu0 <- min(y)
  xi0 <- max(.5 * (1 - (mean(y) - mu0)^2 / var(y)), -.99)
  psi0 <- (mean(y) - mu0) * (1 - xi0)
  txi0 <- -log(1.5 / (1 + xi0) - 1)
  inits <- c(mu0, log(psi0), -log(.36))
  fit <- nlminb(inits, .tpp0, .tpp1, .tpp2, zv = y, w = m, u = u, delta = delta)
  fit$par
}

.d0_pp <- function(pars_mat, likdata) {
  out <- .tppugmrfld0(pars_mat, likdata$u, likdata$uw)
  out <- out + .tppzgmrfld0(pars_mat, likdata$z, likdata$w)
}

.d12_pp <- function(pars_mat, likdata) {
  if (likdata$openmp) {
    gH <- .tppugmrfld12(pars_mat, likdata$u, likdata$uw)
    gH <- gH + .tppzgmrfld12(pars_mat, likdata$z, likdata$w)
  } else {
    gH <- .tppugmrfld12(pars_mat, likdata$u, likdata$uw)
    gH <- gH + .tppzgmrfld12(pars_mat, likdata$z, likdata$w)
  }
  list(g = as.vector(gH[, 1:3]), H = gH[, -c(1:3)]) 
}

.pp_fns <- list(d0 = .d0_pp, d12 = .d12_pp)
.pp_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
.pp_fns$names <- list(link = c('location', 'logscale', 'transshape'),
                    response = c('location', 'scale', 'shape'))

.pp_fns$quantile <- function(p, location, scale, shape) .qgev(p, location, scale, shape)
.pp_fns$quantile0 <- .qgev0
