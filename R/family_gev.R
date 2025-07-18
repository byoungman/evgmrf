## GEV functions

.quick_tgev <- function(y, delta) {
  y <- unlist(y)
  y_NA <- is.na(y)
  if (all(y_NA))
    return(rep(NA, 3))
  y <- y[!y_NA]
  psi0 <- sqrt(6 * var(y, na.rm = TRUE)) / pi
  mu0 <- mean(y, na.rm = TRUE) - 0.57722 * psi0
  inits <- c(mu0, log(psi0), -log(.36))
  nlminb(inits, .tgev0, .tgev1, .tgev2, yv = y, delta = delta)$par
}

.d0_gev <- function(pars_mat, likdata) {
  if (likdata$openmp) {
    out <- .tgevgmrfld0_omp(pars_mat, likdata$z, likdata$w, likdata$threads)
  } else {
    out <- .tgevgmrfld0(pars_mat, likdata$z, likdata$w)
  } 
  out
}  
# 
# .d1_gev <- function(pars_mat, likdata) {
#   .tgevgmrfld1(pars_mat, likdata$z, likdata$w)[, 1]
# }
# 
# .d2_gev <- function(pars_mat, likdata) {
#   .tgevgmrfld2(pars_mat, likdata$z, likdata$w)
# }

.d12_gev <- function(pars_mat, likdata) {
  if (likdata$openmp) {
    gH <- .tgevgmrfld12_omp(pars_mat, likdata$z, likdata$w, likdata$threads)
  } else {
    gH <- .tgevgmrfld12(pars_mat, likdata$z, likdata$w)
  }
  out <- list(g = as.vector(gH[, 1:3]), H = gH[, -c(1:3)]) 
  out
}

.J_gev <- function(pars_mat, likdata) {
.tgevgmrfldJ(pars_mat, likdata$z, likdata$w)
}

# 
# .d2_gevmat <- function(pars_mat, likdata) {
#   .tgevgmrfld2mat(pars_mat, likdata$z, likdata$w)
# }

# .gev_fns <- list(d0 = .d0_gev, d1 = .d1_gev, d2 = .d2_gev, d12 = .d12_gev)
.gev_fns <- list(d0 = .d0_gev, d12 = .d12_gev, J = .J_gev)
.gev_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)

attr(.gev_fns$trans, 'deriv') <- list(function(x) 1 + 0 * x, 
                                      function(x) exp(x), 
                                      function(x) 1.5 * exp(-x)/(1 + exp(-x))^2)

.gev_fns$names <- list(link = c('location', 'logscale', 'transshape'),
                    response = c('location', 'scale', 'shape'))

.qgev <- function(p, loc, scale, shape) {
  shape[shape == 0] <- 1e-06
  shape <- sign(shape) * pmax(abs(shape), 1e-06)
  yp <- -log(p)
  loc - scale * (1 - yp^(-shape))/shape
}

.qgev0 <- function(p, location, logscale, transshape) {
  scale <- exp(logscale)
  shape <- 1.5 / (1 + exp(-transshape)) - 1
#  shape <- sign(shape) * pmax(abs(shape), 1e-06)
  yp <- -log(p)
  location - scale * (1 - yp^(-shape))/shape
}

.qgev0_d0 <- function (p, location, logscale, transshape) {
  .e2 <- exp(-transshape)
  .e3 <- 1 + .e2
  .e5 <- 1.5/.e3 - 1
  .e6 <- -log(p)
  .e7 <- .e6^.e5
  .e8 <- 1 - 1/.e7
  .e9 <- exp(logscale)
  cbind(location = 1 + 0 * location, 
        logscale = -(.e8 * .e9/.e5), 
        transshape = -(1.5 * (.e2 * .e9 * (log(.e6)/.e7 - .e8/.e5)/(.e3^2 * .e5))))
}

attr(.qgev0, 'deriv') <- .qgev0_d0

.gev_fns$quantile <- function(p, location, scale, shape) .qgev(p, location, scale, shape)
.gev_fns$quantile0 <- .qgev0

.tgev0_shrink <- function(pars, yv, delta, pars0, mult) {
  out <- .tgev0(pars, yv, delta)
  out <- out + mult * sum((pars - pars0)^2)
  out
}

.tgev1_shrink <- function(pars, yv, delta, pars0, mult) {
  out <- .tgev1(pars, yv, delta)
  out <- out + mult * 2 * (pars - pars0)
  out
}

.tgev2_shrink <- function(pars, yv, delta, pars0, mult) {
  out <- .tgev2(pars, yv, delta)
  diag(out) <- diag(out) + 2 * mult
  out
}

.quick_tgev_shrink <- function(y, delta, pars0, mult = 10) {
  y <- unlist(y)
  y_NA <- is.na(y)
  if (all(y_NA))
    return(rep(NA, 3))
  y <- y[!y_NA]
  psi0 <- sqrt(6 * var(y, na.rm = TRUE)) / pi
  mu0 <- mean(y, na.rm = TRUE) - 0.57722 * psi0
  inits <- c(mu0, log(psi0), -log(.36))
  nlminb(inits, .tgev0_shrink, .tgev1_shrink, .tgev2_shrink, yv = y, delta = delta, 
         pars0 = pars0, mult = mult)$par
}


