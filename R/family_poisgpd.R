## Point process functions

.quick_tpp <- function(y, m, u, delta) {
  #  psi0 <- sqrt(6 * var(unlist(y), na.rm = TRUE)) / pi
  #  mu0 <- mean(unlist(y), na.rm = TRUE) - 0.57722 * psi0
  #  inits <- c(mu0, log(psi0), -log(.36))
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
  # browser()
  # check
  # n_test <- 2
  # id_test <- as.integer(matrix(1:length(pars_mat), length(likdata$z))[1:n_test, , drop = FALSE])
  # g1 <- numDeriv::grad(function(x) .tppzgmrfld0(matrix(x, 3), likdata$z[1:n_test], likdata$w[1:n_test]), as.vector(pars_mat[, 1:n_test]))
  # out <- t(.tppzgmrfld12(pars_mat[, 1:n_test], likdata$z[1:n_test], likdata$w[1:n_test])[, 1:3])
  # g1 <- numDeriv::grad(function(x) .tppugmrfld0(matrix(x, 3), likdata$u[1:n_test], likdata$uw[1:n_test]), as.vector(pars_mat[, 1:n_test]))
  # g2 <- t(.tppugmrfld12(pars_mat[, 1:n_test], likdata$u[1:n_test], likdata$uw[1:n_test])[, 1:3])
  # evgmrf:::.tppu1(pars_mat[, 1], likdata$u[1], likdata$uw[1])
  out
}

# .d1_pp <- function(pars_mat, likdata) {
#   g <- tppugmrfld1(pars_mat, likdata$u, likdata$m)[, 1]
#   # g <- g + tppzgmrfld1(pars_mat, likdata$z)[, 1]
#   g
# }
# 
# .d2_pp <- function(pars_mat, likdata) {
#   H <- tppugmrfld2(pars_mat, likdata$u, likdata$m)
#   # H <- H + tppzgmrfld2(pars_mat, likdata$z)
#   H
# }

# .d12_pp <- function(pars_mat, likdata) {
#   list(g = .d1_pp(pars_mat, likdata), H = .d2_pp(pars_mat, likdata))
# }

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

# .pp_fns <- list(d0 = .d0_pp, d1 = .d1_pp, d2 = .d2_pp, d12 = .d12_pp)
.pp_fns <- list(d0 = .d0_pp, d12 = .d12_pp)
.pp_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
.pp_fns$names <- list(link = c('location', 'logscale', 'transshape'),
                    response = c('location', 'scale', 'shape'))

