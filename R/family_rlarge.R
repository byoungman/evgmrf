## r-largest functions

.quick_rlarge <- function(y, delta = 0, hessian = FALSE, derivs = 2) {
  # y <- unlist(y)
  # y_NA <- is.na(y)
  # if (all(y_NA))
  #   return(rep(NA, 3))
  # y <- y[!y_NA]
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
  # rlargegmrfld0(pars_mat, likdata$z_cube)
  if (likdata$openmp) {
    out <- .rlargegmrfld0(pars_mat, likdata$z, likdata$w)
  } else {
    out <- .rlargegmrfld0(pars_mat, likdata$z, likdata$w)
  }
  # browser()
  # n_test <- 2
  # id_test <- as.integer(matrix(1:length(pars_mat), length(likdata$z))[1:n_test, , drop = FALSE])
  # t1 <- likdata$z
  # t1 <- lapply(t1, function(x) x[, 1])
  # t2 <- likdata$w
  # t2 <- lapply(t2, function(x) x[, 1])
  # t3 <- lapply(t1, as.matrix)
  # t4 <- lapply(t2, as.matrix)
  # t5 <- .tgevgmrfld0(pars_mat[, 1:n_test], t1[1:n_test], t2[1:n_test])
  # t6 <- .tgevgmrfld12(pars_mat[, 1:n_test], t1[1:n_test], t2[1:n_test])
  # t7 <- .rlargegmrfld0(pars_mat[, 1:n_test], t3[1:n_test], t4[1:n_test])
  # t7 - t5
  # t8 <- .rlargegmrfld12(pars_mat[, 1:n_test], t3[1:n_test], t4[1:n_test])
  # t8 - t6
  # g1 <- numDeriv::grad(function(x) .rlargegmrfld0(matrix(x, 3), likdata$z[1:n_test], likdata$w[1:n_test]), as.vector(pars_mat[, 1:n_test]))
  # g2 <- t(.rlargegmrfld12(pars_mat[, 1:n_test], likdata$z[1:n_test], likdata$w[1:n_test])[, 1:3])
  # g3 <- numDeriv::grad(function(x) .rlargegmrfld0(matrix(x, 3), t3[1:n_test], t4[1:n_test]), as.vector(pars_mat[, 1:n_test]))
  # g4 <- t(.rlargegmrfld12(pars_mat[, 1:n_test], t3[1:n_test], t4[1:n_test])[, 1:3])
  # h1 <- numDeriv::hessian(function(x) .rlargegmrfld0(matrix(x, 3), likdata$z[1:n_test], likdata$w[1:n_test]), as.vector(pars_mat[, 1:n_test]))
  # h2 <- t(.rlargegmrfld12(pars_mat[, 1:n_test], likdata$z[1:n_test], likdata$w[1:n_test])[, -c(1:3)])
  # 
  # plot(g3, g4)
  out
}

# .d1_rlarge <- function(pars_mat, likdata) {
#   rlargegmrfld1(pars_mat, likdata$z)[, 1]
# }
# 
# .d2_rlarge <- function(pars_mat, likdata) {
#   rlargegmrfld2(pars_mat, likdata$z)
# }
# 
# .d2_rlargemat <- function(pars_mat, likdata) {
#   rlargegmrfld2mat(pars_mat, likdata$z)
# }

.d12_rlarge <- function(pars_mat, likdata) {
  # gH <- rlargegmrfld12(pars_mat, likdata$z_cube)
  if (likdata$openmp) {
    gH <- .rlargegmrfld12(pars_mat, likdata$z, likdata$w)
  } else {
    gH <- .rlargegmrfld12(pars_mat, likdata$z, likdata$w)
  }
  list(g = as.vector(gH[, 1:3]), H = gH[, -c(1:3)]) 
}

# .rlarge_fns <- list(d0 = .d0_rlarge, d1 = .d1_rlarge, d2 = .d2_rlarge, d12 = .d12_rlarge)
.rlarge_fns <- list(d0 = .d0_rlarge, d12 = .d12_rlarge)
.rlarge_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
.rlarge_fns$names <- list(link = c('location', 'logscale', 'transshape'),
                    response = c('location', 'scale', 'shape'))

# ## censored
# 
# .d0_rlargec <- function(pars_mat, likdata) {
#   .rlargecgmrfld0(pars_mat, likdata$z, likdata$drop)
# }  
# 
# # .d1_rlarge <- function(pars_mat, likdata) {
# #   rlargegmrfld1(pars_mat, likdata$z)[, 1]
# # }
# # 
# # .d2_rlarge <- function(pars_mat, likdata) {
# #   rlargegmrfld2(pars_mat, likdata$z)
# # }
# 
# .d12_rlargec <- function(pars_mat, likdata) {
#   gH <- .rlargecgmrfld12(pars_mat, likdata$z, likdata$drop)
#   list(g = as.vector(t(gH[, 1:3])), H = gH[, -c(1:3)]) 
# }
# 
# .rlargec_fns <- list(d0 = .d0_rlargec, d12 = .d12_rlargec)
# .rlargec_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
