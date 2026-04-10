## ALD functions

.quick_ald <- function(y, args) {
  psi0 <- sd(y, na.rm = TRUE)
  mu0 <- mean(unlist(y), na.rm = TRUE)
  inits <- c(mu0, log(psi0))
  nlminb(inits, .ald0, .ald1, .ald2, yv = y, C = args$C, tau = args$tau)$par
}

.d0_ald <- function(pars_mat, likdata) {
  .aldgmrfld0(pars_mat, likdata$z, likdata$w, likdata$args$C, likdata$args$tau)
}  

# .d1_ald <- function(pars_mat, likdata) {
#   .aldgmrfld1(pars_mat, likdata$z, likdata$w, likdata$C, likdata$tau)[, 1]
# }
# 
# .d2_ald <- function(pars_mat, likdata) {
#   .aldgmrfld2(pars_mat, likdata$z, likdata$w, likdata$C, likdata$tau)
# }

.d12_ald <- function(pars_mat, likdata) {
  if (likdata$openmp) {
    gH <- .aldgmrfld12(pars_mat, likdata$z, likdata$w, likdata$args$C, likdata$args$tau)
  } else {
    gH <- .aldgmrfld12(pars_mat, likdata$z, likdata$w, likdata$args$C, likdata$args$tau)
  }
  out <- list(g = as.vector(gH[, 1:2]), H = gH[, -c(1:2)]) 
  out
}

.J_ald <- function(pars_mat, likdata) {
  .aldgmrfldJ(pars_mat, likdata$z, likdata$w, likdata$args$C, likdata$args$tau)
}

# .ald_fns <- list(d0 = .d0_ald, d1 = .d1_ald, d2 = .d2_ald)
.ald_fns <- list(d0 = .d0_ald, d12 = .d12_ald, J = .J_ald)
.ald_fns$trans <- list(function(x) x, function(x) exp(x))
.ald_fns$names <- list(link = c('location', 'logscale'),
                    response = c('location', 'scale'))

