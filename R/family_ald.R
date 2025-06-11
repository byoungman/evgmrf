## ALD functions

.quick_ald <- function(y, C, tau) {
  psi0 <- sd(y, na.rm = TRUE)
  mu0 <- mean(unlist(y), na.rm = TRUE)
  inits <- c(mu0, log(psi0))
  nlminb(inits, .ald0, .ald1, .ald2, yv = y, C = C, tau = tau)$par
}

.d0_ald <- function(pars_mat, likdata) {
  .aldgmrfld0(pars_mat, likdata$z, likdata$C, likdata$tau)
}  

.d1_ald <- function(pars_mat, likdata) {
  .aldgmrfld1(pars_mat, likdata$z, likdata$C, likdata$tau)[, 1]
}

.d2_ald <- function(pars_mat, likdata) {
  .aldgmrfld2(pars_mat, likdata$z, likdata$C, likdata$tau)
}

.ald_fns <- list(d0 = .d0_ald, d1 = .d1_ald, d2 = .d2_ald)
.ald_fns$trans <- list(function(x) x, function(x) exp(x))
.ald_fns$names <- list(link = c('location', 'logscale'),
                    response = c('location', 'scale'))

