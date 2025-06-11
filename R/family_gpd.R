## GEV functions

.quick_tgpd <- function(y) {
  y <- y[is.finite(y)]
  inits <- mean(y, na.rm = TRUE)
  inits <- c(log(inits), -log(.36))
  nlminb(inits, .tgpd0, .tgpd1, .tgpd2, yv = y)$par
}

.d0_gpd <- function(pars_mat, likdata) {
  .tgpdgmrfld0(pars_mat, likdata$z)
}  

# .d1_gpd <- function(pars_mat, likdata) {
#   tgpdgmrfld1(pars_mat, likdata$z)[, 1]
# }
# 
# .d2_gpd <- function(pars_mat, likdata) {
#   tgpdgmrfld2(pars_mat, likdata$z)
# }

.d12_gpd <- function(pars_mat, likdata) {
  gH <- .tgpdgmrfld12(pars_mat, likdata$z)
  list(g = as.vector(gH[, 1:2]), H = gH[, -c(1:2)]) 
}

# .d2_gpdmat <- function(pars_mat, likdata) {
#   tgpdgmrfld2mat(pars_mat, likdata$z)
# }

# .gpd_fns <- list(d0 = .d0_gpd, d1 = .d1_gpd, d2 = .d2_gpd, d12 = .d12_gpd)
.gpd_fns <- list(d0 = .d0_gpd, d12 = .d12_gpd)
.gpd_fns$trans <- list(function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
.gpd_fns$names <- list(link = c('logscale', 'transshape'),
                    response = c('scale', 'shape'))
