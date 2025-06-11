## r-largest functions

.d0_rlarge <- function(pars_mat, likdata) {
  # rlargegmrfld0(pars_mat, likdata$z_cube)
  .rlargegmrfld0(pars_mat, likdata$z)
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
  gH <- .rlargegmrfld12(pars_mat, likdata$z)
  list(g = as.vector(gH[, 1:3]), H = gH[, -c(1:3)]) 
}

# .rlarge_fns <- list(d0 = .d0_rlarge, d1 = .d1_rlarge, d2 = .d2_rlarge, d12 = .d12_rlarge)
.rlarge_fns <- list(d0 = .d0_rlarge, d12 = .d12_rlarge)
.rlarge_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
.rlarge_fns$names <- list(link = c('location', 'logscale', 'transshape'),
                    response = c('location', 'scale', 'shape'))

.d0_rlargec <- function(pars_mat, likdata) {
  .rlargecgmrfld0(pars_mat, likdata$z, likdata$drop)
}  

# .d1_rlarge <- function(pars_mat, likdata) {
#   rlargegmrfld1(pars_mat, likdata$z)[, 1]
# }
# 
# .d2_rlarge <- function(pars_mat, likdata) {
#   rlargegmrfld2(pars_mat, likdata$z)
# }

.d12_rlargec <- function(pars_mat, likdata) {
  gH <- .rlargecgmrfld12(pars_mat, likdata$z, likdata$drop)
  list(g = as.vector(t(gH[, 1:3])), H = gH[, -c(1:3)]) 
}

.rlargec_fns <- list(d0 = .d0_rlargec, d12 = .d12_rlargec)
.rlargec_fns$trans <- list(function(x) x, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1)
