## Poisson process functions

.quick_pproc <- function(y) {
  log(mean(y))
}

.d0_pproc <- function(pars_mat, likdata) {
  .tpprocgmrfld0(pars_mat, likdata$z, likdata$mult)
}  

.d12_pproc <- function(pars_mat, likdata) {
  gH <- .tpprocgmrfld12(pars_mat, likdata$z, likdata$mult)
  list(g = as.vector(gH[, 1]), H = gH[, -1]) 
}

.pproc_fns <- list(d0 = .d0_pproc, d12 = .d12_pproc)
.pproc_fns$trans <- list(function(x) exp(x))
.pproc_fns$names <- list(link = c('lograte'),
                    response = c('rate'))
