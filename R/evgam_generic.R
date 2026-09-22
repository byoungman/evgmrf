.evgam.control <- function(inner = NULL, outer = NULL) {
  out <- list(outer = list(steptol = 1e-4, 
                           itlim = 1e2, 
                           fntol = 1e-6, 
                           gradtol = 1e-2, 
                           stepmax = 3, 
                           alpha0 = 10, 
                           dgradtol = 1e-4,
                           rho0 = .5),
              inner = list(steptol = 1e-12, 
                           itlim = 1e2, 
                           fntol = 1e-8, 
                           gradtol = 1e-5, 
                           stepmax = 1e2, 
                           alpha0 = 1, 
                           dgradtol = 1e-6,
                           rho0 = .5))
  
  if (!is.null(inner)) {
    if (!is.list(inner)) {
      stop("control argument 'inner' must be a list")
    } else {
      for (i in names(inner)) {
        out$inner[i] <- inner[i]
      }
    }
  }
  
  if (!is.null(outer)) {
    if (!is.list(outer)) {
      stop("control argument 'outer' must be a list")
    } else {
      for (i in names(outer)) {
        out$outer[i] <- outer[i]
      }
    }
  }
  
  out
  
}
