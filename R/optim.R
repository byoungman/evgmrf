.newton <- function(pars, fn, sfn, ..., steptol = 1e-12, itlim = 1e2, fntol = 1e-12, gradtol = 1e-4, stepmax = 1e4, alpha0 = 1, trace = FALSE, reml = FALSE) {
  
  pars0 <- pars
  
  it <- 1
  okay <- TRUE
  f0 <- fn(pars, ...)
  b0 <- attr(f0, 'beta')
  step1 <- NULL
  
  if (itlim == 1) {
    okay <- FALSE
    report <- 'Fixed parameters.'
    attr(pars, 'beta') <- b0
  }
  
  while (okay) {
    if (it > 1) 
      g0 <- g
    if (!is.null(step1)) {
      step0 <- step1
    } else {
      attr(pars, 'beta') <- b0
      step0 <- sfn(pars, ...)
    }
    vstep0 <- as.vector(step0)
    g <- attr(step0, "gradient")
    H <- attr(step0, "Hessian")
    if (trace) {
      if (it > 1) cat("\n")
      cat("Iteration: ", it - 1, "\n")
      cat("Value: ", as.vector(f0), "\n")
      cat(paste('Gradient range: (', paste0(signif(range(g), 4), collapse = ', '), ')', sep = ''))
      if (reml) {
        cat('\n')
        cat(paste('Parameters: (', paste0(signif(as.vector(pars), 4), collapse = ', '), ')', sep = ''))
        cat('\n')
        cat(paste('Inner max |grad|:', signif(max(abs(attr(f0, 'gradient'))), 4)))
        cat('\n')
      }
    }
    if (max(abs(g)) < gradtol) {
      report <- c("gradient tolerance reached")
      break
    }
    # too_big <- abs(step0[, 1]) > stepmax
    # if (any(too_big)) {
    #   step0[which(too_big)] <- sign(step0[which(too_big)]) * stepmax
    # }
    vstep0 <- sign(vstep0) * pmin(abs(vstep0), stepmax)
    alpha <- alpha0
    report <- NULL
    ls <- TRUE
    while(ls & is.null(report)) {
      step <- alpha * vstep0
      stepokay <- mean(abs(step)) > steptol
      if (!stepokay) {
        report <- c("step tolerance reached")
      } else {
        theta1 <- pars - step
        attr(theta1, 'beta') <- b0
        f1 <- fn(theta1, ...)
        d <- f1 - f0
        if (!is.finite(d)) 
          d <- 10
        if (d < 0) {
          attr(theta1, "beta") <- attr(f1, "beta")
          step1 <- try(sfn(theta1, ...), silent = TRUE)
          if (inherits(step1, "try-error")) 
            d <- 1
          if (any(!is.finite(attr(step1, "gradient")))) 
            d <- 1
          if (any(!is.finite(step1))) 
            d <- 1      
        }  
        if (d < 0) {
          f0 <- f1
          b0 <- attr(f0, 'beta')
          pars <- theta1
          ls <- FALSE
        } else {
          if (d < fntol) 
            report <- c("function tolerance reached")
          alpha <- ifelse(alpha == 1, .1 * alpha, .5 * alpha)
        }
      }
    }
    if (!is.null(report)) 
      break
    it <- it + 1
    if (it == itlim) {
      report <- c("iteration limit reached")
      okay <- FALSE
    }
  }
  
  if (trace) 
    cat(paste("\n ", it, "iterations:", report, "\n"))
  out <- list(pars=as.vector(pars), objective=f0)
  if (reml) {
    out$beta <- attributes(f0)$beta
    out$gradient <- attributes(f0)$gradient
    out$precondHessian <- attributes(f0)$precondHessian
    out$Hessian <- attributes(f0)$Hessian
    out$H0 <- attributes(f0)$H0
    out$cholprecondHessian <- attributes(f0)$cholprecondHessian
    out$diagHessian <- attributes(f0)$diagHessian
    out$idiagHessian <- attributes(f0)$idiagHessian
  } else {
    out$gradient <- g
    out$Hessian <- H
    out$H0 <- attr(step0, 'H0')
    out$precondHessian <- attr(step0, "precondHessian")
    out$cholprecondHessian <- attr(step0, "cholprecondHessian")
    out$diagHessian <- Matrix::diag(attr(step0, "diagHessian"))
    out$idiagHessian <- attr(step0, "idiagHessian")
    out$rankHessian <- attr(step0, "rank")
  }
  out$convergence <- 0
  out$report <- report
  out$iterations <- it
  out$gradconv <- substr(report, 1, 4) == "grad"
  if (!is.null(attr(pars, "beta"))) 
    out$beta <- attr(pars, "beta")
  out
}

.nelder_mead_list <- function(init, f, xtol = 1e-3, ftol = 1e-4, max_iter = 500, trace = 1, ...) {
  n <- length(init)  # Dimension of the problem
  alpha <- 1         # Reflection coefficient
  gamma <- 2         # Expansion coefficient
  rho <- 0.5         # Contraction coefficient
  sigma <- 0.5       # Shrink coefficient
  delta <- 0.5         # initial simplex perturbation
  
  # Initialize the simplex
  simplex <- list()
  simplex[[1]] <- init
  f_values_list <- list()
  f_values_list[[1]] <- f(init, ...)
  
  attr(init, 'beta') <- attr(f_values_list[[1]], 'beta')
  
  if (trace > 0) {
    cat(paste('Iteration:', 0))
    cat('\n')
    cat(paste('Value:', signif(f_values_list[[1]], 8)))
    cat('\n')
    cat(paste('Initial values: (', paste0(signif(init, 4), collapse = ', '), ')', sep = ''))
    cat('\n')
    cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
    cat('\n')
    cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
    cat('\n')
    cat('\n')
  }
  
  for (i in 1:n) {
    xlo <- xhi <- init
    xlo[i] <- xlo[i] - delta  # Perturbation by -delta
    flo <- f(xlo, ...)
    xhi[i] <- xhi[i] + delta  # Perturbation by delta
    fhi <- f(xhi, ...)
    if (flo < fhi) {
      simplex[[i + 1]] <- xlo
      f_values_list[[i + 1]] <- flo
    } else {
      simplex[[i + 1]] <- xhi
      f_values_list[[i + 1]] <- fhi
    }
  }
  # browser()
  # f_values_list0 <- lapply(simplex, function(x) lapply(x, f, ...))
  # Evaluate function at simplex points
  # f_values_list <- lapply(simplex, f, ...)
  f_values <- sapply(f_values_list, as.vector)
  
  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1
    
    # Order simplex points by function values
    order_idx <- order(f_values)
    simplex <- simplex[order_idx]
    f_values <- f_values[order_idx]
    f_values_list <- f_values_list[order_idx]
    b0 <- attr(f_values_list[[1]], 'beta')
    for (i in seq_along(simplex)) attr(simplex[[i]], 'beta') <- b0
    
    # Centroid of all points except the worst
    # centroid <- colMeans(simplex[1:n, ])
    centroid <- rowMeans(do.call(cbind, simplex[1:n]))
    
    if (trace > 0) {
      cat(paste('Iteration:', iter))
      cat('\n')
      cat(paste('Value range: (', paste0(signif(range(f_values), 8), collapse = ', '), ')', sep = ''))
      cat('\n')
      cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
      cat('\n')
      cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
      cat('\n')
      cat(paste('Centroid: (', paste0(signif(simplex[[1]], 4), collapse = ', '), ')', sep = ''))
      cat('\n')
      cat('\n')
    }
    
    # Reflection
    x_reflect <- centroid + alpha * (centroid - simplex[[n + 1]])
    f_reflect <- f(x_reflect, ...)
    
    if (f_reflect < f_values[1]) {
      # Expansion
      x_expand <- centroid + gamma * (x_reflect - centroid)
      f_expand <- f(x_expand, ...)
      if (f_expand < f_reflect) {
        simplex[[n + 1]] <- x_expand
        f_values[n + 1] <- f_expand
        f_values_list[[n + 1]] <- f_expand
      } else {
        simplex[[n + 1]] <- x_reflect
        f_values[n + 1] <- f_reflect
        f_values_list[[n + 1]] <- f_reflect
      }
    } else if (f_reflect < f_values[n]) {
      # Accept reflection
      simplex[[n + 1]] <- x_reflect
      f_values[n + 1] <- f_reflect
      f_values_list[[n + 1]] <- f_reflect
    } else {
      # Contraction
      if (f_reflect < f_values[n + 1]) {
        # Outside contraction
        x_contract <- centroid + rho * (x_reflect - centroid)
      } else {
        # Inside contraction
        x_contract <- centroid + rho * (simplex[[n + 1]] - centroid)
      }
      f_contract <- f(x_contract, ...)
      
      if (f_contract < f_values[n + 1]) {
        simplex[[n + 1]] <- x_contract
        f_values[n + 1] <- f_contract
        f_values_list[[n + 1]] <- f_contract
      } else {
        # Shrink the simplex
        for (i in 2:(n + 1)) {
          simplex[[i]] <- simplex[[1]] + sigma * (simplex[[i]] - simplex[[1]])
          f_values_list[[i]] <- f(simplex[[i]], ...)
          f_values[i] <- f_values_list[[i]]
          
        }
      }
      
    }
    
    # Check convergence
    if ((f_values[n] - f_values[1]) / abs(f_values[1]) < ftol) {
      break
    }
    if (mean(sapply(seq_along(init), function(i) diff(range(sapply(simplex, '[', i))))) < xtol) {
      break
    }
  }
  
  f1 <- f_values[1]
  attr(f1, 'beta') <- attr(simplex[[1]], 'beta')
  list(par = simplex[[1]], objective = f1, iterations = iter, beta = attr(simplex[[1]], 'beta'))
}

.itreport <- function(f, g, it) {
  report <- paste("\n Outer iteration ", it, ":", sep="")
  rep1 <- paste("  Outer max(|grad|):", signif(max(abs(g)), 3))
  rep2 <- paste("  Inner max(|grad|): ", signif(max(abs(attr(f, "gradient"))), 3), ".", sep="")
  report <- c(report, paste(rep1, rep2, sep="; "))
  cat(paste(report, collapse="\n"))
}

.BFGS <- function(pars, fn, gfn, ..., steptol = 1e-12, itlim = 1e2, fntol = 1e-8, gradtol = 1e-2, stepmax = 3, alpha = 1, trace=0) {
  
  it <- 1
  okay <- TRUE
  f0 <- fn(pars, ...)
  g1 <- NULL
  I <- iH <- H <- diag(length(pars))
  
  while (okay) {
    if (it > 1) g0 <- g
    if (!is.null(g1)) {
      g <- g1
    } else {
      attr(pars, "beta") <- attr(f0, "beta")
      attr(pars, "betal") <- attr(f0, "betas")
      g <- gfn(pars, ...)
      iH <- diag(length(g))#diag(rep(1 / max(abs(g)), length(g)))
    }
    if (trace) .itreport(f0, g, it - 1)
    if (mean(abs(g)) < gradtol) {
      report <- c("gradient tolerance reached")
      break
    }
    step0 <- crossprod(iH, g)
    step0 <- sign(step0) * pmin(abs(step0), stepmax)
    report <- NULL
    ls <- TRUE
    while(ls & is.null(report)) {
      step <- alpha * step0
      stepokay <- all(abs(step) > steptol)
      if (!stepokay) {
        report <- c("step tolerance reached")
      } else {
        theta1 <- pars - step
        f1 <- fn(theta1, ...)
        d <- f1 - f0
        if (d < 0) {
          attr(theta1, "beta") <- attr(f1, "beta")
          attr(theta1, "betal") <- attr(f1, "betas")
          g1 <- gfn(theta1, ...)
          if (any(!is.finite(g1))) d <- 1
          yk <- g1 - g
          denom <- sum(- yk * step)
          t1 <- I - tcrossprod(- step, yk) / denom
          t2 <- I - tcrossprod(yk, - step) / denom
          t3 <- tcrossprod(- step) / denom
          iH <- t1 %*% iH %*% t2 + t3
          if (any(!is.finite(iH))) d <- 1
        }
        if (d < 0) {
          f0 <- f1
          pars <- theta1
          ls <- FALSE
        } else {
          if (d < fntol) {
            report <- c("function tolerance reached")
          }
          alpha <- .5 * alpha
        }
      }
    }
    if (!is.null(report)) break
    it <- it + 1
    if (it == itlim) {
      report <- c("iteration limit reached")
      okay <- FALSE
    }
  }
  if (trace) cat(paste("\n ", it, "iterations:", report, "\n"))
  out <- list(par = as.vector(pars), objective = as.vector(f0))
  out$gradient <- g
  out$convergence <- 0
  out$report <- report
  out$iterations <- it
  # if (!is.null(attr(pars, "beta"))) out$beta <- attr(pars, "beta")
  # if (!is.null(attr(pars, "betal"))) out$betal <- attr(pars, "betal")
  out$beta <- attributes(f0)$beta
  out$gradient <- attributes(f0)$gradient
  out$precondHessian <- attributes(f0)$precondHessian
  out$Hessian <- attributes(f0)$Hessian
  out$H0 <- attributes(f0)$H0
  out$cholprecondHessian <- attributes(f0)$cholprecondHessian
  out$diagHessian <- attributes(f0)$diagHessian
  out$idiagHessian <- attributes(f0)$idiagHessian
  out
}

.brent <- function(x0, f, step = 2, tol = 1e-4, max_iter = 100, trace, ...) {
  # Phase 1: Interval expansion to bracket the minimum
  a <- x0
  b <- x0 + step
  fa <- f(a, ...)
  fb <- f(b, ...)
  
  if (fa < fb) {
    # Search in the opposite direction
    a <- x0 - step
    b <- x0
    fa <- f(a, ...)
    fb <- f(b, ...)
  }
  
  # Expand until bracketing interval [a, b] with fa > fb
  while (fa <= fb) {
    step <- 2 * step
    a <- b
    b <- b + step
    fa <- fb
    fb <- f(b, ...)
    if (abs(b) > 1e2) stop("Failed to find bounds; possible divergence.")
  }
  
  # Phase 2: Brent's method within bounds [a, b]
  gr <- (3 - sqrt(5)) / 2  # Golden ratio
  v <- a + gr * (b - a)
  w <- v
  x <- v
  fx <- f(x, ...)
  
  if (trace > 0) {
    cat(paste('Iteration:', 0))
    cat('\n')
    cat(paste('Value:', signif(fx, 8)))
    cat('\n')
    cat(paste('Parameter:', signif(x, 4)))
    cat('\n')
    cat('\n')
  }
  
  fv <- fx
  fw <- fx
  
  d <- 0
  e <- 0
  
  for (iter in 1:max_iter) {
    m <- 0.5 * (a + b)
    tol1 <- tol * abs(x) + tol / 10
    tol2 <- 2 * tol1
    
    if (trace > 0) {
      cat(paste('Iteration:', iter))
      cat('\n')
      cat(paste('Value:', signif(fx, 8)))
      cat('\n')
      cat(paste('Parameter:', signif(x, 4)))
      cat('\n')
      cat('\n')
    }
    
    # Convergence check
    if (abs(x - m) <= tol2 - 0.5 * (b - a)) {
      attr(fx, 'beta') <- attr(x, 'beta')
      return(list(par = x, objective = fx, iterations = iter, beta = attr(x, 'beta')))
    }
    
    # Parabolic fit
    p <- q <- 0
    if (abs(e) > tol1) {
      r <- (x - w) * (fx - fv)
      q <- (x - v) * (fx - fw)
      p <- (x - v) * q - (x - w) * r
      q <- 2 * (q - r)
      if (q > 0) p <- -p
      q <- abs(q)
      if (abs(p) < abs(0.5 * q * e) && p > q * (a - x) && p < q * (b - x)) {
        d <- p / q
        u <- x + d
        if (u - a < tol2 || b - u < tol2) {
          d <- ifelse(x < m, tol1, -tol1)
        }
      } else {
        e <- ifelse(x < m, b - x, a - x)
        d <- gr * e
      }
    } else {
      e <- ifelse(x < m, b - x, a - x)
      d <- gr * e
    }
    
    u <- ifelse(abs(d) >= tol1, x + d, x + ifelse(d > 0, tol1, -tol1))
    attr(u, 'beta') <- attr(x, 'beta')
    fu <- f(u, ...)
    
    # Update points
    if (fu <= fx) {
      if (u < x) b <- x else a <- x
      v <- w; fv <- fw
      w <- x; fw <- fx
      x <- u; fx <- fu
    } else {
      if (u < x) a <- u else b <- u
      if (fu <= fw || w == x) {
        v <- w; fv <- fw
        w <- u; fw <- fu
      } else if (fu <= fv || v == x || v == w) {
        v <- u; fv <- fu
      }
    }
  }
  
  warning("Maximum number of iterations reached")
  #  return(list(minimum = x, value = fx, iterations = max_iter))
  attr(fx, 'beta') <- attr(x, 'beta')
  list(par = x, objective = fx, iterations = iter, beta = attr(x, 'beta'))
}

.newton_diag <- function(fn, gr, hess_diag, x0, max_iter = 100, tol = 2e-4,
                         alpha_init = 1.0, beta = 0.5, c1 = 1e-4, ...) {
  x <- x0
  fx <- fn(x, ...)
  g <- gr(x, ...)
  
  for (k in 1:max_iter) {
    # if (sqrt(sum(g^2)) < tol) {
    if (median(abs(g)) < tol) {
      
      # cat("Converged in", k, "iterations\n")
      return(x)
    }
    
    D <- hess_diag(x, ...)
    if (any(D == 0)) {
      warning("Zero in diagonal Hessian - skipping iteration")
      next
    }
    
    p <- -g / D  # approximate Newton step using diag(H)^-1
    
    # Armijo backtracking line search
    step <- alpha_init
    while (TRUE) {
      x_new <- x + step * p
      f_new <- fn(x_new, ...)
      if (!is.finite(f_new)) {
        step <- step * beta
      } else if (f_new <= fx + c1 * step * sum(g * p)) {
        break
      } else {
        step <- step * beta
      }
      
      if (step < 1e-10) {
        warning("Line search failed - step too small")
        break
      }
    }
    
    x <- x + step * p
    fx <- fn(x, ...)
    g <- gr(x, ...)
  }
  
  warning("Max iterations reached without convergence")
  return(x)
}

# .nelder_mead_discrete_list <- function(init, f, xtol = 1e-3, ftol = 1e-4, max_iter = 500, trace = 1, step_size = 0.2, ...) {
#   n <- length(init)  # Dimension of the problem
#   alpha <- 1         # Reflection coefficient
#   gamma <- 2         # Expansion coefficient
#   rho <- 0.5         # Contraction coefficient
#   sigma <- 0.5       # Shrink coefficient
#   delta <- 0.5       # initial simplex perturbation
#   
#   # Function to round to nearest multiple of step_size
#   round_to_step <- function(x) {
#     return(round(x / step_size) * step_size)
#   }
#   
#   # Round initial point to nearest multiple of step_size
#   init <- round_to_step(init)
#   
#   # Initialize the simplex
#   simplex <- list()
#   simplex[[1]] <- init
#   f_values_list <- list()
#   f_values_list[[1]] <- f(init, ...)
#   
#   attr(init, 'beta') <- attr(f_values_list[[1]], 'beta')
#   
#   if (trace > 0) {
#     cat(paste('Iteration:', 0))
#     cat('\n')
#     cat(paste('Value:', signif(f_values_list[[1]], 8)))
#     cat('\n')
#     cat(paste('Initial values: (', paste0(signif(init, 4), collapse = ', '), ')', sep = ''))
#     cat('\n')
#     cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
#     cat('\n')
#     cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
#     cat('\n')
#     cat('\n')
#   }
#   
#   for (i in 1:n) {
#     xlo <- xhi <- init
#     xlo[i] <- round_to_step(xlo[i] - delta)  # Perturbation by -delta, then round
#     flo <- f(xlo, ...)
#     xhi[i] <- round_to_step(xhi[i] + delta)  # Perturbation by delta, then round
#     fhi <- f(xhi, ...)
#     if (flo < fhi) {
#       simplex[[i + 1]] <- xlo
#       f_values_list[[i + 1]] <- flo
#     } else {
#       simplex[[i + 1]] <- xhi
#       f_values_list[[i + 1]] <- fhi
#     }
#   }
#   
#   # Evaluate function values
#   f_values <- sapply(f_values_list, as.vector)
#   
#   iter <- 0
#   while (iter < max_iter) {
#     iter <- iter + 1
#     
#     # Order simplex points by function values
#     order_idx <- order(f_values)
#     simplex <- simplex[order_idx]
#     f_values <- f_values[order_idx]
#     f_values_list <- f_values_list[order_idx]
#     b0 <- attr(f_values_list[[1]], 'beta')
#     for (i in seq_along(simplex)) attr(simplex[[i]], 'beta') <- b0
#     
#     # Centroid of all points except the worst
#     centroid <- rowMeans(do.call(cbind, simplex[1:n]))
#     
#     # Round centroid to nearest multiple of step_size
#     centroid <- round_to_step(centroid)
#     
#     if (trace > 0) {
#       cat(paste('Iteration:', iter))
#       cat('\n')
#       cat(paste('Value range: (', paste0(signif(range(f_values), 8), collapse = ', '), ')', sep = ''))
#       cat('\n')
#       cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
#       cat('\n')
#       cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
#       cat('\n')
#       cat(paste('Centroid: (', paste0(signif(simplex[[1]], 4), collapse = ', '), ')', sep = ''))
#       cat('\n')
#       cat('\n')
#     }
#     
#     # Reflection - with rounding
#     x_reflect <- round_to_step(centroid + alpha * (centroid - simplex[[n + 1]]))
#     f_reflect <- f(x_reflect, ...)
#     
#     if (f_reflect < f_values[1]) {
#       # Expansion - with rounding
#       x_expand <- round_to_step(centroid + gamma * (x_reflect - centroid))
#       f_expand <- f(x_expand, ...)
#       if (f_expand < f_reflect) {
#         simplex[[n + 1]] <- x_expand
#         f_values[n + 1] <- f_expand
#         f_values_list[[n + 1]] <- f_expand
#       } else {
#         simplex[[n + 1]] <- x_reflect
#         f_values[n + 1] <- f_reflect
#         f_values_list[[n + 1]] <- f_reflect
#       }
#     } else if (f_reflect < f_values[n]) {
#       # Accept reflection
#       simplex[[n + 1]] <- x_reflect
#       f_values[n + 1] <- f_reflect
#       f_values_list[[n + 1]] <- f_reflect
#     } else {
#       # Contraction - with rounding
#       if (f_reflect < f_values[n + 1]) {
#         # Outside contraction
#         x_contract <- round_to_step(centroid + rho * (x_reflect - centroid))
#       } else {
#         # Inside contraction
#         x_contract <- round_to_step(centroid + rho * (simplex[[n + 1]] - centroid))
#       }
#       f_contract <- f(x_contract, ...)
#       
#       if (f_contract < f_values[n + 1]) {
#         simplex[[n + 1]] <- x_contract
#         f_values[n + 1] <- f_contract
#         f_values_list[[n + 1]] <- f_contract
#       } else {
#         # Shrink the simplex - with rounding
#         for (i in 2:(n + 1)) {
#           simplex[[i]] <- round_to_step(simplex[[1]] + sigma * (simplex[[i]] - simplex[[1]]))
#           f_values_list[[i]] <- f(simplex[[i]], ...)
#           f_values[i] <- f_values_list[[i]]
#         }
#       }
#     }
#     
#     # Check convergence
#     if ((f_values[n] - f_values[1]) / abs(f_values[1]) < ftol) {
#       break
#     }
#     
#     # For discrete parameters, we need to check if all simplex points have converged
#     # to the same discrete values or are within the step_size
#     simplex_range <- mean(sapply(seq_along(init), function(i) 
#       diff(range(sapply(simplex, '[', i)))))
#     
#     if (simplex_range < step_size) {
#       break
#     }
#   }
#   
#   f1 <- f_values[1]
#   attr(f1, 'beta') <- attr(simplex[[1]], 'beta')
#   out <- list(par = as.vector(simplex[[1]]))
#   out$objective <- as.vector(f1)
#   out$iterations <- iter
#   out$beta <- attr(simplex[[1]], 'beta')
#   out$gradient <- attributes(f_values_list[[1]])$gradient
#   out$precondHessian <- attributes(f_values_list[[1]])$precondHessian
#   out$Hessian <- attributes(f_values_list[[1]])$Hessian
#   out$cholprecondHessian <- attributes(f_values_list[[1]])$cholprecondHessian
#   out$diagHessian <- attributes(f_values_list[[1]])$diagHessian
#   out$idiagHessian <- attributes(f_values_list[[1]])$idiagHessian
#   out
# }

.nelder_mead_discrete_list <- function(init, f, xtol = 1e-3, ftol = 1e-4, max_iter = 500, trace = 1, step_size = 0.2, ...) {
  n <- length(init)  # Dimension of the problem
  alpha <- 1         # Reflection coefficient
  gamma <- 2         # Expansion coefficient
  rho <- 0.5         # Contraction coefficient
  sigma <- 0.5       # Shrink coefficient
  
  # Adjust initial simplex perturbation based on step_size
  # Using step_size as minimum perturbation to ensure proper exploration
  delta <- max(step_size, 1)
  
  # Function to round to nearest multiple of step_size
  round_to_step <- function(x) {
    return(round(x / step_size) * step_size)
  }
  
  # Round initial point to nearest multiple of step_size
  init <- round_to_step(init)
  
  # Initialize the simplex
  simplex <- list()
  simplex[[1]] <- init
  f_values_list <- list()
  f_values_list[[1]] <- f(init, ...)
  
  attr(init, 'beta') <- attr(f_values_list[[1]], 'beta')
  
  if (trace > 0) {
    cat(paste('Iteration:', 0))
    cat('\n')
    cat(paste('Value:', signif(f_values_list[[1]], 8)))
    cat('\n')
    cat(paste('Initial values: (', paste0(signif(init, 4), collapse = ', '), ')', sep = ''))
    cat('\n')
    cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
    cat('\n')
    cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
    cat('\n')
    cat('\n')
  }
  
  # Create initial simplex with meaningful perturbations
  for (i in 1:n) {
    # Create two candidate points
    xlo <- xhi <- init
    
    # Ensure perturbation is at least one step size
    xlo[i] <- round_to_step(xlo[i] - delta)
    xhi[i] <- round_to_step(xhi[i] + delta)
    
    # If rounding caused no change, force a move of at least one step
    if (all(xlo == init)) xlo[i] <- init[i] - step_size
    if (all(xhi == init)) xhi[i] <- init[i] + step_size
    
    # Evaluate both points
    flo <- f(xlo, ...)
    fhi <- f(xhi, ...)
    
    # Choose the better point for the simplex
    if (flo < fhi) {
      simplex[[i + 1]] <- xlo
      f_values_list[[i + 1]] <- flo
    } else {
      simplex[[i + 1]] <- xhi
      f_values_list[[i + 1]] <- fhi
    }
  }
  
  # Evaluate function values
  f_values <- sapply(f_values_list, as.vector)
  
  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1
    
    # Order simplex points by function values
    order_idx <- order(f_values)
    simplex <- simplex[order_idx]
    f_values <- f_values[order_idx]
    f_values_list <- f_values_list[order_idx]
    b0 <- attr(f_values_list[[1]], 'beta')
    for (i in seq_along(simplex)) attr(simplex[[i]], 'beta') <- b0
    
    # Centroid of all points except the worst
    centroid <- rowMeans(do.call(cbind, simplex[1:n]))
    
    # Round centroid to nearest multiple of step_size
    centroid <- round_to_step(centroid)
    
    if (trace > 0) {
      cat(paste('Iteration:', iter))
      cat('\n')
      cat(paste('Value range: (', paste0(signif(range(f_values), 8), collapse = ', '), ')', sep = ''))
      cat('\n')
      cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
      cat('\n')
      cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
      cat('\n')
      cat(paste('Centroid: (', paste0(signif(centroid, 4), collapse = ', '), ')', sep = ''))
      cat('\n')
      cat('\n')
    }
    
    # Reflection - with rounding
    x_reflect <- round_to_step(centroid + alpha * (centroid - simplex[[n + 1]]))
    
    # Ensure the reflected point is different from centroid
    if (all(x_reflect == centroid)) {
      # Force a move in the direction of reflection
      diff_vector <- centroid - simplex[[n + 1]]
      if (all(diff_vector == 0)) {
        # If direction is undefined, choose a random direction
        diff_vector <- runif(n, -1, 1)
      }
      # Normalize and scale by step_size
      diff_vector <- diff_vector / sqrt(sum(diff_vector^2)) * step_size
      x_reflect <- round_to_step(centroid + diff_vector)
      attr(x_reflect, 'beta') <- b0
    }
    
    f_reflect <- f(x_reflect, ...)
    
    if (f_reflect < f_values[1]) {
      # Expansion - with rounding
      x_expand <- round_to_step(centroid + gamma * (x_reflect - centroid))
      
      # Ensure expanded point is different from reflected point
      if (all(x_expand == x_reflect)) {
        diff_vector <- x_reflect - centroid
        if (all(diff_vector == 0)) {
          diff_vector <- runif(n, -1, 1)
        }
        diff_vector <- diff_vector / sqrt(sum(diff_vector^2)) * step_size
        x_expand <- round_to_step(x_reflect + diff_vector)
        attr(x_expand, 'beta') <- b0
      }
      
      f_expand <- f(x_expand, ...)
      if (f_expand < f_reflect) {
        simplex[[n + 1]] <- x_expand
        f_values[n + 1] <- f_expand
        f_values_list[[n + 1]] <- f_expand
      } else {
        simplex[[n + 1]] <- x_reflect
        f_values[n + 1] <- f_reflect
        f_values_list[[n + 1]] <- f_reflect
      }
    } else if (f_reflect < f_values[n]) {
      # Accept reflection
      simplex[[n + 1]] <- x_reflect
      f_values[n + 1] <- f_reflect
      f_values_list[[n + 1]] <- f_reflect
    } else {
      # Contraction - with rounding
      if (f_reflect < f_values[n + 1]) {
        # Outside contraction
        x_contract <- round_to_step(centroid + rho * (x_reflect - centroid))
      } else {
        # Inside contraction
        x_contract <- round_to_step(centroid + rho * (simplex[[n + 1]] - centroid))
      }
      
      # Ensure contracted point is different from centroid
      if (all(x_contract == centroid)) {
        diff_vector <- if (f_reflect < f_values[n + 1]) (x_reflect - centroid) else (simplex[[n + 1]] - centroid)
        if (all(diff_vector == 0)) {
          diff_vector <- runif(n, -1, 1)
        }
        diff_vector <- diff_vector / sqrt(sum(diff_vector^2)) * step_size
        x_contract <- round_to_step(centroid + diff_vector)
        attr(x_contract, 'beta') <- b0
      }
      
      f_contract <- f(x_contract, ...)
      
      if (f_contract < min(f_values[n + 1], f_reflect)) {
        simplex[[n + 1]] <- x_contract
        f_values[n + 1] <- f_contract
        f_values_list[[n + 1]] <- f_contract
      } else {
        # Shrink the simplex - with rounding
        for (i in 2:(n + 1)) {
          simplex[[i]] <- round_to_step(simplex[[1]] + sigma * (simplex[[i]] - simplex[[1]]))
          
          # Ensure shrinkage produces different points
          if (all(simplex[[i]] == simplex[[1]])) {
            diff_vector <- simplex[[i]] - simplex[[1]]
            if (all(diff_vector == 0)) {
              # If vertices are identical, perturb randomly
              diff_direction <- runif(n, -1, 1)
              diff_direction <- diff_direction / sqrt(sum(diff_direction^2))
              simplex[[i]] <- round_to_step(simplex[[1]] + step_size * diff_direction)
            } else {
              # Maintain the original direction but ensure at least one step_size difference
              diff_direction <- diff_vector / sqrt(sum(diff_vector^2))
              simplex[[i]] <- round_to_step(simplex[[1]] + step_size * diff_direction)
            }
          }
          
          f_values_list[[i]] <- f(simplex[[i]], ...)
          f_values[i] <- f_values_list[[i]]
        }
      }
    }
    
    # Check convergence based on function values
    if ((f_values[n] - f_values[1]) / (abs(f_values[1]) + ftol) < ftol) {
      break
    }
    
    # Check convergence based on simplex geometry
    # For discrete parameters, check if all simplex points have converged
    unique_points <- unique(do.call(rbind, simplex))
    if (nrow(unique_points) == 1) {
      # All points are identical
      break
    }
    
    # Also check if the simplex has become too small
    simplex_range <- mean(sapply(seq_along(init), function(i) 
      diff(range(sapply(simplex, '[', i)))))
    
    if (simplex_range < step_size) {
      break
    }
  }
  
  f1 <- f_values[1]
  attr(f1, 'beta') <- attr(simplex[[1]], 'beta')
  out <- list(par = as.vector(simplex[[1]]))
  out$objective <- as.vector(f1)
  out$iterations <- iter
  out$beta <- attr(simplex[[1]], 'beta')
  out$gradient <- attributes(f_values_list[[1]])$gradient
  out$precondHessian <- attributes(f_values_list[[1]])$precondHessian
  out$Hessian <- attributes(f_values_list[[1]])$Hessian
  out$cholprecondHessian <- attributes(f_values_list[[1]])$cholprecondHessian
  out$diagHessian <- attributes(f_values_list[[1]])$diagHessian
  out$idiagHessian <- attributes(f_values_list[[1]])$idiagHessian
  out
}