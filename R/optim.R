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
    x0 <- x
    x <- round(x / step_size) * step_size
    attributes(x) <- attributes(x0)
    x
  }
  
  # Round initial point to nearest multiple of step_size
  init <- round_to_step(init)
  
  # Initialize the simplex
  simplex <- list()
  simplex[[1]] <- init
  f_values_list <- list()
  f_values_list[[1]] <- f(init, ...)
  
  attr(init, 'beta') <- attr(f_values_list[[1]], 'beta')
  attr(init, 'betal') <- attr(f_values_list[[1]], 'betal')
  
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
    bl <- attr(f_values_list[[1]], 'betal')
    for (i in seq_along(simplex)) {
      attr(simplex[[i]], 'beta') <- b0
      attr(simplex[[i]], 'betal') <- bl
    }
    
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
      attr(x_reflect, 'betal') <- bl
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
        attr(x_expand, 'betal') <- bl
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
        attr(x_contract, 'betal') <- bl
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
  attr(f1, 'betal') <- attr(simplex[[1]], 'betal')
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

