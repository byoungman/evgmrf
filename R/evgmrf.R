#' Fitting extreme value distributions with parameters that vary according to Gaussian Markov random fields.
#'
#' Function \code{evgmrf} fits extreme value distributions with one or more parameters
#' varying with Gaussian Markov random field form.
#'
#' @param z an array, or matrix if W supplied, or array if W supplied and r-largest
#' @param formula to do
#' @param covariates to do
#' @param family to do
#' @param weights to do
#' @param inits to do
#' @param model to do
#' @param W adjacency matrix
#' @param rho0 to do
#' @param alpha to do
#' @param order to do
#' @param lambda0 (log of) initial smoothing parameter values; repeated if necessary
#' @param pp_nper to do
#' @param C to do
#' @param tau to do
#' @param rank to do
#' @param args to do
#' @param control to do
#' @param trace an integer for the level of fitting detail to report. Larger values give more detail. Defaults to 0 (none)
#' @param outer to do
#' @param gamma a scalar that multiplies the log-likelihood (similarly to a constant weights); defaults to 1
#' @param nx,ny grid dimensions, if not deducible from z
#' @param index a matrix of row and column indices, to define grid-based data
#' @param infill infill grid holes values during inference?
#' @param cv use cross-validation mode?
#' @param inits.method to do
#' @param auto.weights to do
#' 
#' @details
#' 
#' To do.
#' 
#' @references 
#' 
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Models. Journal of Statistical Software. \doi{10.18637/jss.v103.i03}
#'
#' @examples
#'
#' # To follow
#'
#' @seealso \link{predict.evgmrf}
#'
#' @return An object of class \code{evgmrf}
#' 
#' @export
#' 
evgmrf <- function(z, formula = ~ -1, covariates, family = 'gev', weights = 1, inits = 'different',
                   model = 'ICAR', W = NULL, rho0 = 1, order = 1, lambda0, # GMRF -1 assumes GMRF for each par
                   alpha = NA,
                   pp_nper,                    # point process
                   rank = -1,
                   args = list(),             # eigen rank 
                   control = list(), # some control parameters for testing
                   trace = 0,
                   outer = 'newton',
                   gamma = 1,
                   nx = NULL, ny = NULL, index = NULL, infill = FALSE,
                   cv = FALSE, 
                   inits.method = 'reml',
                   auto.weights = FALSE
) {
  model <- tolower(model)
  args <- replace(.args0, names(args), args)
  if (family == 'pproc') {
    z <- array(z, c(1, dim(z)))
  }
  w <- 0 * z + weights
  control <- replace(.control.evgmrf(), names(control), control)
  .checks(model, order)
  # some basics
  dz <- dim(z)
  m <- dz[1]
  rlargeish <- substr(family, 1, 6) == 'rlarge'
  if (rlargeish) {
    r <- dz[2]
  }
  holes <- FALSE
  if (infill) {
    if (is.null(index))
      stop("Argument 'index' must be supplied if 'infill = TRUE'.")
    nx <- max(index[, 1])
    ny <- max(index[, 2])
    n <- nx * ny
    here <- matrix(FALSE, nx, ny)
    here[index] <- TRUE
    dz0 <- c(dz[- length(dz)], n)
    dz <- c(dz[- length(dz)], nx, ny)
    zz <- array(NA, dz0)
    zz[ ,, as.logical(here)] <- z
    z <- array(zz, dz)
    index <- NULL
  }  
  # check whether in grid form
  if (length(dz) == 2) {
    # matrix-form
    if (rlargeish)
      stop("z can't be a matrix for r-largest family")
    nx <- nx
    ny <- ny
    zm <- matrix(z, m)
    wm <- matrix(w, m)
    n <- ncol(zm)
    zl <- lapply(1:n, function(i) zm[, i])
    wl <- lapply(1:n, function(i) wm[, i])
    holes <- TRUE
  } else {
    if (rlargeish) {
      if (length(dz) == 4) {
        nx <- dz[3]
        ny <- dz[4]
        n <- nx * ny
        za <- array(z, c(m, r, n))
        wa <- array(w, c(m, r, n))
        zl <- lapply(1:n, function(i) as.matrix(za[ , , i]))
        wl <- lapply(1:n, function(i) as.matrix(wa[ , , i]))
      } else {
        n <- dz[3]
        zl <- lapply(1:n, function(i) as.matrix(z[ , , i]))
        wl <- lapply(1:n, function(i) as.matrix(w[ , , i]))
        holes <- TRUE
      }
    } else {
      nx <- dz[2]
      ny <- dz[3]
      n <- nx * ny
      zm <- matrix(z, m)
      wm <- matrix(w, m)
      zl <- lapply(1:n, function(i) zm[, i])
      wl <- lapply(1:n, function(i) wm[, i])
    }
  }
  if (family == 'pp') {
    wl <- lapply(seq_along(zl), function(i) wl[[i]][order(zl[[i]], decreasing = TRUE, na.last = NA)])
    zl <- lapply(seq_along(zl), function(i) zl[[i]][order(zl[[i]], decreasing = TRUE, na.last = NA)])
    if (is.null(args$u)) {
      if (is.null(args$r)) {
        stop("Supply either args$r or args$u for family = 'pp'.")
      } else {
        args$u <- sapply(zl, function(x) x[args$r])
      }
    }
    u <- as.vector(args$u)
    if (length(u) != length(zl))
      stop("args$u is incompatible with z.")
    # zl <- mapply(function(x, y) x[x > y], zl, u)
  }
  .ld <- list(n = n, z = zl, w = wl, mult = 1 / gamma, args = args, bym4 = FALSE)
  # .ld$z_cube <- array(unlist(.ld$z), c(dim(.ld$z[[1]])[1:2], length(.ld$z)))
  if (family == 'gpd') {
    .lf <- .gpd_fns
    .ld$np <- 2
    if (inits == 'same') {
      p0m1 <- .quick_tgpd(na.omit(unlist(zl)))
      p0m <- t(matrix(p0m1, 2, n))
    } else {
      p0m <- t(sapply(zl, .quick_tgpd))
    }
  }
  if (family == 'pp') {
    .lf <- .pp_fns
    .ld$m <- args$nper
    .ld$np <- 3
    .ld$u <- u
    .ld$uw <- 0 * u + args$nper
    if (is.null(args$u)) {
      .ld$u <- sapply(.ld$z, min, na.rm = TRUE)
    } else {
      .ld$u <- as.vector(args$u)
    }
    if (inits == 'same') {
      p0m1 <- .quick_tpp(unlist(zl), .ld$m, min(.ld$u), args$delta)
      p0m <- t(matrix(p0m1, 3, n))
    } else {
      p0m <- t(mapply(.quick_tpp, zl, .ld$u, m = .ld$m, delta = args$delta))
    }
  }
  if (family == 'gev') {
    .lf <- .gev_fns
    .ld$np <- 3
    p0m1 <- .quick_tgev(na.omit(unlist(zl)), args$delta)
    # p0m <- t(matrix(p0m1, 3, n))
    if (inits == 'different') {
      p0m <- t(sapply(zl, .quick_tgev_shrink, delta = args$delta, pars0 = p0m1, mult = args$mult))
      set_to_mean <- is.na(p0m[, 1])
      if (any(set_to_mean)) {
        infill <- matrix(p0m1, sum(set_to_mean), 3, byrow = TRUE)
        p0m[set_to_mean, ] <- infill
      }
    } else {
      p0m <- t(matrix(p0m1, 3, n))
    }
    if (auto.weights) {
      p1 <- p0m[, 1]
      p2 <- exp(p0m[, 2])
      p3 <- 1.5 / (1 + exp(-p0m[, 3])) - 1
      #zm_exp <- lapply(1:n, function(i) -log(exp(-(1 + p3[i] * (zl[[i]] - p1[i]) / p2[i])^(-1/p3[i]))))
      #z0l_exp2 <- lapply(1:n, function(i) (1 + p3[i] * (zl[[i]] - p1[i]) / p2[i])^(-1/p3[i]))
      zm_exp <- (1 + p3 * (t(zm) - p1) / p2)^(-1/p3)
      zm_exp_agg <- matrix(zm_exp, ncol = 100)
      wt <- (1 / mean(apply(zm_exp_agg, 2, min))) / nrow(zm_exp_agg)
      # wt <- (1 / mean(apply(zm_exp, 2, min))) / n#1 / (n * mean(apply(zm_exp, 1, min)))
      .ld$w <- lapply(.ld$w, function(x) 0 * x + wt)
    }
    # if (inits == 'same') {
    #   p0m1 <- .quick_tgev(na.omit(unlist(zl)), args$delta)
    #   p0m <- t(matrix(p0m1, 3, n))
    # } else {
    #   p0m <- t(sapply(zl, .quick_tgev, delta = args$delta))
    #   set_to_mean <- is.na(p0m[, 1])
    #   if (any(set_to_mean)) {
    #     infill <- matrix(.quick_tgev(unlist(zl), delta = args$delta), sum(set_to_mean), 3, byrow = TRUE)
    #     p0m[set_to_mean, ] <- infill
    #   }
    # }
  }
  if (family == 'rlarge') {
    # if (is.null(args$drop)) {
      .lf <- .rlarge_fns
      # p0m1 <- .quick_tgev(na.omit(unlist(lapply(zl, function(x) x[, 1]))), args$delta)
      p0m1 <- .quick_tgev(na.omit(unlist(lapply(zl, function(x) x))), args$delta)
      p0m <- t(matrix(p0m1, 3, n))
    if (inits == 'different') {
      # p0m <- t(sapply(lapply(zl, function(x) x[, 1]), .quick_tgev_shrink, delta = args$delta, pars0 = p0m1, mult = args$mult))
      p0m <- t(sapply(lapply(zl, function(x) x), .quick_tgev_shrink, delta = args$delta, pars0 = p0m1, mult = args$mult))
      set_to_mean <- is.na(p0m[, 1])
      if (any(set_to_mean)) {
        infill <- matrix(.quick_tgev(unlist(lapply(zl, function(x) x[, 1])), delta = args$delta), sum(set_to_mean), 3, byrow = TRUE)
        p0m[set_to_mean, ] <- infill
      }
    }
    # } else {
    #   .lf <- .rlargec_fns
    #   .ld$drop <- args$drop
    #   if (inits == 'same') {
    #     p0m1 <- .quick_tgev(na.omit(unlist(lapply(zl, function(x) x[, 1]))), delta = args$delta)
    #     p0m <- t(matrix(p0m1, 3, n))
    #   } else {
    #     p0m <- t(sapply(zl, function(x) .quick_tgev(x, delta = args$delta)))
    #     set_to_mean <- is.na(p0m[, 1])
    #     if (any(set_to_mean)) {
    #       infill <- matrix(.quick_tgev(unlist(zl)), sum(set_to_mean), 3, byrow = TRUE)
    #       p0m[set_to_mean, ] <- infill
    #     }
    #   }
    # }
    .ld$np <- 3
  }
  if (family == 'ald') {
    .lf <- .ald_fns
    .ld$np <- 2
    if (is.null(args$tau))
      stop("Must supply args$tau for family = 'ald'.")
    p0m <- t(sapply(zl, .quick_ald, args = args))
  }
  if (family == 'pproc') {
    .lf <- .pproc_fns
    .ld$m <- ifelse(is.null(args$nper), 1, args$nper)
    .ld$mult <- .ld$m# / n
    .ld$np <- 1
    p0m <- .quick_pproc(unlist(zl))
    p0m <- p0m - log(.ld$m)
    p0m <- matrix(p0m, n, 1)
  }
  .ld$np0 <- .ld$np
  if (.ld$np > 1) {
    if (length(model) == 1)
      model <- rep(model, .ld$np)
    if (length(order) == 1)
      order <- rep(order, .ld$np)
    if (length(alpha) == 1)
      alpha <- rep(alpha, .ld$np)
  }  
  .ld$np <- sum(!is.na(model))
  # put together fixed effects stuff
  if (class(formula) == 'formula')
    formula <- lapply(1:.ld$np0, function(i) formula)
  if (missing(covariates)) {
    covariates <- data.frame(id = 1:n)
    if (!is.null(nx) & !is.null(ny)) {
      xy <- expand.grid(x = seq_len(nx), y = seq_len(ny))
      covariates$x <- xy$x
      covariates$y <- xy$y
    }
  } else {
    covariates <- as.data.frame(lapply(covariates, as.vector))
  }
  gmrf <- !is.na(model)
  if (sum(gmrf) == 0)
    stop('Model must have at least one GMRF prior.')
  if (any(!gmrf)) {
    for (i in which(!gmrf)) {
      if (formula[[i]] == ~-1)
        formula[[i]] <- ~ 1
    }
  }
  X1 <- lapply(formula, model.matrix, data = covariates)
  nX1 <- sapply(X1, ncol)
  p0l <- list()
  for (i in 1:.ld$np0) {
    if (nX1[[i]] > 0) {
      b0 <- solve(crossprod(X1[[i]]), crossprod(X1[[i]], p0m[, i]))
      b1 <- as.vector(X1[[i]] %*% b0)
      p0m[, i] <- p0m[, i] - b1
    } else {
      b0 <- numeric(0)
    }
    p0l[[i]] <- b0
  }
  id_bym2 <- lapply(p0l, function(x) rep(FALSE, length(x)))
  for (i in which(gmrf)) {
    p0l[[i]] <- c(p0l[[i]], p0m[, i])
    id_bym2[[i]] <- c(id_bym2[[i]], rep(FALSE, length(p0m[, i])))
    if (model[i] %in% paste('bym', 2:4, sep = '')) {
      id_bym2[[i]] <- rep(c(FALSE, TRUE), c(.ld$n, length(p0l[[i]])))
      p0l[[i]] <- c(numeric(.ld$n), p0l[[i]])
    }
  }
  .ld$psplit <- rep(seq_along(p0l), sapply(p0l, length))
  ## put together GMRF stuff
  if (!is.null(W)) {
    if (nrow(W) != nrow(p0m)) 
      stop(paste('Supplied W dimension not compatible with data size:', nrow(W), '!=', nrow(p0m)))
  }
  Qd0 <- .makeQ_data(nx, ny, model, order, alpha, nX1, W)
  .ld$Xl <- X1
  for (i in which(gmrf)) {
    .ld$Xl[[i]] <- list(.ld$Xl[[i]], Matrix::Diagonal(.ld$n))
    if (model[i] %in% paste('bym', 2:4, sep = ''))
      .ld$Xl[[i]] <- c(.ld$Xl[[i]], Matrix::Diagonal(.ld$n))
  }
  .ld$Xlc <- .ld$Xl # componentwise
  .ld$Xl <- lapply(seq_along(.ld$Xl), function(i) do.call(cbind, .ld$Xl[[i]]))
  Qd0$R <- list()
  for (i in which(gmrf))
    Qd0$R[[i]] <- Matrix::qrR(Matrix::qr(rbind(.ld$Xl[[i]], .ld$Xl[[i]])))
  .ld$X <- Matrix::.bdiag(.ld$Xl)
  .ld$X0 <- Matrix::Diagonal(.ld$n)
  .ld$id_bym2 <- id_bym2
  .ld$openmp <- control$openmp
  .ld$threads <- control$threads
  .ld$control <- control
  # rho0_temp <- as.vector(unlist(lapply(model, .inits_model)))
  # par_var <- -log(20 * apply(p0m, 2, var))
  par_var <- log(apply(p0m, 2, sd))
  lambda0_temp <- unlist(mapply(.inits_model, model, order, alpha, par_var, USE.NAMES = FALSE, SIMPLIFY = FALSE))
  if (missing(lambda0)) {
    lambda0 <- lambda0_temp
  } else {
    if (length(lambda0) < length(lambda0_temp)) {
      lambda0 <- rep(lambda0, length(model))
    }
    if (length(lambda0) != length(lambda0_temp))
      stop('Invalid rho0 supplied')
    names(lambda0) <- names(lambda0_temp)
  }
  p0v <- as.vector(p0m)
  p0v <- unlist(p0l)
  attr(lambda0, 'beta') <- p0v
  attr(lambda0, 'first') <- TRUE
  # if (control$sandwich) {
  #   p0v <- .newton(p0v, .d0_Q, .search_Q, likdata = .ld, likfns = .lf, Q = .mQ(lambda0, Qd0))
  #   G <- .G_Q(p0v$par, likdata = .ld, likfns = .lf, Q = .mQ(lambda0, Qd0))
  #   J <- tcrossprod(G)
  #   # J <- .J_Q(p0v$par, likdata = .ld, likfns = .lf, Q = .mQ(lambda0, Qd0))
  #   # H <- p0v$precondHessian
  #   iHJ <- .cholsolveAB(p0v$Hessian, as.matrix(J))
  #   # 
  #   # lambda <- eigen(iHJ, only.values = TRUE)$values
  #   # return(1 / mean(Re(lambda)))
  #   wt <- mean(diag(iHJ))
  #   .ld$w <- lapply(.ld$w, function(x) x * wt)
  #   p0v <- p0v$par
  #   # print(.ld$w)
  # }
  if ('diag' %in% inits.method)
     p0v <- .newton_diag(.d0_Q, .d1_Q, .d2_Q_diag, p0v, max_iter=1e4, likdata = .ld, likfns = .lf, Q = .mQ(lambda0, Qd0))
  # p0v <- .beta0(rho0, Qd0, .ld, .lf, .mQ, control$it0)$par
  attr(lambda0, 'beta') <- p0v
  # .ld$L0 <- .chol0(p0v, .ld, .lf, .mQ(rho0, Qd0))
  if (control$inner_optim == 'Cholesky')
    .ld$chol0 <- .Cholesky0(lambda0, Qd0, control$super)
  if (is.list(cv)) {
    valid <- cv$valid
    .ldv <- .ldv2 <- lapply(seq_along(valid), function(.) .ld)
    for (j in seq_along(valid)) {
      for (i in valid[[j]]) {
        .ldv[[j]]$z[[i]] <- rep(NA, length(.ldv[[j]]$z[[i]]))
        .ldv[[j]]$w[[i]] <- rep(NA, length(.ldv[[j]]$w[[i]]))
      }
      for (i in seq_along(.ldv[[j]]$z)[-valid[[j]]]) {
        .ldv2[[j]]$z[[i]] <- rep(NA, length(.ldv2[[j]]$z[[i]]))
        .ldv2[[j]]$w[[i]] <- rep(NA, length(.ldv2[[j]]$w[[i]]))
      }
    }
    .ldv$valid <- .ldv2
    out <- .nelder_mead_discrete_list(lambda0, .cv0, likdata = .ldv, likfns = .lf, Qd = Qd0, makeQ = .mQ,
                                     trace = trace, step_size = control$step_size)
    # out <- .BFGS(lambda0, .cv0, .cv1, likdata = .ldv, likfns = .lf, Qd = Qd0, makeQ = .mQ, 
    #             eps = control$cv_eps, direction = control$reml_direction, trace = trace > 0, gradtol = 1e-2, 
    #             stepmax = control$reml_stepmax, steptol = control$cv_steptol)
    # out <- .newton(lambda0, .cv0, .cv_step, likdata = .ldv, likfns = .lf, Qd = Qd0, makeQ = .mQ, reml = TRUE,
    #                eps = control$cv_eps, direction = control$reml_direction, trace = trace > 0, gradtol = control$cv_gradtol, 
    #                stepmax = control$reml_stepmax, steptol = control$reml_steptol, itlim = control$reml_itlim)
  } else {
  if (outer == 'nelder-mead') {
    if (length(lambda0) > 1) {
      # out <- .nelder_mead_list(rho0, .reml0, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ, trace = trace)
      out <- .nelder_mead_discrete_list(lambda0, .reml0, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ, 
                                        trace = trace, step_size = control$step_size)
    } else {
      lambda0[] <- 0
      out <- .brent(lambda0, .reml0, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ, trace = trace)
    }
  } else {
    if (outer == 'newton') {
      out <- .newton(lambda0, .reml0, .reml_step, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ, reml = TRUE,
                    eps = control$reml_eps, direction = control$reml_direction, trace = trace > 0, gradtol = 1e-2, 
                    stepmax = control$reml_stepmax, steptol = control$reml_steptol, itlim = control$reml_itlim)
    } else {
      out <- .BFGS(lambda0, .reml0, .reml1, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ, 
                   eps = control$reml_eps, direction = control$reml_direction, trace = trace > 0, gradtol = 1e-2, 
                   stepmax = control$reml_stepmax, steptol = control$reml_steptol)
    }
  }}
  out$likdata <- .ld
  out$likfns <- .lf
  out$family <- family
  if (is.null(index)) {
    if (!is.null(nx) & !is.null(ny))
      index <- as.matrix(expand.grid(row = 1:nx, col = 1:ny))
  } else {
    nx <- max(index[, 1])
    ny <- max(index[, 2])
    holes <- nrow(index) < nx * ny
  }
  if (control$inner_optim != 'Cholesky')
    out$cholprecondHessian <- Matrix::Cholesky(out$precondHessian, super = control$super, LDL = FALSE)
  out$X <- .ld$Xl
  out$holes <- holes
  out$nx <- nx
  out$ny <- ny
  out$n <- .ld$n
  out$np <- .ld$np
  out$unlink <- .lf$trans
  out$names <- .lf$names
  out$quantile <- .lf$quantile
  out$quantile0 <- .lf$quantile0
  out$gmrf <- gmrf
  out$model <- model
  attr(out$beta, 'split') <- .ld$psplit
  out$index <- index
  out$init <- p0m
  if (control$sandwich) {
    names(out$pars) <- names(lambda0)
    # G <- .G_Q(out$beta, likdata = .ld, likfns = .lf, Q = .mQ(out$pars, Qd0))
    # out$J <- tcrossprod(G)
    G0 <- .G_Q0(out$beta, likdata = .ld, likfns = .lf)
    G_centered <- G0 - rowMeans(G0)
    out$Gc <- G_centered
    out$J <- tcrossprod(G_centered)
  }
  if (auto.weights)
    out$awt <- wt
  class(out) <- 'evgmrf'
  out
}

