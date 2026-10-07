## -----------------------------------------------------------------------------
## Internal helpers
## -----------------------------------------------------------------------------

#' Arrange response data and weights onto the spatial grid
#'
#' Converts `z` (array, compacted matrix/array, or list) into a per-location
#' list, works out the grid dimensions and which cells hold data.
#'
#' @return A list with `z` (array with grid as last two dimensions), `zl`, `wl`
#'   (per-location lists), `n`, `nx`, `ny`, `index`, `here` and `holes`.
#' @noRd
.evgmrf_data <- function(z, family, weights, nx = NULL, ny = NULL,
                         index = NULL, W = NULL) {
  rlargeish <- substr(family, 1, 6) == 'rlarge'
  z_is_list <- is.list(z)
  holes <- TRUE

  if (!z_is_list) {
    w <- 0 * z + weights
    dz <- dim(z)
    m <- dz[1]
    array_style <- length(dz) == 3 + as.integer(rlargeish)
    if (array_style) {
      # full grid: last two dimensions of z are the grid
      holes <- FALSE
      nx <- tail(dz, 2)[1]
      ny <- tail(dz, 2)[2]
      n <- nx * ny
      index <- as.matrix(expand.grid(seq_len(nx), seq_len(ny)))
      if (!rlargeish) {
        zm <- matrix(z, m)
        here <- colSums(!is.na(zm)) > 0
        wm <- matrix(w, m)
        zl <- lapply(seq_len(ncol(zm)), function(i) zm[, i])
        wl <- lapply(seq_len(ncol(zm)), function(i) wm[, i])
      } else {
        za <- array(z, c(dz[1:2], n))
        here <- colSums(!is.na(za[, 1, ])) > 0
        wa <- array(w, c(dz[1:2], n))
        zl <- lapply(seq_len(dim(za)[3]), function(i) as.matrix(za[, , i]))
        wl <- lapply(seq_len(dim(za)[3]), function(i) as.matrix(wa[, , i]))
      }
    } else {
      # compacted: last dimension indexes locations
      n <- tail(dz, 1)
      if (rlargeish) {
        if (length(dz) == 2)
          stop("z can't be a matrix for r-largest family")
        z <- lapply(seq_len(n), function(i) as.matrix(z[, , i]))
        w <- lapply(seq_len(n), function(i) as.matrix(w[, , i]))
      } else {
        z <- lapply(seq_len(n), function(i) z[, i])
        w <- lapply(seq_len(n), function(i) w[, i])
      }
    }
  }

  if (holes) {
    zl <- z
    if (length(weights) == 1) {
      wl <- lapply(zl, function(x) 0 * x + weights)
    } else {
      wl <- if (z_is_list) weights else w
    }
    n <- length(zl)
    grid <- .evgmrf_grid(n, nx, ny, index, W)
    nx <- grid$nx
    ny <- grid$ny
    index <- grid$index
    here <- matrix(FALSE, nx, ny)
    here[index] <- TRUE
    here <- as.logical(here)
    if (rlargeish) {
      z <- .list2array(zl)
      m <- dim(z)[1]
      r <- dim(z)[2]
      w <- .list2array(wl)
      zz <- ww <- array(NA, c(m, r, n))
      zz[, , here] <- z
      ww[, , here] <- w
    } else {
      z <- .list2mat(zl)
      m <- nrow(z)
      w <- .list2mat(wl)
      zz <- ww <- matrix(NA, m, n)
      zz[, here] <- z
      ww[, here] <- w
    }
    dz <- c(dim(zz)[-length(dim(zz))], nx, ny)
    z <- array(zz, dz)
    here <- rowSums(matrix(apply(is.finite(z), -1, any), n)) > 0
    holes <- any(!here)
  }

  list(z = z, zl = zl, wl = wl, n = n, nx = nx, ny = ny, index = index,
       here = here, holes = holes)
}

#' Work out grid dimensions for compacted data
#'
#' Implements the rules in the Details section of `?evgmrf`.
#' @noRd
.evgmrf_grid <- function(n, nx, ny, index, W) {
  if (!is.null(index)) {
    if (is.null(nx)) {
      nx <- max(index[, 1])
    } else if (max(index[, 1]) > nx) {
      stop('index and nx not compatible')
    }
    if (is.null(ny)) {
      ny <- max(index[, 2])
    } else if (max(index[, 2]) > ny) {
      stop('index and ny not compatible')
    }
    return(list(nx = nx, ny = ny, index = index))
  }
  if (is.null(nx) && is.null(ny) && is.null(W))
    stop('Either nx, ny, index or W need supplying.')
  if (is.null(nx) && !is.null(ny)) {
    nx <- n / ny
    if (nx - floor(nx) != 0)
      stop('ny incompatible with data length')
  }
  if (is.null(ny) && !is.null(nx)) {
    ny <- n / nx
    if (ny - floor(ny) != 0)
      stop('nx incompatible with data length')
  }
  if (!is.null(W) && is.null(nx))
    nx <- n
  if (!is.null(W) && is.null(ny))
    ny <- 1
  list(nx = nx, ny = ny,
       index = as.matrix(expand.grid(seq_len(nx), seq_len(ny))))
}

#' Threshold exceedances for the Poisson-GPD family
#'
#' Sorts each location's data, sets thresholds from `args$u` or `args$r`, and
#' keeps only exceedances.
#' @noRd
.evgmrf_poisgpd <- function(zl, wl, args) {
  ord <- lapply(zl, function(x) order(x, decreasing = TRUE, na.last = NA))
  wl <- lapply(seq_along(wl), function(i) wl[[i]][ord[[i]]])
  zl <- lapply(seq_along(zl), function(i) zl[[i]][ord[[i]]])
  if (is.null(args$u)) {
    if (is.null(args$r))
      stop("Supply either args$r or args$u for family = 'poisgpd'.")
    args$u <- sapply(zl, function(x) x[args$r])
  }
  if (length(args$u) == 1)
    args$u <- rep(args$u, length(zl))
  u <- as.vector(args$u)
  if (length(u) != length(zl))
    stop("args$u is incompatible with z.")
  excl <- lapply(seq_along(zl), function(i) zl[[i]] > args$u[i])
  zl <- lapply(seq_along(zl), function(i) zl[[i]][excl[[i]]])
  wl <- lapply(seq_along(wl), function(i) wl[[i]][excl[[i]]])
  list(zl = zl, wl = wl, args = args, u = u)
}

#' Family-specific likelihood setup and initial parameter values
#'
#' @return A list with `likdata` (the `.ld` object), `likfns` (the `.lf`
#'   object) and `inits` (an `n` by `np` matrix of initial values).
#' @noRd
.evgmrf_family <- function(family, zl, wl, n, here, args, gamma, bymfns,
                           inits.method, u = NULL) {
  known <- c('gpd', 'poisgpd', 'gev', 'rlarge', 'ald', 'pproc')
  if (!(family %in% known))
    stop(paste0("Unknown family '", family, "'. Must be one of: ",
                paste(known, collapse = ', ')))
  ld <- list(n = n, z = zl, w = wl, mult = 1, args = args, gamma = gamma)
  ld$bymfns <- bymfns
  inits_list <- list()

  if (family == 'gpd') {
    lf <- .gpd_fns
    ld$np <- 2
    inits_list$same <- .quick_tgpd(na.omit(unlist(zl)))
    if (inits.method == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tgpd(zl[[i]]))
  }
  if (family == 'poisgpd') {
    lf <- .pp_fns
    ld$m <- args$nper
    ld$np <- 3
    ld$u <- u
    ld$uw <- 0 * u + args$nper
    inits_list$same <- .quick_tpp(na.omit(unlist(zl)), ld$m, min(ld$u, na.rm = TRUE), args$delta)
    if (inits.method == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tpp(zl[[i]], ld$u[i], m = ld$m, delta = args$delta))
  }
  if (family == 'gev') {
    lf <- .gev_fns
    ld$np <- 3
    inits_list$same <- .quick_tgev(na.omit(unlist(zl)), args$delta)
    if (inits.method == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tgev_shrink(zl[[i]], delta = args$delta, pars0 = inits_list$same, mult = args$mult))
  }
  if (family == 'rlarge') {
    lf <- .rlarge_fns
    ld$np <- 3
    inits_list$same <- .quick_tgev(na.omit(unlist(zl)), args$delta)
    if (inits.method == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tgev_shrink(zl[[i]], delta = args$delta, pars0 = inits_list$same, mult = args$mult))
  }
  if (family == 'ald') {
    lf <- .ald_fns
    ld$np <- 2
    if (is.null(args$tau))
      stop("Must supply args$tau for family = 'ald'.")
    inits_list$same <- .quick_ald(na.omit(unlist(zl)), args = args)
    if (inits.method == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_ald(zl[[i]], args = args))
  }
  if (family == 'pproc') {
    lf <- .pproc_fns
    ld$m <- ifelse(is.null(args$nper), 1, args$nper)
    ld$mult <- ld$m
    ld$np <- 1
    inits_list$same <- .quick_pproc(na.omit(unlist(zl)))
    inits_list$same <- inits_list$same - log(ld$m)
    p0m <- .quick_pproc(unlist(zl))
    if (inits.method == 'different') {
      inits_list$diff <- sapply(which(here), function(i) .quick_pproc(zl[[i]]))
      inits_list$diff <- inits_list$diff - log(ld$m)
    }
  }

  if (inits.method == 'same')
    inits <- matrix(inits_list$same, ld$np, n)
  if (inits.method == 'different') {
    if (any(!here)) {
      inits <- matrix(inits_list$same, ld$np, n)
      inits[, here] <- inits_list$diff
    } else {
      inits <- inits_list$diff
    }
  }
  # keep shape starting values away from problematic region
  if (family %in% c('gev', 'poisgpd'))
    inits[3, ] <- pmax(inits[3, ], .73)
  if (family %in% c('gpd'))
    inits[2, ] <- pmax(inits[2, ], .73)

  ld$np0 <- ld$np
  list(likdata = ld, likfns = lf, inits = t(inits))
}

#' Expand hyper, order, model and formula to one entry per parameter
#' @noRd
.evgmrf_expand <- function(hyper, order, model, formula, np) {
  if (length(hyper) == 1 && (np > 1 || !is.list(hyper[[1]])))
    hyper <- lapply(seq_len(np), function(.) hyper)
  if (!is.list(order)) {
    if (is.vector(order)) {
      if (length(order) == 1)
        order <- rep(order, np)
      order <- lapply(order, seq_len)
    }
  } else {
    if (length(order) == 1)
      order <- rep(order, np)
  }
  order <- lapply(order, sort)
  if (np > 1 && length(model) == 1)
    model <- rep(model, np)
  if (is.null(formula)) {
    formula <- lapply(model, .model2formula)
  } else if (is(formula, 'formula')) {
    formula <- lapply(model, function(.) formula)
  }
  list(hyper = hyper, order = order, model = model, formula = formula)
}

#' Build and check the covariate data frame
#' @noRd
.evgmrf_covariates <- function(covariates, n, nx, ny) {
  if (is.null(covariates)) {
    covariates <- data.frame(id = seq_len(n))
    if (!is.null(nx) && !is.null(ny)) {
      xy <- expand.grid(x = seq_len(nx), y = seq_len(ny))
      covariates$x <- xy$x
      covariates$y <- xy$y
    }
    return(covariates)
  }
  bad <- paste0('Covariates ', paste(names(covariates), collapse = ' and '),
                ' not compatible with z')
  if (!is.null(nx) && !is.null(ny)) {
    cov_dims <- lapply(covariates, dim)
    compatible <- sapply(cov_dims, function(x) if (length(x) == 2) all(x == c(nx, ny)) else TRUE)
  } else {
    covariates <- lapply(covariates, as.vector)
    compatible <- sapply(covariates, length) == n
  }
  if (any(!compatible))
    stop(bad)
  as.data.frame(lapply(covariates, as.vector))
}

#' Design matrices and starting coefficient vector
#'
#' Fits fixed effects to the initial values by least squares, leaving the
#' residual as the starting GMRF field, and builds the per-parameter design
#' matrices (fixed effects | GMRF | unstructured for BYM2/3).
#' @noRd
.evgmrf_design <- function(formula, covariates, model, inits, n, np0) {
  gmrf <- !is.na(model)
  if (sum(gmrf) == 0)
    stop('Model must have at least one GMRF prior.')
  X1 <- lapply(formula, model.matrix, data = covariates)
  nX1 <- sapply(X1, ncol)
  bym23 <- paste('bym', 2:3, sep = '')

  p0l <- id_bym2 <- par_type <- list()
  for (i in seq_len(np0)) {
    if (nX1[[i]] > 0) {
      b0 <- solve(crossprod(X1[[i]]), crossprod(X1[[i]], inits[, i]))
      b1 <- as.vector(X1[[i]] %*% b0)
      inits[, i] <- inits[, i] - b1
      par_type[[i]] <- c('parametric')
    } else {
      b0 <- numeric(0)
      par_type[[i]] <- character()
    }
    p0l[[i]] <- list(b0)
    if (gmrf[i])
      p0l[[i]][[2]] <- inits[, i]
  }
  for (i in which(gmrf)) {
    id_bym2[[i]] <- rep(c(FALSE, FALSE), sapply(p0l[[i]], length))
    par_type[[i]] <- c(par_type[[i]], 'spatial')
    if (model[i] %in% bym23) {
      id_bym2[[i]] <- c(id_bym2[[i]], rep(TRUE, n))
      p0l[[i]] <- c(p0l[[i]][[1]], .1 * p0l[[i]][[2]], .9 * p0l[[i]][[2]])
      par_type[[i]] <- c(par_type[[i]], 'random')
    }
  }
  p0l <- lapply(p0l, unlist)
  psplit <- rep(seq_along(p0l), sapply(p0l, length))

  # componentwise design matrices
  Xlc <- X1
  for (i in seq_along(Xlc)) {
    if (gmrf[i]) {
      Xlc[[i]] <- list(Xlc[[i]], Matrix::Diagonal(n))
      if (model[i] %in% bym23)
        Xlc[[i]] <- c(Xlc[[i]], Matrix::Diagonal(n))
    } else {
      Xlc[[i]] <- list(Xlc[[i]])
    }
  }
  Xl <- lapply(Xlc, function(x) do.call(cbind, x))

  list(gmrf = gmrf, X1 = X1, nX1 = nX1, inits = inits, p0v = unlist(p0l),
       psplit = psplit, id_bym2 = id_bym2, par_type = par_type,
       Xl = Xl, Xlc = Xlc, X = Matrix::.bdiag(Xl))
}

#' Starting values for the GMRF hyperparameters
#' @noRd
.evgmrf_hyper0 <- function(model, order, inits, bymfns, hyper, hyper_start,
                           np, np0) {
  par_var <- apply(inits, 2, sd)
  if (is.null(bymfns))
    bymfns <- lapply(seq_along(model), function(.) NULL)
  for (i in seq_along(hyper)) {
    if (is.null(hyper[[i]]))
      hyper[[i]] <- list()
  }
  hyper0 <- lapply(seq_len(np), function(i)
    .inits_model(model[i], order[[i]], par_var[i], bymfns[[i]], hyper[[i]]))

  if (!is.null(hyper_start)) {
    if (!is.list(hyper_start[[1]]))
      hyper_start <- lapply(seq_len(np0), function(.) hyper_start)
    for (i in seq_along(hyper_start)) {
      if (!is.null(hyper_start[[i]]))
        hyper0[[i]][names(hyper_start[[i]])] <- hyper_start[[i]]
    }
  } else if ('bym3' %in% model) {
    for (i in which(model %in% 'bym3')) {
      fmls <- formals(bymfns[[i]])[-1]
      if (length(fmls) > 0) {
        fmls_need <- names(fmls)
        fmls_need <- fmls_need[!(names(fmls) %in% names(hyper[[i]]))]
        if (length(fmls_need) > 0)
          stop(paste0('Fixed or starting values missing for BYM3 model parameters ',
                      paste0(fmls_need, collapse = ', ')))
      }
    }
  }
  hyper0
}

#' Starting point for REML and symbolic Cholesky analysis
#'
#' Optionally refines the starting coefficients, converts hyperparameters to
#' the optimisation scale, and analyses the sparsity pattern of the
#' (scaled) Hessian once so later factorisations can reuse it.
#'
#' @return A list with `lambda0` (with `beta` and `first` attributes) and the
#'   updated `likdata` (with `chol_factor`).
#' @noRd
.evgmrf_start <- function(p0v, hyper0, Qd, likdata, likfns, refine) {
  attr(hyper0, 'beta') <- p0v
  attr(hyper0, 'first') <- TRUE
  if (refine)
    p0v <- .newton_step_inner(p0v, .d0_Q, .search_Q, diag = TRUE,
                              likdata = likdata, likfns = likfns,
                              Q = .mQ(hyper0, Qd),
                              control = likdata$control$inner)
  lambda0 <- .hyper2pars(hyper0, Qd$hyper_swap)
  hyper0 <- .pars2hyper(lambda0, hyper0, Qd$hyper_swap)
  attr(lambda0, 'beta') <- p0v
  attr(lambda0, 'first') <- TRUE

  Q0 <- .mQ(hyper0, Qd)
  H0 <- .d12_Q(p0v, likdata, likfns, Q0, hyper0)$H
  d <- pmax(Matrix::diag(H0), 1e-8)
  D <- Matrix::Diagonal(nrow(H0), 1 / sqrt(d))
  H0 <- D %*% H0 %*% D
  likdata$chol_factor <- if (likdata$control$super) {
    .chol_analyze_supernodal(H0)
  } else {
    .chol_analyze_simplicial(H0)
  }
  list(lambda0 = lambda0, likdata = likdata)
}

#' Jitter starting hyperparameters
#'
#' Tries moving each log/probit-scale hyperparameter by +/-1, then takes up
#' to five joint steps in the improving directions.
#' @noRd
.evgmrf_jitter <- function(lambda0, likdata, likfns, Qd) {
  f0 <- .reml0(lambda0, likdata = likdata, likfns = likfns, Qd = Qd, makeQ = .mQ)
  attr(lambda0, 'beta') <- attr(f0, 'beta')
  attr(lambda0, 'first') <- FALSE
  f1 <- 0 * as.vector(lambda0)
  for (i in seq_along(lambda0)) {
    lambda1 <- lambda0
    lambda1[i] <- lambda0[i] + 1
    f1[i] <- .reml0(lambda1, likdata = likdata, likfns = likfns, Qd = Qd, makeQ = .mQ)
  }
  adder <- c(0, 1)[1 + as.numeric(f1 < f0)]
  other_way <- which(f1 > f0)
  if (length(other_way) > 0) {
    for (i in other_way) {
      lambda1 <- lambda0
      lambda1[i] <- lambda0[i] - 1
      f1[i] <- .reml0(lambda1, likdata = likdata, likfns = likfns, Qd = Qd, makeQ = .mQ)
    }
    not_adder <- which(f1[other_way] < f0)
    if (length(not_adder) > 0)
      adder[other_way[not_adder]] <- -1
  }
  attr(lambda0, 'beta') <- attr(f0, 'beta')
  it <- 0
  while (it < 5) {
    lambda1 <- lambda0 + adder
    f1 <- .reml0(lambda1, likdata = likdata, likfns = likfns, Qd = Qd, makeQ = .mQ)
    if (f1 >= f0)
      break
    lambda0 <- lambda1
    f0 <- f1
    attr(lambda0, 'beta') <- attr(f0, 'beta')
    it <- it + 1
  }
  lambda0
}

#' Outer REML optimisation over the hyperparameters
#'
#' @return A list with `out` (the optimiser result) and `likdata` (whose
#'   `control$outer` may have been updated).
#' @noRd
.evgmrf_reml <- function(lambda0, likdata, likfns, Qd, outer, control, trace,
                         jitter) {
  # nothing to estimate: evaluate once at the fixed hyperparameters
  if (length(lambda0) == 0) {
    out0 <- .reml0(lambda0, likdata = likdata, likfns = likfns, Qd = Qd, makeQ = .mQ)
    atts <- attributes(out0)
    attributes(out0) <- NULL
    attributes(out0) <- atts[grep('penalized', names(atts))]
    out <- c(list(objective = out0), atts[-grep('penalized', names(atts))])
    return(list(out = out, likdata = likdata))
  }

  if (jitter)
    lambda0 <- .evgmrf_jitter(lambda0, likdata, likfns, Qd)

  if (outer == 'nelder-mead') {
    if (length(lambda0) > 1) {
      out <- .nelder_mead_discrete_list(lambda0, .reml0, likdata = likdata,
                                        likfns = likfns, Qd = Qd, makeQ = .mQ,
                                        trace = trace, step_size = control$step_size)
    } else {
      lambda0[] <- 0
      out <- .brent(lambda0, .reml0, likdata = likdata, likfns = likfns,
                    Qd = Qd, makeQ = .mQ, trace = trace)
    }
    return(list(out = out, likdata = likdata))
  }

  likdata$control$outer[c('stepmax', 'gradtol', 'steptol', 'itlim', 'dgradtol',
                          'fntol', 'alpha0')] <-
    control[c('reml_stepmax', 'reml_gradtol', 'reml_steptol', 'reml_itlim',
              'reml_dgradtol', 'reml_fntol', 'line_search_mult')]
  if (outer == 'newton') {
    likdata$control$outer$rho0 <- .25
    out <- .newton_step_inner(lambda0, .reml0, .reml_step, likdata = likdata,
                              likfns = likfns, Qd = Qd, makeQ = .mQ,
                              eps = control$reml_eps,
                              direction = control$reml_direction,
                              trace = trace > 0, control = likdata$control$outer,
                              attr2pass = c('beta', 'betal'))
  } else {
    out <- .BFGS(lambda0, .reml0, .reml1, likdata = likdata,
                 likfns = likfns, Qd = Qd, makeQ = .mQ,
                 eps = control$reml_eps,
                 direction = control$reml_direction,
                 trace = trace > 0, control = likdata$control$outer,
                 attr2pass = c('beta', 'betal'))
  }
  list(out = out, likdata = likdata)
}

#' Assemble the fitted `evgmrf` object
#' @noRd
.evgmrf_output <- function(out, likdata, likfns, family, dat, des, model, order,
                           formula, Qd, no_data, call) {
  nx <- dat$nx
  ny <- dat$ny
  index <- dat$index
  gmrf <- des$gmrf
  if (is.null(nx)) {
    nx <- likdata$n
    ny <- 1
  }
  if (is.null(index)) {
    index <- as.matrix(expand.grid(row = 1:nx, col = 1:ny))
  } else {
    nx <- max(index[, 1])
    ny <- max(index[, 2])
  }

  out$likdata <- likdata
  out$likfns <- likfns
  out$family <- family
  out$Hessian <- attr(out$objective, 'Hessian')
  out$X <- likdata$Xl
  out$holes <- dat$holes
  out$nx <- nx
  out$xid <- 1:nx
  out$ny <- ny
  out$yid <- 1:ny
  out$n <- likdata$n
  out$np <- likdata$np0
  out$unlink <- likfns$trans
  out$names <- likfns$names
  out$call <- call
  out$quantile <- likfns$quantile
  out$quantile0 <- likfns$quantile0
  out$gmrf <- gmrf

  order[!gmrf] <- model[!gmrf] <- NA
  names(formula) <- names(order) <- names(model) <- likfns$names$response
  out$formula <- formula
  out$order <- order
  out$model <- model

  attr(out$beta, 'split') <- likdata$psplit
  out$index <- index
  out$init <- des$inits
  out$logLik <- list(unpenalized = -attr(out$objective, 'unpenalized'),
                     penalized = -attr(out$objective, 'penalized'),
                     restricted = -as.numeric(out$objective))
  out$nobs <- sum(!is.na(dat$z))
  out$par_type <- des$par_type

  beta_split <- split(out$beta, attr(out$beta, 'split'))
  out$fixed <- mapply(function(x, y) if (y > 0) x[1:y] else numeric(0), beta_split, des$nX1)
  out$fixed <- mapply(function(x, y) {names(x) <- colnames(y); x}, out$fixed, des$X1)
  out$fixed_id <- mapply('+', c(0, sapply(out$X, ncol)[-length(out$X)]),
                         sapply(des$nX1, seq_len))
  out$Qd <- Qd
  out$no_data <- no_data
  out$supernodal <- likdata$control$super
  class(out) <- 'evgmrf'
  out
}
