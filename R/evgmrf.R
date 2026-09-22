#' Fitting extreme value distributions with parameters that vary according to Gaussian Markov random fields
#'
#' Fits extreme value distributions where one or more parameters vary across space or 
#' structures using a Gaussian Markov random field (GMRF) framework.
#'
#' @param z An array or matrix containing the response data. If `W` is supplied, 
#'   a matrix is expected. If an $r$-largest family is chosen, an array must be provided.
#' @param family A character string specifying the extreme value error distribution. 
#'   See \code{\link{family.evgmrf}} for details on supported families and required 
#'   data formats. Defaults to \code{"gev"}.
#' @param model A character string specifying the GMRF structural framework. See
#'   \code{\link{model.evgmrf}}. Defaults to \code{"icar"}.
#' @param W An optional square adjacency matrix defining the spatial neighbor relationships.
#' @param trace An integer specifying the verbosity level of the optimization output. 
#'   Larger values provide more comprehensive iterations. Defaults to `0` (silent).
#' @param gamma A scalar multiplier applied directly to the log-likelihood function 
#'   representing localized constant weights; defaults to `1`.
#' @param nx,ny Positive integers giving the numbers of grid points in the
#'   first (\code{nx}) and second (\code{ny}) dimensions of the spatial grid. See details.
#' @param index A two-column \code{matrix} of integers identifying row and column
#'   positions; see Details.
#' @param formula A \code{formula} or \code{list} of \code{formula}s for any fixed effects. 
#'   Defaults to ` ~ -1`, which means intercept- and fixed-effect-free forms.
#' @param covariates A \code{list} of covariates with names matching those in \code{formula}
#'   and dimensions compatible with \code{z}.
#' @param weights An array or matrix or scalar of weights for each value in `z`. 
#'   Defaults to `1`.
#' @param inits A character string specifying how initial parameter values should
#'   be chosen: \code{"different"} (Default) uses different values for each points
#'   based on point-wise optima; \code{"same"} uses the same values.
#' @param lambda0 A scalar or vector of initial smoothing parameter values. Defaults
#'   to `1`. Larger values give smoother starting points.
#' @param order A scalar or vector specifying autoregressive order. See 
#'   \code{\link{model.evgmrf}}. Defaults to `1`.
#' @param nx,ny Integers giving numbers of rows and columns in grid; see Details.
#' @param bymfns A \code{list} of functions if model is \code{"bym3"} or 
#'   \code{"bym4"}. See \code{\link{model.evgmrf}}.
#' @param args A \code{list} or arguments required to calculations. See 
#'   \code{\link{family.evgmrf}}.
#' @param control A \code{list} of control parameters. See \code{\link{evgmrf.control}}.
#' @param outer A character string specifying the smoothing parameter optimizer.
#'   One of \code{"newton"}, \code{"nelder-mead"}, or \code{"bfgs"} (Default).
#' @param auto.weights Logical; if \code{TRUE}, weights are automatically calculated
#'   using an effective sample size estimate. EXPERIMENTAL. Defaults to \code{FALSE}.
#' @param gamma EDF multiplier. To do.
#' @param hyper,hyper_start Fixed values and starting values, respectively, for
#'   the hyperparameters of the GMRF. Each is either a named \code{list} or a
#'   \code{list} of such lists, one for each parameter of the extreme value
#'   distribution (in the order given in \code{\link{family.evgmrf}}), with
#'   \code{NULL} for any parameter with nothing to supply. A hyperparameter
#'   that is omitted from \code{hyper}, or set to \code{-1} (the default for
#'   \code{kappa}), is estimated; any other value fixes it. \code{hyper_start}
#'   sets the starting values of those that are estimated, on their natural
#'   scales. See Details.
#' @param hyper.jitter Logical; if \code{TRUE}, start hyperparameter values
#'   are jittered slightly to identify improvement. Defaults to \code{TRUE} unless
#'   \code{hyper_start} is used.
#' 
#' @details
#' In general the last two dimensions of \code{z} will define the grid size in
#' terms of rows and columns; otherwise it is `compacted`. If \code{z} is
#' compacted \code{index} or \code{nx} and \code{ny} can be used to supply
#' grid details, and \code{z} will be a \code{matrix} or a 3-dimensional
#' \code{array} if \code{family == "rlarge"}.
#'
#' For a compacted \code{z}, each location (a list element, or a column or
#' slice of \code{z}) is placed on an \code{nx} by \code{ny} grid as follows.
#' \itemize{
#'   \item If \code{index} is supplied, its \eqn{k}th row gives the row and
#'     column position of the \eqn{k}th location. Grid cells not listed in
#'     \code{index} are treated as missing ("holes"). Here \code{nx} and
#'     \code{ny} default to \code{max(index[, 1])} and \code{max(index[, 2])};
#'     an error is raised if a supplied value is smaller than these.
#'   \item If \code{index} is not supplied, the grid is assumed to be full and
#'     locations are ordered with the first dimension varying fastest, i.e.
#'     as in \code{expand.grid(seq_len(nx), seq_len(ny))}. Then
#'     \code{nx * ny} must equal the number of locations. If only one of
#'     \code{nx} and \code{ny} is given, the other is set to the number of
#'     locations divided by it, and an error is raised if this is not a whole
#'     number.
#'   \item If \code{W} is supplied and neither \code{nx} nor \code{ny} is,
#'     the locations are treated as a single line, with \code{nx} equal to the
#'     number of locations and \code{ny = 1}.
#' }
#' 
#' Hyperparameters can be fixed, or given starting values, for some or all of
#' the parameters of the extreme value distribution using \code{hyper} and
#' \code{hyper_start}. The names available depend on \code{model} (see
#' \code{\link{model.evgmrf}}):
#' \itemize{
#'   \item \code{kappa}: positive precision multiplier, used by all models and
#'     estimated on the log scale.
#'   \item \code{rho}: value in \eqn{(0, 1)} used by \code{"car"} and
#'     \code{"bym2"}, estimated on the probit scale.
#'   \item \code{epsilon}: positive precision of the unstructured component of
#'     \code{"bym"}, estimated on the log scale.
#'   \item \code{nu}: vector of values in \eqn{(0, 1)} weighting the
#'     second- and third-order components when \code{order} includes 2 and/or
#'     3 (\code{"icar"} only), with one element per order above 1 and
#'     estimated on the probit scale.
#' }
#' For \code{"bym3"} the hyperparameters are the arguments of the relevant
#' function in \code{bymfns}, excluding the first. Each must be given a fixed
#' value in \code{hyper} or a starting value in \code{hyper_start}.
#'
#' A list of length one for \code{hyper}, such as the default
#' \code{list(kappa = -1)}, is used for every parameter; otherwise
#' \code{hyper} must be a list of lists. A single named list of any length for
#' \code{hyper_start} is used for every parameter. If the first parameter has
#' no starting values to supply, use \code{list()} rather than \code{NULL} for
#' it in \code{hyper_start}, otherwise the whole argument is treated as a
#' single list. If \code{hyper_start} is missing, starting values are chosen
#' automatically.
#' 
#' @return An object of class \code{evgmrf} containing localized parameter arrays 
#'   and structural optimization summaries.
#' 
#' @references 
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Models. Journal of Statistical Software. \doi{10.18637/jss.v103.i03}
#'
#' @seealso \code{\link{predict.evgmrf}}, \code{\link{COorder}}
#'
#' @examples
#' 
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = "gev")
#' 
#' # Fix the precision parameter kappa at 5 for every parameter
#' m_fixed <- evgmrf(COmxprcp, family = "gev", hyper = list(kappa = 5))
#'
#' # Estimate kappa and rho for a BYM2 model, starting rho at 0.5
#' m_bym2 <- evgmrf(COmxprcp, family = "gev", model = "bym2",
#'                  hyper_start = list(rho = 0.5))
#' 
#' @export
evgmrf <- function(z, 
                   family = 'gev', 
                   model = 'icar', 
                   formula, 
                   covariates, 
                   weights = 1, 
                   inits = 'different',
                   W, 
                   order = 1, 
                   lambda0,
                   args = list(),
                   control = list(),
                   trace = 0,
                   outer = 'bfgs',
                   gamma = 1,
                   nx, ny, index,
                   auto.weights = FALSE,
                   bymfns = NULL,
                   hyper = list(kappa = -1),
                   hyper_start,
                   hyper.jitter = missing(hyper_start)
) {
  model <- tolower(model)
  args <- replace(.args0, names(args), args)
  # if (family == 'poisproc') {
  #   z <- array(z, c(1, dim(z)))
  # }
  holes <- TRUE
  control <- replace(evgmrf.control(), names(control), control)
  .checks(model, order)
  # some basics
  # if (is.list(z)) {
  #   if (!is.null(W)) {
  #     if (nrow(W) == length(z))
  #       infill <- FALSE
  #   }
  # }
  rlargeish <- substr(family, 1, 6) == 'rlarge'
  if (!is.list(z)) {
    w <- 0 * z + weights
    dz <- dim(z)
    m <- dz[1]
    array_style <- length(dz) == 3 + as.integer(rlargeish)
    if (!array_style) {
      n <- tail(dz, 1)
      if (rlargeish) {
        if (length(dz) == 2)
          stop("z can't be a matrix for r-largest family")
        r <- dz[2]
        z <- lapply(1:n, function(i) as.matrix(z[ , , i]))
        w <- lapply(1:n, function(i) as.matrix(w[ , , i]))
      } else {
        z <- lapply(1:n, function(i) z[ , i])
        w <- lapply(1:n, function(i) w[ , i])
      }
    } else {
      holes <- FALSE
      nx <- tail(dz, 2)[1]
      ny <- tail(dz, 2)[2]
      n <- nx * ny
      index <- as.matrix(expand.grid(seq_len(nx), seq_len(ny)))
      if (!rlargeish) {
        zm <- matrix(z, m)
        here <- colSums(!is.na(zm)) > 0
        wm <- matrix(w, m)
        zl <- lapply(1:ncol(zm), function(i) zm[, i])
        wl <- lapply(1:ncol(zm), function(i) wm[, i])
      } else {
        za <- array(z, c(dz[1:2], nx * ny))
        here <- colSums(!is.na(za[,1,])) > 0
        wa <- array(w, c(dz[1:2], nx * ny))
        zl <- lapply(1:dim(za)[3], function(i) as.matrix(za[ , , i]))
        wl <- lapply(1:dim(za)[3], function(i) as.matrix(wa[ , , i]))
      }
    }
    #   
    # 
    # #       
    # #       nx <- max(index[, 1])
    # #       
    # # 
    # # 
    # # not_na <- apply(is.finite(z), tail(seq_along(dim(z)), 2), any)
    # # no_data <- as.vector(apply(!is.finite(z), tail(seq_along(dim(z)), 2), all))
    # # holes <- ifelse(all(not_na), FALSE, TRUE)
    # # if (holes) {
    # #   if (length(dim(z)) != 3 + as.integer(rlargeish)) {
    # #     if (!infill) {
    # #       if (trace >= 0) {
    # #         message('Argument infill changed to TRUE.')
    # #         infill <- TRUE
    # #       }
    # #     }
    # #   }
    # # }
    #   if (missing(index)) {
    #     index <- as.matrix(expand.grid(lapply(tail(dz, 2), seq_len)))
    #   }
    #   nx <- max(index[, 1])
    #   ny <- max(index[, 2])
    #   here <- matrix(FALSE, nx, ny)
    #   here[index] <- TRUE
    #   # if (length(dim(z)) != 3 + as.integer(rlargeish)) {
    #   #   zz <- ww <- matrix(NA, nrow = nrow(z), ncol = length(here))
    #   #   zz[, as.logical(here)] <- z
    #   #   if (any(!here))
    #   #     holes <- TRUE
    #   #   ww[, as.logical(here)] <- w
    #   #   dz <- c(dim(zz)[-length(dim(zz))], nx, ny)
    #   #   z <- array(zz, dz)
    #   #   w <- array(ww, dz)
    #   #   not_na <- apply(is.finite(z), tail(seq_along(dz), 2), any)
    #   #   no_data <- !not_na
    #   #   index <- as.matrix(expand.grid(1:nx, 1:ny))
    #   # }
    #   # here <- here & not_na
    #   # attr(index, 'paddable') <- TRUE
    # # }  
    # # check whether in grid form
    # if (length(dz) == 2) {
    #   # matrix-form
    #   if (rlargeish)
    #     stop("z can't be a matrix for r-largest family")
    #   nx <- nx
    #   ny <- ny
    #   zm <- matrix(z, m)
    #   wm <- matrix(w, m)
    #   zl <- lapply(1:ncol(zm), function(i) zm[, i])
    #   wl <- lapply(1:ncol(zm), function(i) wm[, i])
    #   holes <- TRUE
    # } else {
    #   if (rlargeish) {
    #     if (length(dz) == 4) {
    #       r <- dz[2]
    #       nx <- dz[3]
    #       ny <- dz[4]
    #       za <- array(z, c(dz[1:2], nx * ny))
    #       wa <- array(w, c(dz[1:2], nx * ny))
    #       zl <- lapply(1:dim(za)[3], function(i) as.matrix(za[ , , i]))
    #       wl <- lapply(1:dim(za)[3], function(i) as.matrix(wa[ , , i]))
    #     } else {
    #       zl <- lapply(1:dim(z)[3], function(i) as.matrix(z[ , , i]))
    #       wl <- lapply(1:dim(z)[3], function(i) as.matrix(w[ , , i]))
    #       holes <- TRUE
    #     }
    #   } else {
    #     nx <- dz[2]
    #     ny <- dz[3]
    #     zm <- matrix(z, m)
    #     wm <- matrix(w, m)
    #     zl <- lapply(1:ncol(zm), function(i) zm[, i])
    #     wl <- lapply(1:ncol(zm), function(i) wm[, i])
    #   }
    # }
    # # if (family == 'poisgpd') {
    # #   wl <- lapply(seq_along(zl), function(i) wl[[i]][order(zl[[i]], decreasing = TRUE, na.last = NA)])
    # #   zl <- lapply(seq_along(zl), function(i) zl[[i]][order(zl[[i]], decreasing = TRUE, na.last = NA)])
    # # }
  }
  if (holes) {
    zl <- z
    if (length(weights) == 1)
      wl <- lapply(zl, function(x) 0 * x + weights)
    n <- length(zl)
    if (!missing(index)) {
      if (missing(nx)) {
        nx <- max(index[, 1])
      } else {
        if (max(index[, 1]) > nx)
          stop('index and nx not compatible')
      }
      if (missing(ny)) {
        ny <- max(index[, 2])
      } else {
        if (max(index[, 2]) > ny)
          stop('index and ny not compatible')
      }
    } else {
      if (missing(nx) && missing(ny) && missing(W) && missing(index))
        stop('Either nx, ny, index or W need supplying.')
      if (missing(nx)) {
        if (!missing(ny)) {
          nx <- n / ny
          if (nx - floor(nx) != 0)
            stop('ny incompatible with data length')
        }
      }
      if (missing(ny)) {
        if (!missing(nx)) {
          ny <- n / nx
          if (ny - floor(ny) != 0)
            stop('nx incompatible with data length')
        }
      }
      if (!missing(W) && missing(nx))
        nx <- n
      if (!missing(W) && missing(ny))
        ny <- 1
      index <- as.matrix(expand.grid(seq_len(nx), seq_len(ny)))
    }
    here <- matrix(FALSE, nx, ny)
    here[index] <- TRUE
    here <- as.logical(here)
    if (is.list(z)) {
      if (rlargeish) {
        z <- .list2array(zl)
        m <- dim(z)[1]
        r <- dim(z)[2]
        w <- .list2array(wl)
      } else {
        z <- .list2mat(zl)
        m <- nrow(z)
        w <- .list2mat(wl)
      }
    }
    if (!rlargeish) {
      zz <- ww <- matrix(NA, m, n)
      zz[, here] <- z
      ww[, here] <- w
    } else {
      zz <- ww <- array(NA, c(m, r, n))
      zz[, , here] <- z
      ww[, , here] <- w
    }
    dz <- c(dim(zz)[-length(dim(zz))], nx, ny)
    z <- array(zz, dz)
    here <- rowSums(matrix(apply(is.finite(z), -1, any), n)) > 0
    holes <- any(!here)
    w <- array(ww, dz)
  }
  if (family == 'poisgpd') {
    wl <- lapply(seq_along(zl), function(i) wl[[i]][order(zl[[i]], decreasing = TRUE, na.last = NA)])
    zl <- lapply(seq_along(zl), function(i) zl[[i]][order(zl[[i]], decreasing = TRUE, na.last = NA)])
    if (is.null(args$u)) {
      if (is.null(args$r)) {
        stop("Supply either args$r or args$u for family = 'poisgpd'.")
      } else {
        args$u <- sapply(zl, function(x) x[args$r])
      }
    }
    if (length(args$u) == 1)
      args$u <- rep(args$u, length(zl))
    u <- as.vector(args$u)
    if (length(u) != length(zl))
      stop("args$u is incompatible with z.")
    excl <- lapply(seq_along(zl), function(i) zl[[i]] > args$u[i])
    zl <- lapply(seq_along(zl), function(i) zl[[i]][excl[[i]]])
    wl <- lapply(seq_along(wl), function(i) wl[[i]][excl[[i]]])
  }
  no_data <- sapply(lapply(zl, is.finite), sum) == 0
  zl[no_data] <- NA
  wl[no_data] <- NA
  .ld <- list(n = n, z = zl, w = wl, mult = 1, args = args, gamma = gamma)
  .ld$bymfns <- bymfns
  inits_list <- list()
  if (family == 'gpd') {
    .lf <- .gpd_fns
    .ld$np <- 2
    inits_list$same <- .quick_tgpd(na.omit(unlist(zl)))
    if (inits == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tgpd(zl[[i]]))
  }
  if (family == 'poisgpd') {
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
    inits_list$same <- .quick_tpp(na.omit(unlist(zl)), .ld$m, min(.ld$u, na.rm = TRUE), args$delta)
    if (inits == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tpp(zl[[i]], .ld$u[i], m = .ld$m, delta = args$delta))
  }
  if (family == 'gev') {
    .lf <- .gev_fns
    .ld$np <- 3
    inits_list$same <- .quick_tgev(na.omit(unlist(zl)), args$delta)
    if (inits == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tgev_shrink(zl[[i]], delta = args$delta, pars0 = inits_list$same, mult = args$mult))
  }
  if (family == 'rlarge') {
    .lf <- .rlarge_fns
    inits_list$same <- .quick_tgev(na.omit(unlist(lapply(zl, function(x) x))), args$delta)
    if (inits == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_tgev_shrink(zl[[i]], delta = args$delta, pars0 = inits_list$same, mult = args$mult))
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
    inits_list$same <- .quick_ald(na.omit(unlist(zl)), args = args)
    if (inits == 'different')
      inits_list$diff <- sapply(which(here), function(i) .quick_ald(zl[[i]], args = args))
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
  if (inits == 'same') {
    inits <- matrix(inits_list$same, .ld$np, n)
  }
  if (inits == 'different') {
    if (any(!here)) {
      inits <- matrix(inits_list$same, .ld$np, n)
      inits[, here] <- inits_list$diff
    } else {
      inits <- inits_list$diff
    }
  }
  inits <- t(inits)
  .ld$np0 <- .ld$np
  if (length(hyper) == 1 && .ld$np0 > 1)
    hyper <- lapply(1:.ld$np0, function(.) hyper)
  if (length(hyper) == 1 && !is.list(hyper[[1]]))
    hyper <- lapply(1:.ld$np0, function(.) hyper)
  if (!is.list(order)) {
    if (is.vector(order)) {
      if (length(order) == 1)
        order <- rep(order, .ld$np)
      order <- lapply(order, seq_len)
    }
  } else {
    if (length(order) == 1)
      order <- rep(order, .ld$np)
  }
  order <- lapply(order, sort)
  if (.ld$np > 1) {
    if (length(model) == 1)
      model <- rep(model, .ld$np)
  }
  if (missing(formula)) {
    formula <- lapply(model, .model2formula)
  } else {
    if (is(formula, 'formula'))
      formula <- lapply(model, function(.) formula)
  }
  .ld$np <- sum(!is.na(model))
  if (missing(covariates)) {
    covariates <- data.frame(id = 1:n)
    if (!is.null(nx) & !is.null(ny)) {
      xy <- expand.grid(x = seq_len(nx), y = seq_len(ny))
      covariates$x <- xy$x
      covariates$y <- xy$y
    }
  } else {
    cov_dims <- lapply(covariates, dim)
    if (!is.null(nx) && !is.null(ny)) {
      compatible <- sapply(cov_dims, function(x) if (length(x) == 2) all(x == c(nx, ny)) else TRUE)
      if (any(!compatible))
        stop(paste0('Covariates ', paste(names(covariates), collapse = ' and '), ' not compatible with z'))
    } else {
      covariates <- lapply(covariates, as.vector)
      compatible <- sapply(covariates, length) == length(zl)
      if (any(!compatible))
        stop(paste0('Covariates ', paste(names(covariates), collapse = ' and '), ' not compatible with z'))
    }
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
  fixed_names <- lapply(X1, colnames)
  p0l <- id_bym2 <- par_type <- list()
  for (i in 1:.ld$np0) {
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
  # if (!is.null(bymfns))
  #   model[!sapply(bymfns, is.null) & is.na(model)] <- 'bym4'
  # id_bym2 <- lapply(p0l, function(x) rep(FALSE, length(x)))
  for (i in which(gmrf)) {
    # p0l[[i]] <- c(p0l[[i]], p0m[, i])
    id_bym2[[i]] <- rep(c(FALSE, FALSE), sapply(p0l[[i]], length))
    par_type[[i]] <- c(par_type[[i]], 'spatial')
    if (model[i] %in% paste('bym', 2:3, sep = '')) {
      id_bym2[[i]] <- c(id_bym2[[i]], rep(TRUE, .ld$n))
      p0l[[i]] <- c(p0l[[i]][[1]], .1 * p0l[[i]][[2]], .9 * p0l[[i]][[2]])
      par_type[[i]] <- c(par_type[[i]], 'random')
    }
  }
  p0l <- lapply(p0l, unlist)
  .ld$psplit <- rep(seq_along(p0l), sapply(p0l, length))
  ## put together GMRF stuff
  if (missing(W))
    W <- NULL
  if (!is.null(W)) {
    nW <- if (isS4(W)) W@Dim[1] else nrow(W)
    if (nW != n)
      stop(paste('Supplied W dimension not compatible with data size:', nW, '!=', n))
  }
  Qd0 <- .makeQ_data(nx, ny, model, order, nX1, W, bymfns, hyper)
  .ld$Xl <- X1
  for (i in which(gmrf)) {
    .ld$Xl[[i]] <- list(.ld$Xl[[i]], Matrix::Diagonal(.ld$n))
    if (model[i] %in% paste('bym', 2:3, sep = ''))
      .ld$Xl[[i]] <- c(.ld$Xl[[i]], Matrix::Diagonal(.ld$n))
  }
  .ld$Xlc <- .ld$Xl # componentwise
  for (i in seq_along(.ld$Xl)) {
    if (is.list(.ld$Xl[[i]])) 
      .ld$Xl[[i]] <- do.call(cbind, .ld$Xl[[i]])
  }
  # .ld$Xl <- lapply(seq_along(.ld$Xl), function(i) do.call(cbind, .ld$Xl[[i]]))
  # Qd0$R <- list()
  # for (i in which(gmrf))
  #   Qd0$R[[i]] <- Matrix::qrR(Matrix::qr(rbind(.ld$Xl[[i]], .ld$Xl[[i]])))
  .ld$X <- Matrix::.bdiag(.ld$Xl)
  .ld$X0 <- Matrix::Diagonal(.ld$n)
  .ld$id_bym2 <- id_bym2
  .ld$openmp <- control$openmp
  .ld$threads <- control$threads
  .ld$control <- c(.evgam.control(), control)
  # lambda0_temp <- as.vector(unlist(lapply(model, .inits_model)))
  # par_var <- -log(20 * apply(p0m, 2, var))
  par_var <- apply(inits, 2, sd)
  if(is.null(bymfns)) 
    bymfns <- lapply(seq_along(model), function(.) NULL)
  for (i in seq_along(hyper)) {
    if (is.null(hyper[[i]]))
      hyper[[i]] <- list()
  }
  hyper0 <- lapply(1:.ld$np, function(i) .inits_model(model[i], order[[i]], par_var[i], bymfns[[i]], hyper[[i]]))
  if (!missing(hyper_start)) {
    if (!is.list(hyper_start[[1]]))
      hyper_start <- lapply(1:.ld$np0, function(.) hyper_start)
    for (i in seq_along(hyper_start)) {
      if (!is.null(hyper_start[[i]]))
        hyper0[[i]][names(hyper_start[[i]])] <- hyper_start[[i]]
    }
  } else {
    if ('bym3' %in% model) {
      bym3s <- which(model %in% 'bym3')
      for (i in bym3s) {
        fmls <- formals(bymfns[[i]])[-1]
        if (length(fmls) > 0) {
          fmls_need <- names(fmls)
          fmls_need <- fmls_need[!(names(fmls) %in% names(hyper[[i]]))]
          if (!missing(hyper_start))
            fmls_need <- fmls_need[names(fmls) %in% names(hyper_start[[i]])]
          if (length(fmls_need) > 0)
            stop(paste0('Fixed or starting values missing for BYM3 model parameters ', paste0(fmls_need, collapse = ', ')))
        }
      }
    }
  }
  p0v <- as.vector(inits)
  p0v <- unlist(p0l)
  attr(hyper0, 'beta') <- p0v
  attr(hyper0, 'first') <- TRUE
  if (control$refine.inits)
     p0v <- .newton_step_inner(p0v, .d0_Q, .search_Q, diag = TRUE,
                         likdata = .ld, likfns = .lf, Q = .mQ(hyper0, Qd0),
                         control = .ld$control$inner)
  lambda0 <- .hyper2pars(hyper0, Qd0$hyper_swap)
  attr(lambda0, 'beta') <- p0v
  attr(lambda0, 'first') <- TRUE
  if (control$inner_optim == 'Cholesky')
    .ld$chol0 <- .Cholesky0(lambda0, Qd0, control$super)
  if (length(lambda0) == 0) {
    out0 <- .reml0(lambda0, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ)
    atts <- attributes(out0)
    attributes(out0) <- NULL
    attributes(out0) <- atts[grep('penalized', names(atts))]
    out <- list(objective = out0)
    out <- c(out, atts[-grep('penalized', names(atts))])
  } else {
    if (hyper.jitter) {
      lambda00 <- lambda0
      f0 <- .reml0(lambda0, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ)
      attr(lambda0, "beta") <- attr(f0, 'beta')
      attr(lambda0, "first") <- FALSE
      f1 <- 0 * as.vector(lambda0)
      for (i in seq_along(lambda0)) {
        lambda1 <- lambda0
        lambda1[i] <- lambda0[i] + 1
        f1[i] <- .reml0(lambda1, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ)
      }
      adder <- c(0, 1)[1 + as.numeric(f1 < f0)]
      other_way <- which(f1 > f0)
      if (length(other_way) > 0) {
        for (i in other_way) {
          lambda1 <- lambda0
          lambda1[i] <- lambda0[i] - 1
          f1[i] <- .reml0(lambda1, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ)
        }
        not_adder <- which(f1[other_way] < f0)
        if (length(not_adder) > 0)
          adder[other_way[not_adder]] <- -1
      }
      attr(lambda0, "beta") <- attr(f0, 'beta')
      cond <- TRUE
      it <- 0
      while(cond & it < 5) {
        lambda1 <- lambda0 + adder
        f1 <- .reml0(lambda1, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ)
        if (f1 < f0) {
          lambda0 <- lambda1
          f0 <- f1
          attr(lambda0, "beta") <- attr(f0, 'beta')
          it <- it + 1
        } else {
          cond <- FALSE
        }
      }
    }
    if (outer == 'nelder-mead') {
      if (length(lambda0) > 1) {
        out <- .nelder_mead_discrete_list(lambda0, .reml0, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ, 
                                          trace = trace, step_size = control$step_size)
      } else {
        lambda0[] <- 0
        out <- .brent(lambda0, .reml0, likdata = .ld, likfns = .lf, Qd = Qd0, makeQ = .mQ, trace = trace)
      }
    } else {
      .ld$control$outer[c('stepmax', 'gradtol', 'steptol', 'itlim', 'dgradtol', 
                          'fntol', 'alpha0')] <- control[c('reml_stepmax', 'reml_gradtol', 'reml_steptol', 
                                                           'reml_itlim', 'reml_dgradtol', 'reml_fntol', 'line_search_mult')]
      if (outer == 'newton') {
        .ld$control$outer$rho0 <- .25
        out <- .newton_step_inner(lambda0, .reml0, .reml_step, likdata = .ld, 
                                  likfns = .lf, Qd = Qd0, makeQ = .mQ, eps = control$reml_eps, 
                                  direction = control$reml_direction, trace = trace > 
                                    0, control = .ld$control$outer, attr2pass = c('beta', 'betal'))
      } else {
        out <- .BFGS(lambda0, .reml0, .reml1, likdata = .ld, 
                     likfns = .lf, Qd = Qd0, makeQ = .mQ, eps = control$reml_eps, 
                     direction = control$reml_direction, trace = trace > 
                       0, control = .ld$control$outer, attr2pass = c('beta', 'betal'))
        
      }
    }
  }
  out$likdata <- .ld
  out$likfns <- .lf
  out$family <- family
  if (is.null(nx)) {
    nx <- .ld$n
    ny <- 1
  }
  if (is.null(index)) {
    if (!is.null(nx) & !is.null(ny))
      index <- as.matrix(expand.grid(row = 1:nx, col = 1:ny))
  } else {
    nx <- max(index[, 1])
    ny <- max(index[, 2])
    # holes <- nrow(index) < nx * ny
  }
  if (control$inner_optim != 'Cholesky') {
    out$cholprecondHessian <- suppressWarnings(try(Matrix::Cholesky(out$precondHessian, super = control$super, LDL = FALSE), silent = TRUE))
    if (inherits(out$cholprecondHessian, 'try-error')) {
      out$precondHessian <- .perturb_super(out$precondHessian)
      out$cholprecondHessian <- Matrix::Cholesky(out$precondHessian, LDL = FALSE)
    }
  }
  out$Hessian <- attr(out$objective, 'Hessian')
  out$X <- .ld$Xl
  out$holes <- holes
  out$nx <- ifelse(is.null(nx), .ld$n, nx)
  out$xid <- 1:nx
  out$ny <- ifelse(is.null(ny), 1, ny)
  out$yid <- 1:ny
  out$n <- .ld$n
  out$np <- .ld$np0
  out$unlink <- .lf$trans
  out$names <- .lf$names
  out$call <- sys.call()
  out$quantile <- .lf$quantile
  out$quantile0 <- .lf$quantile0
  out$gmrf <- gmrf
  order[!gmrf] <- model[!gmrf] <- NA
  names(formula) <- names(order) <- names(model) <- .lf$names$response
  out$formula <- formula
  out$order <- order
  out$model <- model
  attr(out$beta, 'split') <- .ld$psplit
  out$index <- index
  out$init <- inits
  out$logLik <- list(unpenalized = -attr(out$objective, 'unpenalized'),
                     penalized = -attr(out$objective, 'penalized'),
                     restricted = -as.numeric(out$objective))
  out$nobs <- sum(!is.na(z))
  out$par_type <- par_type
  beta_split <- split(out$beta, attr(out$beta, 'split'))
  out$fixed <- mapply(function(x, y) if (y > 0) x[1:y] else numeric(0), beta_split, nX1)
  out$fixed <- mapply(function(x, y) {names(x) <- colnames(y); x}, out$fixed, X1)
  out$fixed_id <- mapply('+', c(0, sapply(out$X, ncol)[-length(out$X)]), sapply(nX1, seq_len))
  out$Qd <- Qd0
  out$no_data <- no_data
  class(out) <- 'evgmrf'
  out
}

