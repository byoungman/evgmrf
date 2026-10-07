#' Fitting extreme value distributions with parameters that vary according to Gaussian Markov random fields
#'
#' Fits extreme value distributions where one or more parameters vary across space or
#' structures using a Gaussian Markov random field (GMRF) framework.
#'
#' @param z An array or matrix containing the response data. If `W` is supplied,
#'   a matrix is expected. If an \eqn{r}-largest family is chosen, an array must be provided.
#' @param family A character string specifying the extreme value error distribution.
#'   See \code{\link{family.evgmrf}} for details on supported families and required
#'   data formats. Defaults to \code{"gev"}.
#' @param model A character string specifying the GMRF structural framework. See
#'   \code{\link{model.evgmrf}}. Defaults to \code{"icar"}.
#' @param formula A \code{formula} or \code{list} of \code{formula}s for any fixed effects.
#'   Defaults to ` ~ -1`, which means intercept- and fixed-effect-free forms.
#' @param covariates A \code{list} of covariates with names matching those in \code{formula}
#'   and dimensions compatible with \code{z}.
#' @param weights An array or matrix or scalar of weights for each value in `z`.
#'   Defaults to `1`.
#' @param inits.method A character string specifying how initial parameter values should
#'   be chosen: \code{"different"} (Default) uses different values for each points
#'   based on point-wise optima; \code{"same"} uses the same values.
#' @param W An optional square adjacency matrix defining the spatial neighbor relationships.
#' @param order A scalar or vector specifying autoregressive order. See
#'   \code{\link{model.evgmrf}}. Defaults to `1`.
#' @param lambda0 A scalar or vector of initial smoothing parameter values. Defaults
#'   to `1`. Larger values give smoother starting points.
#' @param args A \code{list} of additional arguments required by some families. See
#'   \code{\link{family.evgmrf}}.
#' @param control A \code{list} of control parameters. See \code{\link{evgmrf.control}}.
#' @param trace An integer specifying the verbosity level of the optimization output.
#'   Larger values provide more comprehensive iterations. Defaults to `0` (silent).
#' @param outer A character string specifying the smoothing parameter optimizer.
#'   One of \code{"newton"}, \code{"nelder-mead"}, or \code{"bfgs"} (Default).
#' @param gamma A scalar multiplier applied directly to the log-likelihood function
#'   representing localized constant weights; defaults to `1`.
#' @param nx,ny Positive integers giving the numbers of grid points in the
#'   first (\code{nx}) and second (\code{ny}) dimensions of the spatial grid. See Details.
#' @param index A two-column \code{matrix} of integers identifying row and column
#'   positions; see Details.
#' @param auto.weights Logical; if \code{TRUE}, weights are automatically calculated
#'   using an effective sample size estimate. EXPERIMENTAL. Defaults to \code{FALSE}.
#' @param bymfns A \code{list} of functions if model is \code{"bym3"} or
#'   \code{"bym4"}. See \code{\link{model.evgmrf}}.
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
#' @param supernodal Logical; activates CHOLMOD supernodal sparse matrix factorization
#'   settings. Defaults to `FALSE`.
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
#' # select top order statistic per year, and 20 x 15 subgrid
#' COmxprcp <- COorder$prcp[1:20, 1, 1:16, 1:14]
#' m_gev <- evgmrf(COmxprcp, family = "gev")
#'
#' \dontrun{
#' # Fix the precision parameter kappa at 5 for every parameter
#' m_fixed <- evgmrf(COmxprcp, family = "gev", hyper = list(kappa = 5))
#'
#' # Estimate kappa and rho for a BYM2 model, starting rho at 0.5
#' m_bym2 <- evgmrf(COmxprcp, family = "gev", model = "bym2",
#'                  hyper_start = list(rho = 0.5))
#' }
#'
#' @export
evgmrf <- function(z,
                   family = 'gev',
                   model = 'icar',
                   formula,
                   covariates,
                   weights = 1,
                   inits.method = 'different',
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
                   hyper.jitter = missing(hyper_start),
                   supernodal = FALSE
) {
  # hyper.jitter's default depends on missing(hyper_start), so evaluate it
  # before missing arguments are replaced by NULL below
  force(hyper.jitter)
  call <- sys.call()
  
  ## standardise arguments ----------------------------------------------------
  if (missing(formula)) formula <- NULL
  if (missing(covariates)) covariates <- NULL
  if (missing(W)) W <- NULL
  if (missing(nx)) nx <- NULL
  if (missing(ny)) ny <- NULL
  if (missing(index)) index <- NULL
  if (missing(hyper_start)) hyper_start <- NULL
  model <- tolower(model)
  args <- replace(.args0, names(args), args)
  control <- replace(evgmrf.control(), names(control), control)
  .checks(model, order)
  
  ## data ---------------------------------------------------------------------
  dat <- .evgmrf_data(z, family, weights, nx, ny, index, W)
  zl <- dat$zl
  wl <- dat$wl
  u <- NULL
  if (family == 'poisgpd') {
    pg <- .evgmrf_poisgpd(zl, wl, args)
    zl <- pg$zl
    wl <- pg$wl
    args <- pg$args
    u <- pg$u
  }
  no_data <- sapply(lapply(zl, is.finite), sum) == 0
  zl[no_data] <- NA
  wl[no_data] <- NA
  
  ## likelihood setup and initial values --------------------------------------
  fam <- .evgmrf_family(family, zl, wl, dat$n, dat$here, args, gamma, bymfns,
                        inits.method, u)
  .ld <- fam$likdata
  .lf <- fam$likfns
  
  ## per-parameter model specification ----------------------------------------
  spec <- .evgmrf_expand(hyper, order, model, formula, .ld$np0)
  hyper <- spec$hyper
  order <- spec$order
  model <- spec$model
  formula <- spec$formula
  .ld$np <- sum(!is.na(model))
  
  ## design matrices and starting coefficients --------------------------------
  covariates <- .evgmrf_covariates(covariates, dat$n, dat$nx, dat$ny)
  des <- .evgmrf_design(formula, covariates, model, fam$inits, .ld$n, .ld$np0)
  .ld$psplit <- des$psplit
  
  ## GMRF precision structure -------------------------------------------------
  if (!is.null(W)) {
    nW <- if (isS4(W)) W@Dim[1] else nrow(W)
    if (nW != dat$n)
      stop(paste('Supplied W dimension not compatible with data size:', nW, '!=', dat$n))
  }
  Qd0 <- .makeQ_data(dat$nx, dat$ny, model, order, des$nX1, W, bymfns, hyper)
  
  .ld$Xl <- des$Xl
  .ld$Xlc <- des$Xlc
  .ld$X <- des$X
  .ld$X0 <- Matrix::Diagonal(.ld$n)
  .ld$id_bym2 <- des$id_bym2
  .ld$openmp <- control$openmp
  .ld$threads <- control$threads
  .ld$control <- c(.evgam.control(), control)
  .ld$control$super <- supernodal
  
  ## hyperparameter starting values -------------------------------------------
  hyper0 <- .evgmrf_hyper0(model, order, des$inits, bymfns, hyper, hyper_start,
                           .ld$np, .ld$np0)
  
  ## starting point and sparse Cholesky analysis ------------------------------
  st <- .evgmrf_start(des$p0v, hyper0, Qd0, .ld, .lf, control$refine.inits)
  .ld <- st$likdata
  
  ## REML optimisation --------------------------------------------------------
  fit <- .evgmrf_reml(st$lambda0, .ld, .lf, Qd0, outer, control, trace,
                      hyper.jitter)
  
  ## assemble output ----------------------------------------------------------
  .evgmrf_output(fit$out, fit$likdata, .lf, family, dat, des, model, order,
                 formula, Qd0, no_data, call)
}
