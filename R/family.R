#' Supported Families for evgmrf Models
#'
#' @description
#' Documentation for the extreme value distributions (families) and data
#' structure expectations supported by the \code{\link{evgmrf}} fitting
#' function.
#'
#' @name family.evgmrf
#'
#' @details
#' Argument \code{family} of \code{\link{evgmrf}} is a character string giving
#' the distribution to be fitted, and defaults to \code{"gev"}. The following
#' families are available:
#' \itemize{
#'   \item \code{"gev"}: Generalized extreme value distribution for block
#'     maxima. Years are typically the blocks, so the data are annual maxima.
#'     Its parameters are the location, scale and shape.
#'   \item \code{"rlarge"}: \eqn{r}-largest order statistics model, an
#'     extension of the GEV in which the \eqn{r} largest values in each block
#'     are modelled, so that \eqn{r = 1} is the GEV. It has the same
#'     parameters as the GEV, and the same quantiles. The value of \eqn{r} is
#'     taken from the data; see below.
#'   \item \code{"gpd"}: Generalized Pareto distribution for excesses of a high
#'     threshold. The threshold is not supplied to \code{\link{evgmrf}}: the
#'     data must already be excesses, i.e. with the threshold subtracted. Its
#'     parameters are the scale and shape.
#'   \item \code{"poisgpd"}: Poisson-GPD point process model, an extension of
#'     the GPD that also uses the threshold, and is parametrized by the GEV's
#'     location, scale and shape. The data are the values themselves, not
#'     excesses. Requires \code{args$nper} and either \code{args$u} or
#'     \code{args$r}.
#'   \item \code{"ald"}: Asymmetric Laplace distribution, used for quantile
#'     regression, and hence for estimating thresholds for the GPD and
#'     Poisson-GPD. Its parameters are the location, which is the quantile
#'     being estimated, and the scale. Requires \code{args$tau}.
#' }
#'
#' Where the shape parameter of a GEV, \eqn{r}-largest or GPD model is
#' estimated to be within \eqn{10^{-6}} of zero, it is set to \eqn{10^{-6}} to
#' avoid numerical instability. For the GEV, \eqn{r}-largest and Poisson-GPD
#' families, the scale is modelled on the log scale and the shape on a
#' transformed scale, which \code{plot} and \code{predict} label
#' \code{"logscale"} and \code{"transshape"}. Setting \code{type = "response"}
#' gives \code{"scale"} and \code{"shape"}.
#'
#' Data \code{z} can be supplied as an array with dimensions indexing time (or
#' block), then grid row, then grid column, e.g. year, longitude and latitude.
#' For \code{"rlarge"} it is a four-dimensional array with dimensions indexing
#' block, order statistic (largest first), grid row and grid column, with
#' \eqn{r} the length of the second dimension. Alternatively, in a compacted
#' form, \code{z} can be a matrix with one column for each location (a
#' three-dimensional array of block, order statistic and location for
#' \code{"rlarge"}), or a list with one element for each location. If a matrix
#' or list is supplied, \code{nx} and \code{ny}, or \code{index}, describe the
#' grid, and \code{W} can be supplied instead. Lists allow numbers of values to
#' vary by location, as they do for excesses. Missing values are ignored, and
#' locations where all values are missing are part of the GMRF but do not
#' contribute to the likelihood.
#'
#' Family-specific arguments are supplied in \code{args}:
#' \itemize{
#'   \item \code{args$tau}: For \code{"ald"}, the quantile to be estimated. For
#'     example, \code{args = list(tau = .5)} gives a median estimate.
#'   \item \code{args$nper}: For \code{"poisgpd"}, the length of the period
#'     covered by the data, in units of time such as years. Estimates are then
#'     on that scale, e.g. annual, so are comparable with those from annual
#'     maxima.
#'   \item \code{args$u}: For \code{"poisgpd"}, the threshold above which
#'     values are used, either a single value or one for each location.
#'   \item \code{args$r}: For \code{"poisgpd"}, used instead of \code{args$u}
#'     to give a threshold for each location, which is its \eqn{r}th largest
#'     value. Data that have already been restricted to the largest
#'     \eqn{r} values at each location are then equivalent to
#'     \code{args$u} being their minimum at each location.
#' }
#'
#' An estimated quantile from the \code{"ald"} family can be used as a
#' threshold, by supplying \code{predict(m_ald)$location}, where \code{m_ald}
#' is the fitted model, as \code{args$u} to \code{"poisgpd"}.
#'
#' @references
#' Coles, S. G. (2001). \emph{An Introduction to Statistical Modeling of
#' Extreme Values}. Springer-Verlag. \doi{10.1007/978-1-4471-3675-0}
#'
#' Oh, H.-S., Lee, T. C. M. and Nychka, D. W. (2011). Fast nonparametric
#' quantile regression with arbitrary smoothing methods. \emph{Journal of
#' Computational and Graphical Statistics}, 20(2), 510-526.
#' \doi{10.1198/jcgs.2010.10063}
#'
#' Yu, K. and Moyeed, R. A. (2001). Bayesian quantile regression.
#' \emph{Statistics & Probability Letters}, 54(4), 437-447.
#' \doi{10.1016/s0167-7152(01)00124-9}
#'
#' @examples
#'
#' \donttest{
#' data(COorder)
#'
#' # GEV distribution for annual maxima
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp, family = "gev")
#'
#' # r-largest order statistics model, using the top five values in each year
#' m_rlarge <- evgmrf(COorder$prcp, family = "rlarge")
#'
#' # GPD for excesses of 20mm. The excesses are supplied as a list, without
#' # missing values, with nx and ny giving the grid
#' data(COexc)
#' excess <- COexc$prcp - 20
#' nx <- length(COexc$x)
#' ny <- length(COexc$y)
#' ind <- expand.grid(seq_len(nx), seq_len(ny))
#' excess_list <- lapply(1:nrow(ind),
#'                       function(i) na.omit(excess[, ind[i, 1], ind[i, 2]]))
#' m_gpd <- evgmrf(excess_list, family = "gpd", nx = nx, ny = ny)
#'
#' # Poisson-GPD with a threshold of 20mm, for 56 years of data
#' m_pp1 <- evgmrf(COexc$prcp, family = "poisgpd",
#'                 args = list(u = 20, nper = 56))
#'
#' # Poisson-GPD to the top 50 values at each grid cell, which are used to
#' # give cell-specific thresholds
#' data(COtop50)
#' m_pp2 <- evgmrf(COtop50$prcp, family = "poisgpd",
#'                 args = list(r = 50, nper = 56))
#'
#' # Estimate the median of the top 50 values, roughly the 99.88th percentile
#' # of daily precipitation, and use it as the threshold for a Poisson-GPD
#' m_ald <- evgmrf(COtop50$prcp, family = "ald", args = list(tau = .5))
#' m_pp3 <- evgmrf(COtop50$prcp, family = "poisgpd",
#'                 args = list(nper = 56, u = predict(m_ald)$location))
#' }
#'
#' @seealso \code{\link{evgmrf}}, \code{\link{model.evgmrf}}
#' @keywords models
NULL