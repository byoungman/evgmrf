#' Supported Models for evgmrf Models
#'
#' @description 
#' Documentation for the GMRF structures supported by the \code{\link{evgmrf}} 
#' fitting function.
#'
#' @name model.evgmrf
#' 
#' @details
#' The following models are available:
#' \itemize{
#'   \item \code{"icar"}: Intrinsic Conditional Autoregressive model.
#'   \item \code{"car"}: Conditional Autoregressive model.
#'   \item \code{"bym2"}: Besag-York-Mollie (BYM) 2 model.
#'   \item \code{"bym3"}: As above but with (not necessarily Gaussian) random 
#'   effect structure with parameters that need estimating. See Details.
#'   \item \code{"ald"}: As above but with all parameters fixed.
#' }
#' 
#' For \code{'bym3'} and \code{'bym4'}, \code{'bymfns'} must be supplied to 
#' \code{\link{evgmrf}}. See examples below.
#' 
#' Argument \code{order} can be supplied to \code{\link{evgmrf}} to specify the
#' order of autoregressive forms. The default of `1` assumes a first-order NSEW
#' neighborhood structure. Order `2` extends this to diagonals, and second-order
#' NESW, and so on. Order `1:2` fits first and second order. If a \code{list} each
#' element gives the order for each distribution parameter. Relative variances
#' can be fixed with \code{alpha}. If \code{order = 1:2} and \code{alpha = .5}, 
#' and lambda is smoothing parameter, then first-order term has smoothing parameter
#' lambda, second-order has alpha lambda, third-order has alpha^2 lambda, etc.
#' 
#' @examples
#' 
#' \dontrun{
#' data(COtop5prcp)
#' COmxprcp <- COtop5prcp$prcp[, 1, , ]
#' m_gev_car <- evgmrf(COmxprcp, family = 'gev', model = c('car', 'icar', 'icar'))
#' # As above, but with constant shape parameter, i.e. no GMRF
#' m_gev_car2 <- evgmrf(COmxprcp, family = 'gev', model = c('car', 'icar', NA))
#' 
#' m_gev_bym <- evgmrf(COmxprcp, family = 'gev', model = c('bym2', 'icar', 'icar'))
#' 
#' fn3 <- function(x, par1)  - dnorm(x, 0, exp(par1), log = TRUE)
#' attr(fn3, 'inits') <- c(-2)
#' m_gev_bym3 <- evgmrf(COmxprcp, family = 'gev', model = c('bym3', 'icar', 'icar'), 
#'                      bymfns = list(fn3, NULL, NULL))
#'                      
#' fn4 <- function(x)  - dnorm(x, log = TRUE)
#' m_gev_bym4 <- evgmrf(COmxprcp, family = 'gev', model = c('bym4', 'icar', 'icar'), 
#'                      bymfns = list(fn4, NULL, NULL))
#' }
#' 
#' @seealso \code{\link{evgmrf}}
#' @keywords datasets
NULL
