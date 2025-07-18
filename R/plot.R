#' Plots from a fitted \code{evgmrf} object.
#'
#' Produces plots of parameters or quantiles from a \code{evgmrf} fit.
#'
#' @param x a fitted \code{evgmrf} object
#' @param type a character string giving the type of prediction sought; see Details. Defaults to \code{"link"}
#' @param prob a scalar or vector of probabilities for quantiles to be estimated if \code{type == "quantile"}; defaults to 0.5
#' @param lims to do
#' @param nlev to do
#' @param edge.drop to do
#' @param index to do
#' @param xid to do
#' @param yid to do
#' @param nrc to do
#' @param ... arguments passed to \code{lattice::levelplot}
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
#' @seealso \link{evgmrf} \link{predict.evgmrf}
#'
#' @return A plot
#' 
#' @export
#' 
plot.evgmrf <- function(x, type = 'link', prob = NULL, 
                        lims = NULL, nlev = NULL, edge.drop, index = x$index, 
                        xid = 1:x$nx, yid = 1:x$ny, nrc = NULL, ...) {
  out <- predict.evgmrf(x, type = type, index = index, prob = prob, xid = xid, yid = yid)
  nms <- names(out)
  if (!missing(edge.drop))
    out <- lapply(out, .drop.edge, edrop = edge.drop)
  if (is.null(nrc))
    nrc <- rev(n2mfrow(length(out)))
  if (x$holes & is.null(index))
    stop("Can't plot an object with holes if index not supplied here or to evgmrf().")
  if (is.null(lims)) {
    lims <- lapply(out, range, na.rm = TRUE)
  }
  if (is.null(nlev)) {
    nlev <- rep(15, length(out))
  }
  if (length(nlev) == 1)
    nlev <- rep(nlev, length(out))
  plots <- lapply(1:length(out), function(i) 
    lattice::levelplot(out[[i]], at = pretty(lims[[i]], nlev), main = nms[i], ...))
  gridExtra::grid.arrange(
    grobs = plots,
    ncol = nrc[2], # Number of columns
    nrow = nrc[1]  # Number of rows
  )
}

.drop.edge <- function(mat, edrop) {
  if (length(edrop) == 1)
    edrop <- rep(edrop, 4)
  if (length(edrop) != 4)
    stop('Wrong length edge.drop supplied.')
  if (edrop[1] > 0) {
    mat[1:edrop[1], ] <- NA
  }
  if (edrop[2] > 0) {
    mat[, 1:edrop[2]] <- NA
  }
  if (edrop[3] > 0) {
    mat[(nrow(mat) - edrop[3] + 1):nrow(mat), ] <- NA
  }
  if (edrop[4] > 0) {
    mat[, (ncol(mat) - edrop[4] + 1):ncol(mat)] <- NA
  }
  mat
}

# .plot_evgmrf <- function(object, type = 'link') {
#   out <- predict.evgmrf(object, type)
#   par(mfrow = rev(n2mfrow(length(out))))
#   for (i in 1:length(out)) image(out[[i]])
# }
# 
# .plot_evgmrf2 <- function(object, type = 'link') {
#   out <- predict.evgmrf(object, type)
#   nrc <- rev(n2mfrow(length(out)))
#   plots <- lapply(out, lattice::levelplot)
#   gridExtra::grid.arrange(
#     grobs = plots,
#     ncol = nrc[2], # Number of columns
#     nrow = nrc[1]  # Number of rows
#   )
# }
# 
# .plot_evgmrf_lattice <- function(object, type = 'link') {
#   out <- .predict_evgmrf(object, type)
#   ij <- expand.grid(i = 1:nrow(out[[1]]), j = 1:ncol(out[[1]]))
#   dfs <- lapply(out, function(x) data.frame(i = ij[, 1], j = ij[, 2], value = as.vector(x)))
#   plots <- lapply(dfs, function(x) lattice::levelplot(value ~ i * j, x))
#   gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
# }
