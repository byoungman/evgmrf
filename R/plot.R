#' Plot Method for Fitted \code{evgmrf} Objects
#'
#' Produces spatial or line plots of estimated parameters, linear predictors, or
#' predicted quantiles from a fitted extreme value Markov random field (\code{evgmrf}) model.
#'
#' @param x An object of class \code{evgmrf}, typically representing a fitted
#'   spatial extreme value model.
#' @param which A numeric vector of indices or character vector of names specifying
#'   which parameter estimates, linear predictors, or quantiles to plot. Defaults to
#'   plotting all components returned by \code{\link{predict.evgmrf}}.
#' @param polygons An optional polygon dataset (such as spatial polygons or boundary geometries)
#'   used to map predicted values spatially. If omitted, predictions are displayed as grid maps
#'   via \code{\link[lattice]{levelplot}} or 1D lines via \code{\link[lattice]{xyplot}}.
#' @param lims An optional list specifying the plot limits (e.g., color bar ranges) for
#'   each plot component. Defaults to the range of values in each plotted component.
#' @param nlev An integer or vector of integers giving the number of color levels/intervals
#'   used in image or polygon plots. Defaults to \code{15L}.
#' @param edge.drop An optional integer or logical specification passed internally to
#'   drop outer edge pixels/nodes from spatial grid outputs prior to plotting.
#' @param mfrow A two-element integer vector \code{c(rows, cols)} specifying the grid
#'   layout dimensions for arranging multiple panel plots by row.
#' @param mfcol A two-element integer vector \code{c(rows, cols)} specifying the grid
#'   layout dimensions for arranging multiple panel plots by column.
#' @param show Logical; if \code{TRUE} plots are shown in graphics device and if \code{FALSE} a 
#'   \code{list} is returned. Defaults to \code{TRUE}.
#' @param plot.args A \code{list} of additional arguments passed to the underlying
#'   plotting functions (e.g., \code{\link[lattice]{levelplot}} or \code{\link[lattice]{xyplot}}).
#'   Defaults to \code{list(asp = 1)}.
#' @param ... Arguments passed on to \code{\link{predict.evgmrf}} (e.g., \code{type},
#'   \code{prob}, or \code{decompose}).
#'
#' @details
#' This function computes model predictions by calling \code{\link{predict.evgmrf}} with
#' any extra arguments supplied via \code{...}, then renders visual displays of the output:
#' \itemize{
#'   \item If \code{polygons} is supplied, values are rendered onto spatial polygon regions.
#'   \item If \code{polygons} is missing and a prediction component is 2D, a heat map grid is plotted using \code{\link[lattice]{levelplot}}.
#'   \item If a prediction component is 1D (e.g., single row or column), a line plot is rendered using \code{\link[lattice]{xyplot}}.
#' }
#' Multiple plot panels are automatically arranged in a grid using \code{\link[gridExtra]{grid.arrange}}.
#'
#' @references
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Models. Journal of Statistical Software. \doi{10.18637/jss.v103.i03}
#'
#' @examples
#' data(COorder)
#' COmxprcp <- COorder$prcp[, 1, , ]
#' m_gev <- evgmrf(COmxprcp)
#' plot(m_gev)
#' plot(m_gev, type = 'response')
#' plot(m_gev, prob = .99)
#'
#' @seealso \code{\link{evgmrf}}, \code{\link{predict.evgmrf}}, \code{\link[lattice]{levelplot}}
#'
#' @return Invisibly returns the arranged grid of plot objects (a \code{grob} layout from \code{gridExtra}).
#'
#' @export
plot.evgmrf <- function(x, which, polygons, lims, nlev = 15L, edge.drop, 
                        mfrow, mfcol, plot.args = list(asp = 1), show = TRUE, ...) 
{
  out <- do.call(predict.evgmrf, list(object = x, ...))
  if (!all(unique(sapply(out, class)) %in% c('matrix', 'array')))
    out <- unlist(out, recursive = FALSE)
  if (missing(which)) {
    which <- seq_along(out)
  } else {
    if (!is.numeric(which))
      which <- unlist(lapply(which, grep, names(out)))
  }
  out <- out[which]
  nms <- names(out)
  has_se <- grep('se.', nms)
  if ((x$holes && is.null(x$index)) && (x$holes && is.null(x$index)))
    stop("Can't plot an object with holes if index not supplied here or to evgmrf().")
  if (missing(polygons)) {
    if (!missing(edge.drop))
      out <- lapply(out, .drop.edge, edrop = edge.drop)
  }
  if (missing(lims)) {
    lims <- lapply(out, range, na.rm = TRUE)
  }
  if (length(nlev) == 1)
    nlev <- rep(nlev, length(out))
  plots <- list()
  if(length(plot.args) == 1)
    plot.args <- lapply(seq_along(out), function(.) plot.args)
  plot.args <- lapply(plot.args, function(x) {
    if (is.null(x$asp)) x$asp <- 1
    x
  })
  if (length(has_se)) {
    for (i in has_se) {
      if (is.null(plot.args[[i]]$col.regions))
        plot.args[[i]]$col.regions <- hcl.colors(100)
    }
  }
  if (missing(polygons)) {
    for (i in 1:length(out)) {
      plot.args.i <- plot.args[[i]]
      plot.args.i$at <- pretty(lims[[i]], nlev[i])
      plot.args.i$main <- nms[i]
      if (1 %in% dim(out[[i]])) {
        if (dim(out[[i]])[1] == 1) {
          xlab == 'column'
        } else {
          xlab = 'row'
        }
        plot.args.i$x <- data.frame(estimate = as.vector(out[[i]]), 
                                    index = seq_len(length(out[[i]])))
        plot.args.i$formula <- estimate ~ index
        plot.args.i$xlab <- xlab
        plot.args.i$type <- 'l'
        plots[[i]] <- do.call(lattice::xyplot, plot.args.i)
      } else {
        plot.args.i$x <- out[[i]]
        plots[[i]] <- do.call(lattice::levelplot, plot.args.i)
      }
    }
  } else {
    for (i in 1:length(out)) {
      plot.args.i <- plot.args[[i]]
      plot.args.i$n <- nlev[i]
      plot.args.i$main <- nms[i]
      plot.args.i$polys <- polygons
      plot.args.i$values <- out[[i]]
      plots[[i]] <- do.call(.plot.polygons, plot.args.i)
    }
  }
  if (!('decompose' %in% ...names())) {
    decompose <- FALSE
  } else {
    decompose <- list(...)$decompose
  }
  if (missing(mfrow) && missing(mfcol)) {
    nrc <- rev(n2mfrow(length(out)))
    lm <- matrix(1:(nrc[1] * nrc[2]), nrc[1], byrow = length(has_se))
  } else {
    if (!missing(mfrow)) {
      lm <- matrix(1:(mfrow[1] * mfrow[2]), mfrow[1], byrow = TRUE)
    } else {
      lm <- matrix(1:(mfcol[1] * mfcol[2]), mfcol[1])
    }
  }
  if (show) {
    gridExtra::grid.arrange(grobs = plots, layout_matrix = lm)
  }
  invisible(plots)
}

.plot.polygons <- function(polys, values, n = 15, main = '', ...) {
  
  x.rg <- range(unlist(sapply(polys, '[[', 1)))
  y.rg <- range(unlist(sapply(polys, '[[', 2)))
  brks <- pretty(values, n = n)
  rmp <- lattice::trellis.par.get("regions")$col
  cols <- lattice::level.colors(x = values, at = brks, col.regions = rmp)
  lattice::xyplot(1 ~ 1, 
         xlim = x.rg, ylim = y.rg,
         xlab = "x", ylab = "y", main = main,

         legend = list(
           right = list(
             fun = lattice::draw.colorkey,
             args = list(key = list(at = brks, col = rmp, ticks = list(at = brks))))
         ),
         
         panel = function(...) {
           for (i in 1:length(polys)) {
             lattice::panel.polygon(x = polys[[i]][[1]], y = polys[[i]][[2]], 
                           col = cols[i], border = "darkgrey", lwd = 2)
           }
         }, ...
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
