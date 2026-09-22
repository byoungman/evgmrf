#' Colorado daily precipitation annual order statistics
#'
#' The five largest daily precipitation values in each year, for each cell of
#' a 28 by 16 longitude/latitude grid covering Colorado, for 56 years
#' (1970-2025), together with the elevation at each grid cell. The order
#' statistics are suitable for fitting r-largest order statistics models with
#' \code{\link{evgmrf}}.
#'
#' The precipitation data are simulated, and taken from the North American
#' CORDEX (NA-CORDEX) archive. They come from the historical experiment of the
#' WRF regional climate model driven by the GFDL-ESM2M global climate model, on
#' the NAM-22i grid (approximately 0.25 degree resolution, regular
#' longitude/latitude). The raw (not bias-corrected) daily precipitation
#' output was used, assuming a 365-day calendar. The grid is the rectangle of
#' grid cells spanning the cells whose centres lie in Colorado, found using
#' \code{maps::map.where}; no cells are set to \code{NA}.
#'
#' Elevation was obtained by bilinear interpolation of a 1 km global elevation
#' data set onto the centres of the grid cells, so it is a point value at each
#' cell centre rather than a cell average.
#'
#' @format A list with four components:
#' \describe{
#'   \item{x}{Numeric vector of length 28. Longitudes (degrees east) of the
#'     grid columns.}
#'   \item{y}{Numeric vector of length 16. Latitudes (degrees north) of the
#'     grid rows.}
#'   \item{prcp}{Numeric array with dimensions \code{c(56, 5, 28, 16)}, indexed
#'     by year, rank, longitude and latitude. \code{prcp[t, k, i, j]} is the
#'     \code{k}th largest daily precipitation in year \code{1969 + t} at
#'     longitude \code{x[i]} and latitude \code{y[j]}, so that
#'     \code{prcp[, 1, , ]} holds the annual maxima. Units are those of the
#'     NA-CORDEX \code{prec} variable (kg m-2 s-1 in the raw files; multiply
#'     by 86400 for mm/day).}
#'   \item{elevation}{Numeric matrix with dimensions \code{c(28, 16)}.
#'     \code{elevation[i, j]} is the elevation in metres above sea level at
#'     longitude \code{x[i]} and latitude \code{y[j]}.}
#' }
#'
#' @source NA-CORDEX archive, file
#'   \code{prec.hist.GFDL-ESM2M.WRF.day.NAM-22i.raw.nc}.
#'   \url{https://na-cordex.org}
#'   
#' @references
#' Mearns, L. O., et al. (2017). The NA-CORDEX dataset, version 1.0.
#' NCAR Climate Data Gateway. \doi{10.5065/D6SJ1JCH}
#'
#' @docType data
#' @keywords datasets
#' @name COorder
#' @usage data(COorder)
#'
#' @examples
#' library(evgmrf)
#' data(COorder)
#'
#' # Mean annual maximum of daily precipitation
#' COorder$z <- colMeans(COorder$prcp[, 1, , ])
#' image(COorder)
#'
#' # Elevation on the same grid, and its relationship with the mean annual
#' # maximum
#' image(COorder$x, COorder$y, COorder$elevation)
#' plot(as.vector(COorder$elevation), as.vector(COorder$z),
#'      xlab = "Elevation (m)", ylab = "Mean annual maximum")
#'
#' # r-largest model using the top five values in each year
#' m_r5 <- evgmrf(COorder$prcp, family = "rlarge")
#' plot(m_r5)
NULL

#' Colorado daily precipitation above 20mm
#'
#' A list of grid longitude (\code{x}) and latitude (\code{y}) values and an 
#' array (\code{prcp}) of daily precipitation amounts above 20mm (\code{threshold})
#' for 56 years (dimension 1), for a 28 by 16 grid.
#'
#' @format A list of length three
#' 
#' @docType data
#' @keywords datasets
#' @name COexc
#' @usage data(COexc)
#' @examples
#' 
#' library(evgmrf)
#' data(COexc)
#'
#' COexc$z <- colMeans(COexc$prcp, na.rm = TRUE)
#' image(COexc)
#'
NULL

#' Colorado daily precipitation top 50
#'
#' The 50 largest daily precipitation values over the 56-year simulation
#' period (1970-2025), stored in decreasing order, for each cell of a 28 by 16
#' longitude/latitude grid covering Colorado. They are suitable for fitting
#' peaks-over-threshold models based on the largest values with
#' \code{\link{evgmrf}}. The grid is the same as that of \code{\link{COorder}}.
#'
#' The precipitation data are simulated, and taken from the North American
#' CORDEX (NA-CORDEX) archive. They come from the historical experiment of the
#' WRF regional climate model driven by the GFDL-ESM2M global climate model, on
#' the NAM-22i grid (approximately 0.25 degree resolution, regular
#' longitude/latitude). The raw (not bias-corrected) daily precipitation
#' output was used. The grid is the rectangle of grid cells spanning the cells
#' whose centres lie in Colorado, found using \code{maps::map.where}; no cells
#' are set to \code{NA}.
#'
#' @format A list with three components:
#' \describe{
#'   \item{x}{Numeric vector of length 28. Longitudes (degrees east) of the
#'     grid columns.}
#'   \item{y}{Numeric vector of length 16. Latitudes (degrees north) of the
#'     grid rows.}
#'   \item{prcp}{Numeric array with dimensions \code{c(50, 28, 16)}, indexed
#'     by rank, longitude and latitude. \code{prcp[k, i, j]} is the
#'     \code{k}th largest daily precipitation over the whole period at
#'     longitude \code{x[i]} and latitude \code{y[j]}. Units are those of the
#'     NA-CORDEX \code{prec} variable (kg m-2 s-1 in the raw files; multiply
#'     by 86400 for mm/day).}
#' }
#'
#' @source NA-CORDEX archive, file
#'   \code{prec.hist.GFDL-ESM2M.WRF.day.NAM-22i.raw.nc}.
#'   \url{https://na-cordex.org}
#'
#' @references
#' Mearns, L. O., et al. (2017). The NA-CORDEX dataset, version 1.0.
#' NCAR Climate Data Gateway. \doi{10.5065/D6SJ1JCH}
#'
#' @docType data
#' @keywords datasets
#' @name COtop50
#' @usage data(COtop50)
#'
#' @examples
#' library(evgmrf)
#' data(COtop50)
#'
#' # Mean of the top 50 daily values
#' COtop50$z <- colMeans(COtop50$prcp)
#' image(COtop50)
#'
#' # Fit the model to the full array of top-50 values
#' m_co <- evgmrf(COtop50$prcp, family = "poisgpd",
#'                args = list(r = 50, nper = 56))
#' plot(m_co)
NULL

#' Colorado county annual maxima of daily precipitation
#'
#' Annual maxima of daily precipitation for the 64 counties of Colorado over
#' 56 years, together with the county boundaries and the adjacency structure
#' of the counties, for use with \code{\link{evgmrf}}.
#'
#' The precipitation data are simulated, and taken from the North American
#' CORDEX (NA-CORDEX) archive. They come from the historical experiment of the
#' WRF regional climate model driven by the GFDL-ESM2M global climate model, on
#' the NAM-22i grid (approximately 0.25 degree resolution). The raw (not
#' bias-corrected) daily precipitation output was used, assuming a 365-day
#' calendar. For each year (1970-2025) the maximum of the daily values was
#' found at every grid cell whose centre lies in Colorado. Each such grid cell
#' was then assigned to a county using \code{maps::map.where}, and the value
#' for a county is the largest annual maximum over the grid cells assigned to
#' it. County boundaries are from the \pkg{maps} package. Two counties are
#' adjacent if their boundaries share at least one point (queen contiguity, as
#' found by \code{spdep::poly2nb}).
#'
#' @format A list with three components:
#' \describe{
#'   \item{prcp}{Numeric matrix with 56 rows and 64 columns. Row \code{i} is
#'     year \code{1969 + i} (1970 to 2025) and the columns, which are named,
#'     are the counties. Units are those of the NA-CORDEX \code{prec} variable
#'     (kg m-2 s-1 in the raw files; multiply by 86400 for mm/day).}
#'   \item{polygons}{Named list of length 64, one element per county, in the
#'     same order as the columns of \code{prcp}. Each element is a list with
#'     components \code{x} (longitude) and \code{y} (latitude) giving the
#'     coordinates of the county boundary in degrees.}
#'   \item{adjacency}{Sparse 64 x 64 symmetric matrix (a \pkg{Matrix}
#'     \code{"CsparseMatrix"}) with rows and columns in the same order as
#'     \code{polygons} and the columns of \code{prcp}. Element \eqn{(i, j)} is
#'     1 if counties \eqn{i} and \eqn{j} are adjacent and 0 otherwise, and the
#'     diagonal is zero. It can be passed as \code{W} to \code{\link{evgmrf}}.}
#' }
#'
#' @source NA-CORDEX archive, file
#'   \code{prec.hist.GFDL-ESM2M.WRF.day.NAM-22i.raw.nc}.
#'   \url{https://na-cordex.org}. County boundaries from
#'   \code{maps::map("county", "colorado")}.
#'
#' @references
#' Mearns, L. O., et al. (2017). The NA-CORDEX dataset, version 1.0.
#' NCAR Climate Data Gateway. \doi{10.5065/D6SJ1JCH}
#'
#' @docType data
#' @keywords datasets
#' @name COcnty
#' @usage data(COcnty)
#'
#' @examples
#' library(evgmrf)
#' data(COcnty)
#'
#' # Map of the mean annual maximum by county
#' COcnty_prcp_mean <- colMeans(COcnty$prcp)
#' coords <- do.call(rbind, lapply(COcnty$polygons,
#'                                 function(x) rbind(cbind(x$x, x$y), NA)))
#' plot(coords, type = "n")
#' cols <- rev(grey(ppoints(COcnty_prcp_mean)))[rank(COcnty_prcp_mean)]
#' polygon(coords, col = cols)
#'
#' # GEV model with a GMRF over the counties, defined by the adjacency matrix
#' m_poly <- evgmrf(COcnty$prcp, family = "gev", W = COcnty$adjacency)
#' plot(m_poly, polygons = COcnty$polygons)
NULL

#' Top 50 daily precipitation values over Washington State
#'
#' Simulated daily precipitation extremes over Washington State, taken from
#' the North American CORDEX (NA-CORDEX) archive. For each grid cell, the 50
#' largest daily precipitation values over the 56-year simulation period
#' (1970-2025) are stored in decreasing order. Grid cells that fall outside
#' Washington State are set to \code{NA}.
#'
#' The data come from the historical experiment of the WRF regional climate
#' model driven by the GFDL-ESM2M global climate model, on the NAM-22i grid
#' (approximately 0.25 degree resolution, regular longitude/latitude). The
#' raw (not bias-corrected) precipitation output was used.
#'
#' @format A list with three components:
#' \describe{
#'   \item{x}{Numeric vector of length 31. Longitudes (degrees east) of the
#'     grid columns, ranging from -124.625 to -117.125.}
#'   \item{y}{Numeric vector of length 14. Latitudes (degrees north) of the
#'     grid rows, ranging from 45.625 to 48.875.}
#'   \item{prcp}{Numeric array with dimensions \code{c(50, 31, 14)}, indexed
#'     by rank, longitude and latitude. \code{prcp[k, i, j]} is the
#'     \code{k}th largest daily precipitation at longitude \code{x[i]} and
#'     latitude \code{y[j]}. Cells outside Washington State are \code{NA}.
#'     Units are those of the NA-CORDEX \code{prec} variable
#'     (kg m-2 s-1 in the raw files; multiply by 86400 for mm/day).}
#' }
#'
#' @source NA-CORDEX archive, file
#'   \code{prec.hist.GFDL-ESM2M.WRF.day.NAM-22i.raw.nc}.
#'   \url{https://na-cordex.org}
#'
#' @references
#' Mearns, L. O., et al. (2017). The NA-CORDEX dataset, version 1.0.
#' NCAR Climate Data Gateway. \doi{10.5065/D6SJ1JCH}
#'
#' @docType data
#' @keywords datasets
#' @name WAprcp
#' @usage data(WAprcp)
#'
#' @examples
#' library(evgmrf)
#' data(WAprcp)
#'
#' # Fit the model to the full array of top-50 values
#' m_wa1 <- evgmrf(WAprcp$prcp, family = "poisgpd",
#'                 args = list(r = 50, nper = 56))
#' plot(m_wa1, set2NA = TRUE)
#'
#' # Fit the model using only the grid cells with data, with an explicit
#' # neighbourhood matrix W
#' some_prcp <- apply(is.finite(WAprcp$prcp), 2:3, any)
#' some_prcp_index <- which(some_prcp, arr.ind = TRUE)
#' W <- index2W(some_prcp_index)
#' WAprcp_list <- apply(some_prcp_index, 1,
#'                      function(idx) WAprcp$prcp[, idx[1], idx[2]],
#'                      simplify = FALSE)
#' m_wa2 <- evgmrf(WAprcp_list, family = "poisgpd", W = W,
#'                 args = list(r = 50, nper = 56))
NULL