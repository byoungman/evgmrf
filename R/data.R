#' Colorado daily precipitation annual order statistics
#'
#' A list of grid longitude (\code{x}) and latitude (\code{y}) values and an 
#' array (\code{prcp}) of precipitation values for 56 years (dimension 1), the top
#' five order statistics (dimension 2) from 1970 to 2025 for a 28 by 16 grid.
#'
#' @format A list of lenght three
#' 
#' @docType data
#' @keywords datasets
#' @name COtop5prcp
#' @usage data(COtop5prcp)
#' @examples
#' 
#' library(evgmrf)
#' data(COtop5prcp)
#'
#' COtop5prcp$z <- colMeans(COtop5prcp$prcp, dims = 2)
#' image(COtop5prcp)
#'
NULL

#' Colorado daily precipitation above 20mm
#'
#' A list of grid longitude (\code{x}) and latitude (\code{y}) values and an 
#' array (\code{prcp}) of daily precipitation amounts above 20mm for 56 years 
#' (dimension 1), for a 28 by 16 grid.
#'
#' @format A list of length three
#' 
#' @docType data
#' @keywords datasets
#' @name COprcp20
#' @usage data(COprcp20)
#' @examples
#' 
#' library(evgmrf)
#' data(COprcp20)
#'
#' COprcp20$z <- colMeans(COprcp20$prcp, na.rm = TRUE)
#' image(COprcp20)
#'
NULL

#' Colorado daily precipitation top 50
#'
#' A list of grid longitude (\code{x}) and latitude (\code{y}) values and an 
#' array (\code{prcp}) of top 50 daily precipitation amounts for 56 years 
#' (dimension 1), for a 28 by 16 grid.
#'
#' @format A list of length three
#' 
#' @docType data
#' @keywords datasets
#' @name COprcp50
#' @usage data(COprcp50)
#' @examples
#' 
#' library(evgmrf)
#' data(COprcp50)
#'
#' COprcp50$z <- colMeans(COprcp50$prcp, na.rm = TRUE)
#' image(COprcp50)
#'
NULL

#' Colorado county annual maxima of daily precipitation
#'
#' Three data files. \code{COcnty.prcp} a 56 x 64 matrix of 56 annual maxima of
#' daily precipitation for the 64 Colorado counties. \code{COcnty.polygons} a list
#' of 64 polygons giving the county boundaries as a list with longitude (\code{x}) 
#' and latitude (\code{y}) coordinates. \code{COcnty.adjacency} a 64 x 64 
#' adjacency matrix for with rows and columns matching the order of 
#' \code{COcnty.polygons} and the columns of \code{COcnty.prcp}.
#'
#' @format A matrix (\code{COcnty.prcp}), list \code{COcnty.polygons} and a matrix
#' \code{COcnty.adjacency}.
#' 
#' @docType data
#' @keywords datasets
#' @name COcnty
#' @aliases COcnty.adjacency COcnty.polygons COcnty.prcp
#' @usage data(COcnty)
#' @examples
#' 
#' library(evgmrf)
#' data(COcnty)
#'
#' COcnty.prcp_mean <- colMeans(COcnty.prcp)
#' coords <- do.call(rbind, lapply(COcnty.polygons, function(x) rbind(cbind(x$x, x$y), NA)))
#' plot(coords, type = 'n')
#' cols <- rev(grey(ppoints(COcnty.prcp_mean)))[rank(COcnty.prcp_mean)]
#' polygon(coords, col = cols)
#'
NULL

#' Colorado excess precipitation data
#'
#' Description of what COexcprcp contains goes here.
#'
#' @format Describe the format (e.g., A data frame with X rows and Y variables)
#' @docType data
#' @keywords datasets
#' @name COexcprcp
#' @usage data(COexcprcp)
NULL
