#' Reducing array dimensions to drop NA-filled dimensions
#'
#' Function \code{tighten} takes a three-dimensional array with NAs in, and reduces 
#' dimension, if possible.
#'
#' @param x an array, or matrix if W supplied, or array if W supplied and r-largest
#' 
#' @details
#' 
#' To do.
#' 
#' @examples
#'
#' # To follow
#'
#' @return An \code{array}
#' 
#' @export
#' 
tighten <- function(x) {
  x <- apply(x, 2:3, sort, decreasing = TRUE, na.last = TRUE)
  end <- max(apply(x, 2:3, function(x) min(which(!is.finite(x)))))
  x[1:end, , ]
}
