#' CDF_trapz
#'
#' Computing the area under curve using trapezoidal method
#'
#' @param f A function of curve.
#' @param a A starting point value.
#' @param b An ending point value.
#' @param n.pts A number of subdivisions of trapezoid shape under the curve.
#' @returns A value.
#' @examples
#' f <- exp(-x^2)
#' cdf_trapz(f=f, a=0, b=1, n.pts=15)

cdf_trapz <- function(f,a=0,b,n.pts=1001){
  h <- (b-a)/n.pts
  x.vec <- seq(a,b,by=h)
  f.vec <- sapply(x.vec,f)
  return(h*(f.vec[1]/2 + sum(f.vec[2:n.pts]) + f.vec[n.pts+1]/2))
}
