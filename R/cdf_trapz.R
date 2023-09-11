cdf_trapz <- function(f,a=0,b,n.pts=1001){
  h <- (b-a)/n.pts
  x.vec <- seq(a,b,by=h)
  f.vec <- sapply(x.vec,f)
  return(h*(f.vec[1]/2 + sum(f.vec[2:n.pts]) + f.vec[n.pts+1]/2))
}