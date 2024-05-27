## calculate distrance between two points
len.py <- function(x1, x2, y1, y2) {
  a <- abs((x1-x2))
  b <- abs((y1-y2))
  c <- sqrt(a^2 + b^2) 
  return(c)
}