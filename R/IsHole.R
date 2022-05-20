#' Determines if a set of points is a hole in spatstat
#' @param x vector of x-coordinates of points
#' @param y vector of y-coordinates of points
#' @return boolean
#' @import spatstat
#' @noRd
IsHole = function(x, y){
  if(length(x) != length(y))
    stop("length of x must equal length of y")

  if(!(is.numeric(x)) || !(is.numeric(y)))
    stop("x and y must be numeric vectors")

  if(length(x) < 3){
    warning("length of point list is less than three")
    return(FALSE)
  }



  testwindow = try(owin(poly = list(x = x, y = y)),
                   silent = T)

  if(class(testwindow) == "try-error"){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
