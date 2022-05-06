#' Determines if a set of points is a hole in spatstat
#' @param x vector of x-coordinates of points
#' @param y vector of y-coordinates of points
#' @return boolean
#' @noRd
IsHole = function(x, y){
  testwindow = try(owin(poly = list(x = x, y = y)),
                   silent = T)

  if(class(testwindow) == "try-error"){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
