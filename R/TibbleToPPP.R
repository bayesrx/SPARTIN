#' Converts a tibble with columns CentroidX, CentroidY,
#' and PredictedClass into a point process (ppp) object
#' with a square window determined by the data.
#' Points with Class = 1 are tumor cells, points with
#' Class = 2 are immune cells.
#'
#' @param t A tibble with columns CentroidX, CentroidY, and PredictedClass
#'
#' @return an object of class "ppp"
TibbleToPPP = function(t){
  square_ppp = ppp(x = t$CentroidX, y = t$CentroidY, marks = t$PredictedClass,
                   window = owin(c(min(t$CentroidX), max(t$CentroidX)),
                                 c(min(t$CentroidY), max(t$CentroidY))))

  return(square_ppp)
}
