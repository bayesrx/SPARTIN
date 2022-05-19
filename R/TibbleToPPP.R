#' Converts a tibble with columns CentroidX, CentroidY,
#' and Mark into a point process (ppp) object
#' with a square window determined by the data.
#'
#' @param t A tibble with columns CentroidX, CentroidY, and Mark
#' @return an object of class "ppp"
#' @import spatstat
#' @noRd
TibbleToPPP = function(t){
  square_ppp = ppp(x = t$CentroidX, y = t$CentroidY, marks = t$Mark,
                   window = owin(c(min(t$CentroidX), max(t$CentroidX)),
                                 c(min(t$CentroidY), max(t$CentroidY))))

  return(square_ppp)
}
