#' Converts a marked point process (ppp) object to a tibble with
#' three columns: CentroidX, CentroidY, and PredictedClass
#'
#' @param p marked point process
#' @return tibble with columns "CentroidX", "CentroidY", and "PredictedClass"
PPPToTibble = function(p){
  tib = tibble(CentroidX = p$x, CentroidY = p$y, PredictedClass = p$marks)
  return(tib)
}
