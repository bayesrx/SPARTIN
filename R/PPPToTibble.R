#' Converts a marked point process (ppp) object to a tibble with
#' three columns: CentroidX, CentroidY, and Mark
#'
#' @param p marked point process
#' @return tibble with columns "CentroidX", "CentroidY", and "Mark"
#' @noRd
PPPToTibble = function(p){
  tib = tibble(CentroidX = p$x, CentroidY = p$y, Mark = p$marks)
  return(tib)
}
