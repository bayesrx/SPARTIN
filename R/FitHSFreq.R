#' Essentially a wrapper for spatstat's ppm function; fits a
#' special case of the Hieararchical Strauss model that assumes
#' no interaction between points of the same type.
#'
#' @param p spatstat ppp object with two qualitative mark levels, 1 and 2.
#' @param r radius of interaction for model fitting.
#' @param quad.spacing space between points in quadrature used in estimation.
#' @param correction correction used in fitting; see spatstat documentation of
#' "ppm.ppp" for more details.
#' @return ppm object.
#' @import spatstat.model
#' @export
FitHSFreq = function(p, r, quad.spacing, correction = "Ripley"){
  mat = matrix(c(NA, r, r, NA), nrow = 2)
  mod = spatstat.model::ppm(p ~ marks,
            spatstat.model::HierStrauss(radii=mat, archy = c(1,2)),
            eps = quad.spacing,
            correction = correction)

  return(mod)
}
