#' Essentially a wrapper for spatstat's ppm function; fits a
#' special case of the Hieararchical Strauss model that assumes
#' no interaction between points of the same type.
#'
#' @param p spatstat ppp object with two qualitative mark levels, 1 and 2
#' @param r radius of interaction for model fitting
#' @param quad.spacing space between points in quadrature used in estimation
#' @param correction correction used in fitting; see spatstat documentation of
#' "ppm.ppp" for more details
#' @return ppm object

FitHSPPPFull = function(p, r = 30, quad.spacing = 30, correction = "Ripley"){
  mat = matrix(c(NA, r, r, NA), nrow = 2)
  mod = ppm(p ~ marks,
            HierStrauss(radii=mat, archy = c(1,2)),
            eps = quad.spacing,
            correction = correction)

  return(mod)
}
