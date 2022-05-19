#' A basic function for simulating data
#'
#' @param n1 Number of points of type 1 to be generated.
#' @param n2 Number of points of type 2 to be generated.
#' @param phi Interaction parameter, a real number in [-1, 1]. -1 indicates most
#' negative Interaction; 1 indicates most positive interaction.
#' @param winX Horizontal width of the simulation window.
#' @param winY Vertical height of the simulation window.
#' @param r Radius of interaction; only necessary for phi > 0, i.e. positive
#' interaction between points.
#' @import purrr
#' @import spatstat
#' @export
SimulateData = function(n1, n2, phi, winX, winY, r = NULL){
  if(phi < -1 || phi > 1)
    stop("phi must be in [-1, 1]")

  if(phi <= 0){
    # Negative Simulation
    t1.x = runif(n1, 0, winX/2*(2 + phi))
    t1.y = runif(n1, 0, winY)

    t2.x = runif(n2, winX/2*(-1 * phi), winX)
    t2.y = runif(n2, 0, winY)

    final.x = c(t1.x, t2.x)
    final.y = c(t1.y, t2.y)
    final.marks = factor(c(rep(1, n1), rep(2, n2)))

    final.pp = ppp(x = final.x, y = final.y, marks = final.marks,
                   window = owin(c(0, winX), c(0, winY)))

  }else{
    if(is.null(r))
      stop("r must be provided for positive simulations")
    t1.x = runif(n1, 0, winX)
    t1.y = runif(n1, 0, winY)


    # Generate "nearness" indicators, i.e. whether each til will be near tumor cell
    # or randomly generated
    t2.near.t1 = rbernoulli(n2, p = phi)

    # Random til x values
    t2.not.near.x = runif(sum(!t2.near.t1), 0, winX)
    t2.not.near.y = runif(sum(!t2.near.t1), 0, winY)

    # Which tumor cell should each til be near (if applicable)
    t1.centers = sample(1:n1, sum(t2.near.t1), replace = T)
    # Generate offsets from tumor centers within radius r, using
    # polar coordinates
    t2.near.r = runif(sum(t2.near.t1), 0, r)
    til.near.theta = runif(sum(t2.near.t1), 0, 2*pi)
    # Convert back to normal x/y coords
    til.near.x = t2.near.r*cos(til.near.theta) + t1.x[t1.centers]
    til.near.y = t2.near.r*sin(til.near.theta) + t1.y[t1.centers]

    final.x = c(t1.x, t2.not.near.x, til.near.x)
    final.y = c(t1.y, t2.not.near.y, til.near.y)
    final.marks = factor(c(rep(1, n1), rep(2, n2)))

    final.pp = ppp(x = final.x, y = final.y, marks = final.marks,
                   window = owin(c(0, winX), c(0, winY)))
  }

  return(final.pp)
}
