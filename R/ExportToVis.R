#' Given an object output from the TessellateBiopsy function, formats it and
#' outputs it into a `.json` file that can be visualized using the companion
#' web application:
#'
#' https://nateosher.github.io/SPARTIN.html
#'
#' @param tes Output from the TessellateBiopsy function- the data you want to
#' visualize.
#' @param val_list List of values associated with each tile you'd like to
#' visualize. Should be a list the same length as tes$tiles, such that each
#' list entry has the same set of attributes each with associated values that
#' you'd like to visulize.
#' @param path Path where resulting .json file should be written to
#' @import jsonlite
#' @import purrr
#' @export
ExportToVis = function(tes, val_list, path, write=TRUE){
  if(length(tes$tiles != length))
    stop("value list should be same length as tile list")

  final_list = vector('list', length(tes$tiles))
  for(i in 1:length(tes$tiles)){
    final_list[[i]]$polygons = map(tes$tiles[[i]]$window$bdry, function(b){
      c(map2(b$x, b$y, ~ c(.x, .y)), list(c(b$x[1], b$y[1])))
    })
    final_list[[i]]$is_hole = map_lgl(tes$tiles[[i]]$window$bdry, function(b){
      IsHole(b$x, b$y)
    })
    final_list[[i]]$vals = val_list[[i]]
    final_list[[i]]$id = i
    final_list[[i]]$point_data = map(1:length(tes$tiles[[i]]$x), function(j){
      c(tes$tiles[[i]]$x[j], tes$tiles[[i]]$y[j], tes$tiles[[i]]$marks[j])
    })
  }

  output = toJSON(final_list, auto_unbox = TRUE)


  sink(path)
  cat(output)
  sink()
}
