#' Plots a tessellation, colored according to the value vector passed. See
#' "SPARTIN-General-Usage" vignette for example.
#' @param tes Tessellation from `TessellateBiopsy` function.
#' @param val Vector of values according to which the tiles should be colored.
#' Should be the same length as `tes$tiles`, i.e. one for each tile.
#' @param val_name String; this is what the legend will be labeled in the plot.
#' purely aesthetic.
#' @return gg/ggplot object.
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @export
PlotTessellation = function(tes, val, val_name = ""){
  tes_tib = tibble(x = numeric(0), y = numeric(0), ids = numeric(0), val = numeric(0))

  cur_id = 1

  if(length(tes$tiles) != length(val))
    stop("length of 'val' must be same as number of tiles in tessellation")

  for(i in 1:length(tes$tiles)){
    for(j in 1:length(tes$tiles[[i]]$window$bdry)){
      cur_bdry = tes$tiles[[i]]$window$bdry[[j]]

      cur_n = length(cur_bdry$x)

      if(IsHole(cur_bdry$x, cur_bdry$y)){
        cur_val = NA
      }else{
        cur_val = val[i]
      }

      cur_tib = tibble(x = cur_bdry$x, y = cur_bdry$y,
                       ids = cur_id, val = cur_val)


      tes_tib = bind_rows(tes_tib, cur_tib)
      cur_id = cur_id + 1
    }
  }

  plt = ggplot(tes_tib) +
    geom_polygon(aes(x = x, y = y, group = ids, fill = val)) +
    labs(fill = val_name) +
    scale_fill_continuous(na.value = "white") +
    theme_bw()

  return(plt)
}
