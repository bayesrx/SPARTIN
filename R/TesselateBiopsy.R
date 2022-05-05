#' Tessellates a biopsy into sub-windows that the models are ultimately fit on.
#'
#' @param full.tib Full tibble representing all cells in the biopsy to be
#' tessellated. Should have three columns: CentroidX, CentroidY, and PredictedClass
#' @param sigma Radius of kernel smoothing used in intensity thresholding. Larger
#' radius means less white space.
#' @param eps Size of "pixel" that window is divided into for intensity
#' thresholding. Note that this parameter tends to be the "bottleneck"- it's
#' best to start larger, and see how small a value is necessary and/or
#' computationally feasible.
#' @param threshold Minimum value of pixels to keep in intensity threshold. Too
#' large of a value will result in the tesselation being "choppy," i.e. too
#' much of the space will be thresholded out, and too small of a value
#' will result in a tesselation with too much white space.
#' @param clust.size Target number of type 1 points in each cluster. In all
#' likelihood this won't be exactly satisfied.
#' @param max.clust.size Maximum number of type 1 points allowed in a cluster.
#' clusters with more than this number of points will be recursively split
#' into smaller clusters until they have fewer points than this bound.
#' @param min.clust.size Minimum number of type 1 points allowed in a cluster.
#' Clusters with fewer than this number of points will not be included in the
#' final partition.
#' @param enforce.contiguity Determines whether each resulting tile should be
#' contiguous. If true, keeps largest contiguous portion of each tile,
#' provided it has enough points to meet the threshold.
#' @param clust.method Determines how dirichlet tesselation
#' tiles are grouped into larger tiles. Must be either "greedy" or "kmeans."
#' While both work, "kmeans" is strongly recommended at the moment.
#' @param progress Boolean; if "TRUE," prints updates to the console.
#' @return a list with four attributes. "tiles" is a list where each entry is
#' a spatstat ppp object corresponding to each tile the space has been
#' partitioned into. "pp" is a spastat ppp object corresponding to the overall
#' intensity thresholded space and points. "pixels" is a matrix of the "pixels"
#' that the window is divided into for intensity thresholding, as well as the
#' values assigned to each pixel by the kernel smoothing. "window" is the
#' overall window that results from the tesselation; this is equivalent to
#' $pp$window.

TesselateBiopsy = function(full.tib, sigma, eps,
                           threshold, clust.size,
                           max.clust.size = NULL, min.clust.size = NULL,
                           enforce.contiguity = TRUE,
                           clust.method = "kmeans",
                           progress = TRUE){

  if(is.null(max.clust.size)) max.clust.size = max.clust.size
  if(is.null(min.clust.size)) min.clust.size = 1

  if(progress) print("Intensity thresholding...")

  pp = TibbleToPPP(full.tib)
  tk = density.ppp(pp, sigma = sigma, eps = eps)

  test.v = apply(tk$v, 2, rev)

  pos.int = which(tk$v > threshold, arr.ind = T)

  x.min = min(pp$x)
  y.min = min(pp$y)

  x.step = tk$xstep
  y.step = tk$ystep

  pos.r = round(pos.int[,1])
  pos.c = round(pos.int[,2])

  final.win = owin(xrange = c(x.min + (pos.c[1] - 1)*x.step, x.min + pos.c[1]*x.step),
                   yrange = c(y.min + (pos.r[1] - 1)*y.step, y.min + pos.r[1]*y.step))

  final.win.list = vector('list', length(pos.r) + 1)
  for(i in 1:length(pos.r)){
    cur.win = owin(xrange = c(x.min + (pos.c[i] - 1)*x.step, x.min + pos.c[i]*x.step),
                   yrange = c(y.min + (pos.r[i] - 1)*y.step, y.min + pos.r[i]*y.step))
    # final.win = union.owin(final.win, cur.win)
    final.win.list[[i + 1]] = cur.win
  }

  if(progress) print("Constructing window...")
  final.win = do.call(union.owin, final.win.list)

  final.pp = ppp(full.tib$CentroidX, full.tib$CentroidY,
                 marks = factor(full.tib$PredictedClass),
                 window = final.win)

  tum.set = full.tib %>%
    filter(PredictedClass == 1)

  tums.pp = ppp(tum.set$CentroidX, tum.set$CentroidY,
                window = final.win)

  # Applying dirichlet tesselation
  if(progress) print("Applying dirichlet tesselation...")
  dirch.tess = dirichlet(tums.pp)

  # Seeing if sorting is necessary
  xy.mat = matrix(c(tums.pp$x, tums.pp$y), ncol = 2)
  dirch.points = dirch.tess$tiles
  dirch.points = unname(dirch.points)

  full.max.x = max(tums.pp$window$xrange)
  full.max.y = max(tums.pp$window$yrange)

  if(progress) print("Creating clusters...")
  finished = FALSE
  safety = 1
  max.iter = 10000
  # TODO: edge case where number of tumor cells is divisible by cluster size
  full.window = owin(c(0, 100000), c(0, 100000))
  if(clust.method %in% c("greedy", "Greedy")){
    pointset.list = vector('list', length = ceiling(tums.pp$n/clust.size))
    while(!finished && safety < max.iter){
      # print(paste("Tile", safety, "out of", length(pointset.list)))
      # Pick starting point as whichever one has smallest x-coordinate
      start.index = which.min(xy.mat[1,])
      cur.point = xy.mat[start.index,]
      cur.x = xy.mat[,1]
      cur.y = xy.mat[,2]

      if(length(dirch.points) >= clust.size){
        nn.test = FindNeighborDists(cur.point, cur.x, cur.y)
        subgroup.indices = which(nn.test %in% sort(nn.test)[1:clust.size])
        subgroup.x = cur.x[subgroup.indices]
        subgroup.y = cur.y[subgroup.indices]
      }

      else{
        subgroup.indices = 1:length(dirch.points)
        subgroup.x = cur.x
        subgroup.y = cur.y
      }

      cur.xranges = unlist(map(dirch.points[subgroup.indices], ~{.x$xrange}))
      # old \/
      # cur.xranges = unlist(map(dirch.points[subgroup.indices], ~{.x$xrange}),
      #                     use.names = F)
      maxx = max(cur.xranges)
      minx = min(cur.xranges, 0)

      cur.yranges = unlist(map(dirch.points[subgroup.indices], ~{.x$yrange}))
      # old \/
      # cur.yranges = unlist(map(dirch.points[subgroup.indices], ~{.x$yrange}),
      #                      use.names = F)
      maxy = max(cur.yranges)
      miny = min(cur.yranges, 0)


      tile.subset = dirch.points[subgroup.indices]
      cur.win = do.call(union.owin, tile.subset)
      # cur.win = intersect.owin(cur.win, full.window)

      # tile.subset = map(tile.subset, ~{.x$bdry})
      # tile.subset = unlist(tile.subset, recursive = F)

      dirch.points = dirch.points[-subgroup.indices]

      xy.mat = xy.mat[-subgroup.indices,]


      # cur.win = owin(xrange = c(minx, maxx),
      #                yrange = c(miny, maxy),
      #                poly = tile.subset)


      cur.pointset = ppp(x = pp$x, y = pp$y,
                         marks = pp$marks,
                         window = cur.win)

      cur.pointset = as.ppp(cur.pointset)

      pointset.list[[safety]] = cur.pointset

      safety = safety + 1
      if(is.na(nrow(xy.mat)) || is.null(nrow(xy.mat)) || nrow(xy.mat) == 0){
        finished = TRUE
      }
      if(safety >= max.iter){
        print("Something went badly wrong")
      }
    }
  }
  else if(clust.method %in% c("kmeans", "k", "Kmeans")){
    points.km = kmeans(xy.mat, centers = floor(nrow(xy.mat)/clust.size))
    xy.km = cbind(xy.mat, points.km$cluster)

    cluster.counts = table(points.km$cluster)
    clusters = as.numeric(names(cluster.counts))
    names(cluster.counts) = NULL

    too.large = which(cluster.counts > max.clust.size)

    max.clust.val = max(clusters)

    while(length(too.large) > 0){
      # print(length(too.large))
      large.clusters = clusters[too.large]
      for(i in 1:length(large.clusters)){
        cur.submat.rows = which(as.numeric(xy.km[,3]) == large.clusters[i])
        cur.submat = xy.km[cur.submat.rows,]
        new.km = kmeans(cur.submat, centers = 2)
        cur.km.clusters = new.km$cluster
        # Number the new clusters according to the overall matrix; assign
        # points in cluster 1 to the old cluster index, and in cluster 2
        # to the current maximum + 1
        cur.km.clusters = map_dbl(cur.km.clusters, ~{ifelse(.x == 1,
                                                            large.clusters[i],
                                                            max.clust.val + 1)})

        max.clust.val = max.clust.val + 1
        xy.km[as.numeric(xy.km[,3]) == large.clusters[i], 3] = cur.km.clusters
      }

      # Check for oversize clusters
      cluster.counts = table(as.numeric(xy.km[,3]))
      clusters = as.numeric(names(cluster.counts))
      names(cluster.counts) = NULL

      too.large = which(cluster.counts > max.clust.size)
    }

    if(length(too.large > 0)) stop("Too large > 0")
    if(any(cluster.counts > max.clust.size)) stop("cluster.counts > max.size")

    # return(xy.km)


    # SHOULD BE xy.km, not points.km$cluster- FIX
    pointset.list = vector('list', length = max(points.km$cluster))
    for(i in 1:max(as.numeric(xy.km[,3]))){
      subgroup.indices = which(points.km$cluster == i)
      # TODO: Use to throw out small tiles?
      if(length(subgroup.indices) < min.clust.size) next
      cur.xranges = unlist(map(dirch.points[subgroup.indices], ~{.x$xrange}))
      maxx = max(cur.xranges)
      minx = min(cur.xranges, 0)

      cur.yranges = unlist(map(dirch.points[subgroup.indices], ~{.x$yrange}))
      maxy = max(cur.yranges)
      miny = min(cur.yranges, 0)


      tile.subset = dirch.points[subgroup.indices]
      cur.win = do.call(union.owin, tile.subset)

      # To get rid of extraneous borders
      cur.win = intersect.owin(cur.win, full.window)

      # If windows should be contiguous, keep sub-tile with the most points
      if(enforce.contiguity && length(cur.win$bdry) > 1){
        sub.tile.point.counts = numeric(length(cur.win$bdry))

        for(ec in 1:length(sub.tile.point.counts)){
          if(IsHole(cur.win$bdry[[ec]]$x, cur.win$bdry[[ec]]$y)){
            sub.tile.point.counts[ec] = -1
            next
          }
          cur.sub.tile.window = owin(poly = list(x = cur.win$bdry[[ec]]$x,
                                                 y = cur.win$bdry[[ec]]$y))

          cur.sub.tile = ppp(x = pp$x, y = pp$y,
                             marks = pp$marks,
                             window = cur.sub.tile.window)

          cur.sub.tile = as.ppp(cur.sub.tile)

          sub.tile.point.counts[ec] = cur.sub.tile$n
        }

        max.sub.tile = which.max(sub.tile.point.counts)

        cur.sub.tile.window = owin(poly = list(x = cur.win$bdry[[max.sub.tile]]$x,
                                               y = cur.win$bdry[[max.sub.tile]]$y))

        final.poly.list = list(list(x = cur.win$bdry[[max.sub.tile]]$x,
                                    y = cur.win$bdry[[max.sub.tile]]$y))

        # Now check which "holes" are in the window itself

        hole.indices = which(sub.tile.point.counts == -1)

        for(hi in hole.indices){
          if(sum(!inside.owin(w = cur.sub.tile.window,
                              x = cur.win$bdry[[ec]]$x,
                              y = cur.win$bdry[[ec]]$y)) == 0){
            final.poly.list[[length(final.poly.list) + 1]] = list(
              x = cur.win$bdry[[hi]]$x,
              y = cur.win$bdry[[hi]]$y
            )
          }
        }

        cur.win = owin(poly = final.poly.list)
      }


      cur.pointset = ppp(x = pp$x, y = pp$y,
                         marks = pp$marks,
                         window = cur.win)

      cur.pointset = as.ppp(cur.pointset)

      pointset.list[[i]] = cur.pointset
    }

  }

  if(progress) print("Done")
  ret.obj = {}
  ret.obj$tiles = compact(pointset.list)
  ret.obj$pp = final.pp
  ret.obj$pixels = test.v
  ret.obj$window = dirch.tess$window
  return(ret.obj)
}
