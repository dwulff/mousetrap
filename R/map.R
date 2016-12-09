#' Map trajectories.
#'
#' \code{mt_map} maps trajectories onto a predefined set of prototype
#'   trajectories.
#'
#' @inheritParams mt_distmat
#' @param proto_list list specifying the list of prototypes the trajectories are
#'   mapped to. Each prototype must be a n x 2 matrix, where the columns correspond
#'   to the to be mapped trajectory dimensions.
#' @param point_wise boolean specifying the way dissimilarity between the
#'   trajectories is measured. \code{point-wise = TRUE} measures the average
#'   dissimilarity and then sums the results. \code{point-wise = FALSE} measures
#'   dissimilarity once (by treating the various points as independent dimensions).
#' @param minkowski_p an integer specifying the distance metric. \code{minkowski_p}
#'   = 1 computes the city-block distance, \code{minkowski_p} = 2 computes the
#'   Euclidian distance, and so on.
#'
#' @param A mousetrap data object (see \link{mt_example}) with an additional
#'   object (by default called \code{map}) containing a vector or list indiciating
#'   set of mapped prototypes. If a trajectory array was provided directly as
#'   \code{data}, only the vector or list containing the results will be returned.
#'
#' @examples
#' proto_list <- readRDS(prototypes.RDS)
#' mt_example <- mt_map(data=mt_example,
#'   proto_list = proto_list,
#'   dimensions = c('xpos','ypos'),
#'   save_as="map")
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   Jonas M. B. Haslbeck  <jonas.haslbeck@gmail.com>
#'
#' @export



mt_map = function(
  data,
  use = 'trajectories',
  dimensions = c('xpos','ypos'),
  save_as = 'map',

  # prototype arguments
  proto_list = 'default', # list object containing prototypes in a n x 2 -matrix

  # distance arguments
  pointwise = TRUE,
  minkowski_p = 2
  )

{

  # ---- tests
  if(!length(dimensions) %in% c(2,3)){
    stop('Dimensions must of length 2 or 3!')
    }
  if(is.list(data)){
    trajectories = data[[use]]
    } else {
    if(is.array(data)){
      trajectories = data
      } else {
      stop('Data must be array or list!')
      }
    }
  if(!all(dimensions %in% dimnames(trajectories)[[2]])) stop('Not all dimensions exist')
  if(is.character(proto_list)){
    if(proto_list == 'default'){
      proto_list = load('~/Dropbox (2.0)/Work/Software/mousetrap/data/prototypes.rda')
      }
    }

  # rescale prototypes
  n_points = dim(trajectories)[3]
  n_proto  = length(proto_list)
  joint_array = array(dim = c(dim(trajectories)[1] + n_proto,
                            length(dimensions),
                            n_points),
                            dimnames = list(c(paste0('proto_',1:n_proto),dimnames(trajectories)[[1]]),
                                            dimensions,
                                            NULL))

  for(i in 1:n_proto){
    joint_array[i,,] = mt_spatialize(proto_list[[i]],
                                       n_points = n_points,
                                       dimensions = dimnames(proto_list[[i]])[[2]])
    }

  joint_array[(n_proto + 1) : dim(joint_array)[1],dimensions,] = trajectories[,dimensions,]


  # ---- compute distance & closest prototype
  distm =  mt_distmat(joint_array,
                      dimensions = dimensions,
                      pointwise = pointwise,
                      minkowski_p = minkowski_p)
  dists = distm[1:n_proto,-c(1:n_proto)]
  prototypes = apply(dists,2,function(x) which(x == min(x)))


  # ---- save data
  if(is.list(data)){
    if(save_as %in% names(data)){
      data[[save_as]]$prototypes = prototypes
      } else {
      result = data.frame(mt_id = dimnames(trajectories)[[1]],prototyping = prototypes)
      data[[save_as]] = result
      }
    } else {
    result = data.frame(mt_id = dimnames(trajectories)[[1]],prototyping = prototypes)
    data = result
    }

  return(data)
  }
