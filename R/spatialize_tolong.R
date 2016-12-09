#' Spatialize trajectories to long
#'
#' \code{spatialize} re-represents the trajectory spatially so that
#'   adjacent points on the trajectory become equidistant to each other.
#'
#' @param data a mousetrap data object created using one of the mt_import
#'   functions (see \link{mt_example} for details). Alternatively, a trajectory
#'   array can be provided directly (in this case \code{use} will be ignored).
#' @param use a character string specifying which trajectory data should be
#'   used.
#' @param dimensions a character string specifying which trajectory variables
#'   should be used. Can be of length 2 or 3 for 2-dimensional or 3-dimensional
#' @param save_as a character string specifying where the resulting trajectory
#'   data should be stored.
#' @param n_points an integer vector specifying the number of points used to
#'   represent the spatially rescaled trajectories. Length must be 1 or the number
#'   number of trajectories in the data. The values in \code{n_points} may vary.
#'
#' @return A mousetrap data object (see \link{mt_example}) with an additional
#'   array (by default called \code{spatialized}) containing the
#'   aligned trajectories. If a trajectory array was provided directly as
#'   \code{data}, only the aligned trajectories will be returned.
#'
#' @examples
#' mt_example <- mt_spatialize(data=mt_example,
#'   dimensions = c('xpos','ypos'), save_as="rescaled_trajectories",
#'   n_points = 20, format = 'wide')
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   Jonas M. B. Haslbeck  <jonas.haslbeck@gmail.com>



mt_spatialize_tolong = function(data,
                                use = 'trajectories',
                                dimensions = c('xpos','ypos'),
                                save_as = 'spatialized_trajectories',
                                n_points = 20
                                ){
  
  # ---- tests
  if(!length(dimensions) %in% c(2,3)){
    stop('Dimensions must of length 2 or 3')
  }
  
  trajectories <- mousetrap:::extract_data(data,use) 
  
  # if(is.list(data)){
  #   if(!use %in% names(data)){
  #     stop(paste0(use,' does not exist!'))
  #     }
  #   trajectories = data[[use]]
  #   } else {
  #   if(is.array(data)){
  #     trajectories = data
  #     }
  #   }
  
  if(!all(dimensions %in% dimnames(trajectories)[[2]])) stop('Not all dimensions exist')
  

  if(length(dimensions) == 2){
    if(nrow(trajectories) == 1){
      dim_1 = matrix(trajectories[,dimensions[1],],nrow=1)
      dim_2 = matrix(trajectories[,dimensions[2],],nrow=1)
      } else {
      dim_1 = trajectories[,dimensions[1],]
      dim_2 = trajectories[,dimensions[2],]
      }
    spatialized_trajectories = spatializeArrayToLong(dim_1,dim_2,n_points)
    }
  if(length(dimensions) == 3){
    if(nrow(trajectories) == 1){
      dim_1 = matrix(trajectories[,dimensions[1],],nrow=1)
      dim_2 = matrix(trajectories[,dimensions[2],],nrow=1)
      dim_3 = matrix(trajectories[,dimensions[3],],nrow=1)
      } else {
      dim_1 = trajectories[,dimensions[1],]
      dim_2 = trajectories[,dimensions[2],]
      dim_3 = trajectories[,dimensions[3],]
      }
    spatialized_trajectories = spatializeArrayToLong3d(dim_1,dim_2,dim_3,n_points)
    }
  
  
  # ---- save data
  return(spatialized_trajectories)
  }
  