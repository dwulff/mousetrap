#' Detect change of mind trajectories
#' 
#' \code{mt_com} detects change-of-mind (com) trials based on the 
#'   trajectory's degree of deflection and the development of trajectory
#'   angles.
#' 
#' @param data a mousetrap data object created using one of the mt_import 
#'   functions (see \link{mt_example} for details). Alternatively, a trajectory 
#'   array can be provided directly (in this case \code{use} will be ignored).
#' @param use a character string specifying which trajectory data should be 
#'   used.
#' @param dimensions a character string specifying which trajectory variables 
#'   should be used. Can be of length 2 or 3 for 2-dimensional or 3-dimensional
#'   assessment respectively.
#' @param save_as a character string specifying where the resulting trajectory 
#'   label data should be stored.
#' @param deflection_threshold an integer vector specifying the thresholds for the
#'   trajectory's deflection. Deflection is measured in terms of movement in the 
#'   direction opposite of the chosen option concerning the first, second, ..., 
#'   nth dimension. Deflection is expressed in terms of percentiles measured from
#'   the midpoint between the options.
#' @param angle_threshold an integer specifying the threshold for the trajectory's
#'   angles. Angles are determined for the 2nd, 3rd, ..., nth-1 point by measuring,
#'   for instance, for the 2nd point the angle between the line defined by the 1st
#'   and 2nd point and the line defined by the 2nd and 3rd point.      
#'
#' @return A mousetrap data object (see \link{mt_example}) with added columns in 
#'   the data frame specified by \code{save_as} or a new data frame of that name. 
#' 
#' @examples
#' mt_example <- detect_CoM(data=mt_example,
#'   dimensions = c('xpos','ypos'), save_as="data",
#'   deflection_thresholds  = 30, angle_threshold = 50)
#'
#' @author Dirk U. Wulff
#'         
#' @export

mt_com = function(data,
                      use        = 'trajectories',
                      dimensions = c('xpos','ypos'),
                      save_as    = 'data',
                      deflection_thresholds  = 30,
                      angle_threshold = 50
                      ) {

  # ---- Tests
  if(!use %in% names(data)) stop(paste0(use,' does not exist!'))
  if(!all(dimensions %in% dimnames(data[[use]])[[2]])) stop(paste0('Not all dimensions exist!'))
  
    
  # ---- norm data
  data_tmp = trajectory_align(data, 
                              use = use, 
                              dimensions = dimensions, 
                              coordinates = 'norm')
  
  # ---- distance based detection
  deflection_CoM = logical()
  for(i in 1:length(deflection_tresholds)){
    max_deflection = apply(data_tmp[[use]][,dimensions[i],],1,function(x) min(x))
    deflection_CoM = deflection_CoM  | max_deflection > -threshold_dist/100
    }
  
  # ---- angle based detection
  if(!all(c('anglesP') %in% dimnames(data_tmp[[use]]))){
    data_tmp = add_angles(data_tmp,
                          dimensions = dimensions,
                          use = use,
                          overwrite = FALSE)
                          }
  min_angles = apply(data_tmp[[use]][,'anglesP',],1,function(x) min(x*(180/pi),na.rm=T))
  angle_CoM  = min_angles < threshold_angle
  
  
  # ---- save data
  if(save_as %in% names (data)){
    data[[save_as]]$dist_CoM = dist_CoM
    data[[save_as]]$angle_CoM = angle_CoM 
    data[[save_as]]$max_deflection  = max_deflection
    data[[save_as]]$min_angle  = min_angle
    } else {
    data[[save_as]] = data.frame(dist_CoM,angle_CoM,max_deflection,min_angle)
    }

  return(data)
  }


