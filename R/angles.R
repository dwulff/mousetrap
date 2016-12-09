#' Calculate movement angles.
#' 
#' Calculate point-based and vertical-based angles for each of the n-2 points in
#' the movement trajectory. Point-based angles are the angle defined by three
#' subsequent points on the trajectory. Vertical-based angles are the angles
#' between two subsequent points and the vertical axis.
#' 
#' By default, angles are reported in radian, the alternative is degree. For the
#' first point in a trajectory, the angle values are always not defined
#' (\code{NA})
#' 
#' For vertical-based angles (\code{angle_v}), positive values indicate a 
#' movement to the left of the vertical, negative values to the right of the 
#' vertical. If there was no movement across two consecutive points, 
#' \code{angle_v} is not defined and, by default, \code{NA} is returned. If 
#' \code{na_replace} is \code{TRUE}, the next existing angle value is reported
#' instead.
#'   
#' For point-based angles (\code{angle_p}), angles indicate changes of movement 
#' within three consecutive time steps. Always the smaller angle is reported. Pi
#' (= 3.14) (radians) / 180 (degree) indicate a constant movement direction, 0 
#' (radians) / 0 (degree) a complete reversal. If there was no movement across
#' two consecutive points, \code{angle_p} is not defined and, by default,
#' \code{NA} is returned. If \code{na_replace} is \code{TRUE}, the next existing
#' angle value is reported instead.
#'   
#'   
#' @inheritParams mt_time_normalize
#' @param dimensions a character string specifying which trajectory variables 
#'   should be used. Must be of length 2.
#' @param na_replace logical specifying whether \code{NA}s in the angle values 
#'   should be replaced using the next existing angle value (see Details).
#'   Defaults to \code{FALSE}.
#' @param unit character specifying the unit for the angles. Default is 
#'   "radian", alternative is "degree".
# 
#' @author Dirk U. Wulff
#'  
#' @return A mousetrap data object (see \link{mt_example}) with point-based and 
#'   vertical-based angles added as additional columns to the trajectory array
#'   (called \code{angle_p} and \code{angle_v}). If a trajectory array was
#'   provided directly as \code{data}, only the trajectory array will be
#'   returned.
#'
#' @examples
#' # Calculate movement angles
#' mt_example <- mt_angles(mt_example)
#'
#' # Calculate movement angles (in degree)
#' # and replace NAs with next existing value
#' mt_example <- mt_angles(mt_example,
#'   unit="degree", na_replace=TRUE)
#' 
#' @export
mt_angles = function(data,
                     use = 'trajectories',
                     dimensions = c('xpos','ypos'),
                     save_as = use,
                     na_replace = FALSE,
                     unit = "radian",
                     verbose = FALSE
                     ){
  # Extract data
  trajectories <- extract_data(data,use) 
  
  # Tests
  n_dim <- length(dimensions)
  if(!n_dim %in% c(2,3)) stop('Dimensions must of length 2') 
  if(!unit %in% c("radian","degree")) stop('Unit must either be "radian" or "degree"')
  if(!all(dimensions %in% dimnames(trajectories)[[2]])) stop(paste0('Not all dimensions exist!'))
  
  # Get angles
  anglesP <- getAnglesP(trajectories[,dimensions[1],],
                       trajectories[,dimensions[2],])
  anglesV <- getAnglesV(trajectories[,dimensions[1],],
                       trajectories[,dimensions[2],])
  

  # Replace NAs (where angle could not be calculated, see documentation)
  if(na_replace == TRUE){
    cleanAngles(anglesP) 
    cleanAngles(anglesV)
  } else {
    anglesP[anglesP == -100] <- NA
    anglesV[anglesV == -100] <- NA
  }
  
  # Set NAs for cases where no dimension values where found in the data
  anglesP[anglesP == -90] <- NA
  anglesV[anglesV == -90] <- NA

  
  # Convert to degree (optional)
  if (unit == "degree") {
    anglesP <- (anglesP * 180) / (pi)
    anglesV <- (anglesV * 180) / (pi)
  } 

  # Add angles to trajectories
  trajectories <- mt_add_variables(trajectories,
                           variables=list('angle_p'=anglesP,
                                          'angle_v'=anglesV))

  return(create_results(data=data, results=trajectories, use=use, save_as=save_as))
}



