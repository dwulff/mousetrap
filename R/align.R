#' Align trajectories
#' 
#' \code{trajectory_align} aligns trajectories to common start point, end point, 
#'   and coordinate system.
#' 
#' @param data a mousetrap data object created using one of the mt_import 
#'   functions (see \link{mt_example} for details). Alternatively, a trajectory 
#'   array can be provided directly (in this case \code{use} will be ignored).
#' @param use a character string specifying which trajectory data should be 
#'   used.
#' @param dimensions a character string specifying which trajectory variables 
#'   should be used. Can be of length 2 or 3 for 2-dimensional or 3-dimensional
#'   alignment respectively.
#' @param save_as a character string specifying where the resulting trajectory 
#'   data should be stored.
#' @param start boolean specifying wheter the trajectorie's start points should 
#'   be aligned. 
#' @param end boolean specifying wheter the trajectorie's end points should 
#'   be aligned. 
#' @param end either a numeric vector specying the xstart, ystart, xend, yend
#'   coordinates of the trajectory start and end points, or 'mt' to set the 
#'   coordinates to c(0,0,-1,1.5) or 'norm' to set coordinates to c(0,0,1,1). 
#'   In the 3-dimensional case use coordinates is a vector of size 6.
#'   
#' @return A mousetrap data object (see \link{mt_example}) with an additional
#'   array (by default called \code{aligned_trajectories}) containing the 
#'   aligned trajectories. If a trajectory array was provided directly as
#'   \code{data}, only the aligned trajectories will be returned.
#' 
#' @examples
#' mt_example <- trajectory_align(mt_example,
#'   save_as="aligned_trajectories",  dimensions = c('xpos','ypos'),
#'   start = TRUE,end = TRUE,coordinates = 'mt')
#' 
#' @author Dirk U. Wulff
#' 
#' @export
   
mt_align = function(data,
                    use = 'trajectories',
                    dimensions = c('xpos','ypos'),
                    save_as = 'al_trajectories',
                    start = TRUE, 
                    end = TRUE,
                    coordinates = 'norm'
                    ){
  
  # ---- tests
  if(!length(dimensions) %in% c(2,3)) stop('Dimensions must of length 2 or 3') 
  if(is.character(coordinates) & 
     !any(coordinates == 'mt') & 
     !any(coordinates == 'norm') & 
     !is.numeric(coordinates)){
    stop("Dimensions must be numeric or 'mt' or 'norm'")
    }
  if(any(coordinates == 'mt')){
    if(length(dimensions) == 2){
      coordinates = c(0,0,-1,1.5)
      } else {
      coordinates = c(0,0,0,-1,1.5,-1)
      }
    }
  if(any(coordinates == 'norm')){
    if(length(dimensions) == 2){
      coordinates = c(0,0,1,1)
      }else{
      coordinates = c(0,0,0,1,1,1)
      }
    }
  if(is.list(data)){
    if(!use %in% names(data)) stop(paste0(use,' does not exist!'))
    trajectories = data[[use]][,dimensions,]
    } else {
    if(is.array(data)){
      trajectories = data[,dimensions,]
      }  
    }
  if(!all(dimensions %in% dimnames(trajectories)[[2]])) stop('Not all dimensions exist')


  # ---- flip if necessary
  for(i in dimensions){
    starts = trajectories[,i,1]
    ends   = apply(trajectories[,i,],1,function(x) x[max(which(!is.na(x)))])
    test   = ends <= starts
    if(!mean(test) %in% c(0,1)){
      trajectories[test,i,] =  ((trajectories[test,i,] - starts[test])*-1)+starts[test]
      message(paste('Flipped trajectories to one side with regard to',i))
      }
    }

    
  # ---- replace NAs
  reset_NAs = FALSE
  if(any(is.na(trajectories))){ 
    reset_NAs = TRUE
    trajectories[is.na(trajectories)] =  -3.141592653589793
    }
  
  # ---- aligning data   
  if(length(dimensions) == 2){
    if(length(coordinates) != 4) stop('Coordinates must be a numeric vector of length 4')
    al_trajectories = trajAlign(trajectories[,dimensions[1],],
                                trajectories[,dimensions[2],],
                                start = start,end = end,coordinates = coordinates)
                                }
  if(length(dimensions) == 3){
    if(length(coordinates) != 6) stop('Coordinates must be a numeric vector of length 6')
    al_trajectories = trajAlign3d(trajectories[,dimensions[1],],
                                  trajectories[,dimensions[2],],
                                  trajectories[,dimensions[3],],
                                  start = start,end = end,coordinates = coordinates)
                                  }
  
  # ---- reset NAs
  if(reset_NAs){
  for(i in 1:length(dimensions)){
    al_trajectories[[i]][al_trajectories[[i]] == -3.141592653589793] =  NA
    }
  }
  

  # ---- add aligned data
  if(is.list(data)){
    if(use == save_as){
      for(i in 1:length(dimensions)){
        data[[use]][,dimensions[i],] = al_trajectories[[i]]
        }
      } else {
        data[[save_as]] = data[[use]]
        for(i in 1:length(dimensions)){
          data[[save_as]][,dimensions[i],] = al_trajectories[[i]]
          }
        }
    } else {
      for(i in 1:length(dimensions)){
        data[,dimensions[i],] = al_trajectories[[i]]
        }
      }
  

  return(data)
  }





