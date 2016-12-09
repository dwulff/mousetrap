#' Clean trajectories
#' 
#' Identify (and remove) unusal trajectories based on (a) their most extreme 
#' points and (b) their similarity to other trajectories. Specifically, the 
#' function evaluates whether the trajectories most extreme points exceed 
#' specified limits on the specified dimensions, as well as whether a 
#' trajectory's average similarity to all other trajectories lies below a
#' given similarity threshold.
#' 
#' @param data a mousetrap data object created using one of the mt_import 
#'   functions (see \link{mt_example} for details). Alternatively, a trajectory 
#'   array can be provided directly (in this case \code{use} will be ignored).
#' @param use a character string specifying which trajectory object should be
#'   used. 
#' @param dimensions a character string specifying which trajectory variables 
#'   should be used. Can be of length 2 or 3 for 2-dimensional or 3-dimensional
#'   trajectories.
#' @param save_as a character string specifying how the clean variables should be
#'   stored. If save_as is an appropriately-sized data frame existing in data, 
#'   the variables are added, otherwise a new data frame will be created.
#' @param box_limits a numeric vector matching the length of dimensions that 
#'   specifying the maximum allowed distance from the start point (measured
#'   in terms of the distance between the trajectorie's start and end point).
#' @param sim_limit a numeric specifying the maximum allowed z-value
#'   of a trajectory relative to the distribution of average trajectory 
#'   similarities.
#' @param aligned_sim boolean specifying whether the trajectories specified
#'   in \code{use} should be aligned (see \link{trajectory_align} for details) 
#'   or not before assessing similarity between trajectorie.
#' @param rescaled_sim boolean specifying whether the trajectories specified
#'   in \code{use} should be spatially normalized (see \link{spatial_rescale} 
#'   for details) or not before assessing similarity between trajectories. If 
#'   FALSE \code{n_rescale} is ignored.
#' @param minkowski_p numeric specifying the distance function used to
#'   assess similarity between trajectories. \code{minkowski_p} of 1 and 2 
#'   translate into city-block and euclidian distances.
#' @param clean character specifying the variable used to clean the data. Can
#'   be "box", "sim", "union, and "intersection", which will remove trials 
#'   identified by the box_limit assessment, the similarity_limit, one of them,
#'   or both, respectively.
#'   
#' @return A mousetrap data object (see \link{mt_example}) with added columns in 
#'   the data frame specified by \code{save_as} or a new data frame of that name. 
#' 
#' @examples
#' mt_example <- clean_trajectories(mt_example,
#'   save_as="data", dimensions = c('xpos','ypos'),
#'   box_limits = c(2,2), sim_limit  = 2, rescaled_sim = TRUE,
#'   n_rescale = 20, minkowski_p = 2, clean = FALSE)
#'   
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   
#' @export

mt_clean = function(data,
                            use = 'trajectories',
                            dimensions = c('xpos','ypos'),
                            save_as = 'measures',
                            method = 'intersect',
                            box_limits = c(-2,2,-.2,2),
                            sim_limit  = 2,
                            minkowski_p = 2,
                            clean = TRUE
                            ){
  
  # ---- tests
  if(is.list(data)){
    if(!use %in% names(data)) stop(paste0(use, ' does not exist!'))
    trajectories = data[[use]]
    } else {
    if(is.array(data)){
      trajectories = data
      } else {
      stop('Argument data must be either list or array!')  
      }  
    }
  if(!all(dimensions %in% dimnames(trajectories)[[2]])) stop('Not all dimensions exist')
  
  
  # ---- initialize variables
  box_clean = NULL
  sim_clean = NULL
  max_dist  = NULL
  mean_sim  = NULL
  
  
  # ---- box clean
  # Remove trajectories that have points beyond box.xlim and box.ylim
  if(method %in% c('box','intersect','union')){
    box_clean = rep(FALSE,dim(trajectories)[1])
    for(i in 1:length(dimensions)){
      max_pos    = apply(trajectories[,dimensions[i],],1,function(x) c(min(x,na.rm=TRUE),max(abs(x),na.rm=TRUE)))
      box_clean  = box_clean | max_pos[1,] < box_limits[i*2-1] | max_pos[2,] > box_limits[i*2]
      }
    }
  
  # ---- sim clean
  # determine average similarity of trajectory to all other trajectories
  if(method %in% c('sim','intersect','union')){
    dist_mat  = mt_distmat(trajectories, dimensions = dimensions, use = use, power = minkowski_p)      
    mean_sim  = scale(-apply(dist_mat,1,mean))[,1]
    sim_clean = mean_sim < -sim_limit
    }


  # ---- remove dirty trajectories or save indices
  if(clean == TRUE){
    if(method == 'box')            cleaner = box_clean
    if(method == 'sim')            cleaner = sim_clean
    if(method == 'union')          cleaner = box_clean | sim_clean
    if(method == 'intersection')   cleaner = sim_clean & box_clean
    for(nam in names(data)){
      if(length(dim(data[[nam]])) == 2) data[[nam]] = data[[nam]][!cleaner,]  
      if(length(dim(data[[nam]])) == 3) data[[nam]] = data[[nam]][!cleaner,,]       
      }
    cat(paste0(sum(cleaner), ' trials removed!'),'\n')
    } else {
    if(save_as %in% names (data)){
      data[[save_as]]$box_clean = box_clean
      data[[save_as]]$sim_clean = sim_clean 
      data[[save_as]]$max_dist  = max_dist
      data[[save_as]]$mean_sim  = mean_sim
      } else {
      data[[save_as]] = data.frame(box_clean,sim_clean,max_dist,mean_sim)
      }
    }
  
  return(data)
  }

