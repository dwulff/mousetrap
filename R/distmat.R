#' Compute distance matrix.
#'
#' Computes the point- or vector-wise dissimilarity between
#'   each pair of trajectories.
#'
#' @inheritParams mt_time_normalize
#' @param dimensions a character vector specifying which trajectory variables
#'   should be used. Can be of length 2 or 3 for 2-dimensional or 3-dimensional
#'   trajectories respectively.
#' @param pointwise boolean specifying the way dissimilarity between the 
#'   trajectories is measured. If \code{TRUE} (the default), \code{mt_distmat}
#'   measures the average dissimilarity and then sums the results. If
#'   \code{FALSE}, \code{mt_distmat}  measures dissimilarity once (by treating
#'   the various points as independent dimensions).
#' @param minkowski_p an integer specifying the distance metric. 
#'   \code{minkowski_p = 1} computes the city-block distance, \code{minkowski_p 
#'   = 2} (the default) computes the Euclidian distance, \code{minkowski_p = 3} 
#'   the cubic distance, etc.
#'
#' @return A mousetrap data object (see \link{mt_example}) with an additional 
#'   object added (by default called \code{distmat}) containing the distance
#'   matrix. If a trajectory array was provided directly as \code{data}, only
#'   the distance matrix will be returned.
#'
#' @examples
#' # Spatialize trajectories
#' mt_example <- mt_spatialize(mt_example)
#'  
#' # Compute distance matrix
#' mt_example <- mt_distmat(mt_example, use="sp_trajectories")
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   Jonas M. B. Haslbeck  <jonas.haslbeck@gmail.com>
#'
#' @export
mt_distmat = function(data,
                      use = 'sp_trajectories',
                      save_as = 'distmat',
                      dimensions = c('xpos','ypos'),
                      pointwise = TRUE,
                      minkowski_p = 2
                      ){
  
  # Extract data
  trajectories <- extract_data(data,use) 
  
  # Tests
  if(!length(dimensions) %in% c(2,3)) stop('Dimensions must be of length 2 or 3')
  if(!all(dimensions %in% dimnames(trajectories)[[2]])) stop('Not all dimensions exist')
  
  # Ensure that there are no NAs
  if(any(is.na(trajectories[,dimensions,]))) {
    stop("Missing values in trajectories not allowed for mt_distmat ",
         "as all trajectories must have the same number of observations.")
  }

  # Get distances
  if(length(dimensions) == 2){
    if(pointwise == TRUE){
        dmat <- distMat(trajectories[,dimensions[1],],
                       trajectories[,dimensions[2],],
                       power = minkowski_p)
      } else {
        dmat <- distMatV(trajectories[,dimensions[1],],
                        trajectories[,dimensions[2],],
                        power = minkowski_p)
      }
    } else {
      if(pointwise == TRUE){
          dmat <- distMat3d(trajectories[,dimensions[1],],
                           trajectories[,dimensions[2],],
                           trajectories[,dimensions[3],],
                           power = minkowski_p)
        } else {
          dmat <- distMat3dV(data[[use]][,dimensions[1],],
                            data[[use]][,dimensions[2],],
                            data[[use]][,dimensions[3],],
                            minkowski_p = power)
      }
    }
  
  # Set row and colnames
  rownames(dmat) <- dimnames(trajectories)[[1]]
  colnames(dmat) <- dimnames(trajectories)[[1]]

  # Save and return data
  if(is_mousetrap_data(data)){
    data[[save_as]] = dmat
    } else {
    data = dmat
    }
  return(data)
}


