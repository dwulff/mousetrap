#' Cluster trajectories.
#'
#' \code{mt_cluster} performs trajectory clustering.
#' 
#' ToDo: Provide details on different clustering methods and the resulting data.
#'
#' @inheritParams mt_time_normalize
#' @param method character string specifiying the type of clustering procedure. 
#'   Either \link[fastcluster]{hclust} (the default) or \link[stats]{kmeans}.
#' @param pointwise boolean specifying the way dissimilarity between the 
#'   trajectories is measured. If \code{TRUE} (the default), \code{mt_distmat} 
#'   measures the average dissimilarity and then sums the results. If 
#'   \code{FALSE}, \code{mt_distmat}  measures dissimilarity once (by treating 
#'   the various points as independent dimensions). Only relevant if
#'   \code{method} is "hclust".
#' @param minkowski_p an integer specifying the distance metric. 
#'   \code{minkowski_p = 1} computes the city-block distance, \code{minkowski_p 
#'   = 2} (the default) computes the Euclidian distance, \code{minkowski_p = 3} 
#'   the cubic distance, etc. Only relevant if \code{method} is "hclust".
#' @param hclust_method character string specifying the linkage criterion used. 
#'   Passed on to the \code{method} argument of \link[stats]{hclust}. Only 
#'   relevant if \code{method} is "hclust".
#' @param kmeans_nstart integer specifying the number of reruns of the kmeans 
#'   procedure. Large numbers minimize the risk of finding local minima.  Passed
#'   on to the \code{nstart} argument of \link[stats]{kmeans}. Only relevant if 
#'   \code{method} is "kmeans".
#'
#' @return A mousetrap data object (see \link{mt_example}) with an additional 
#'   \link{data.frame} added to it (by default called \code{clustering}) that
#'   contains the cluster assignments. If a trajectory array was provided
#'   directly as \code{data}, only the clustering data.frame will be returned.
#'
#' @examples
#' # Spatialize trajectories
#' mt_example <- mt_spatialize(mt_example)
#'  
#' # Cluster trajectories
#' mt_example <- mt_cluster(mt_example, use="sp_trajectories")
#' 
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   Jonas M. B. Haslbeck  <jonas.haslbeck@gmail.com>
#'
#' @export

mt_cluster = function(data,
                      use = 'sp_trajectories',
                      save_as = 'clustering',
                      dimensions = c('xpos','ypos'),
                      n_cluster = 5, # = k
                      method    = 'hclust',
                      
                      # distance arguments
                      pointwise   = TRUE,
                      minkowski_p = 2,

                      # cluster arguments
                      hclust_method   = 'ward.D',
                      kmeans_nstart   = 10
                      ){

  # Extract data
  trajectories <- extract_data(data,use) 
  
  # Tests
  if(!length(dimensions) %in% 2) stop('dimensions must of length 2!')
  if(!all(c(dimensions) %in% dimnames(trajectories)[[2]])) stop('Not all dimensions exist')
  if(!method %in% c("hclust","kmeans")) stop('method must either be "hclust" or "kmeans"')
  
  # Ensure that there are no NAs
  if(any(is.na(trajectories[,dimensions,]))) {
    stop("Missing values in trajectories not allowed for mt_distmat ",
         "as all trajectories must have the same number of observations.")
  }

  # Cluster trajectories
  
  # ... using hclust method
  if(method == 'hclust') {
    distm =  mt_distmat(
      trajectories,
      dimensions = dimensions,
      pointwise = pointwise,
      minkowski_p = minkowski_p)
    distm = stats::as.dist(distm)

    # clustering
    cl_obj = fastcluster::hclust(distm, method = hclust_method)
    cl_ass = stats::cutree(cl_obj, n_cluster)
    
    # ... using kmeans method
    } else {

    # rearrange data structura for clustering input
    rearranged_trajectories = cbind(trajectories[,dimensions[1],],trajectories[,dimensions[2],])

    # k-means
    cl_obj <- stats::kmeans(
      rearranged_trajectories, centers = n_cluster, nstart = kmeans_nstart)
    cl_ass <- cl_obj$cluster

    }


  # Save data
  return(create_results(
    data=data, results=data.frame(cluster=cl_ass), use=use, save_as=save_as,
    ids=rownames(trajectories), overwrite=TRUE))
  

}
