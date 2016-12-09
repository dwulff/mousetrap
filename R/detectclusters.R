#' Detect number of clusters
#'
#' \code{mt_detectclusters} estimates the number of clusters using various
#'   methods.
#'
#' @inheritParams mt_cluster
#' @param compute character vector specifying the to be computed measures. Can
#'   be any subset of \code{c('stability','gap','jump','slope')}.
#' @param method character string specifiying the type of clustering procedure
#'   for the stability-based method. Either \code{hclust} or \code{kmeans}.
#' @param kseq ToDo
#' @param n_bootstrap ToDo
#' @param model_based ToDo
#' @param n_gap ToDo
#'
#' @return ToDo.
#'
#' @examples
#' # Spatialize trajectories
#' mt_example <- mt_spatialize(mt_example)
#'  
#' # Detect clusters
#' results <- mt_detectclusters(mt_example, use="sp_trajectories")
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   Jonas M. B. Haslbeck  <jonas.haslbeck@gmail.com>
#'
#' @export
mt_detectclusters = function(
  data,
  use = 'sp_trajetcories',
  dimensions = c('xpos','ypos'),
  
  # k-selection type
  compute  = c('stability','gap','jump','slope'),
  kseq    = 2:15, # considered k-sequence

  # distance arguments
  pointwise = TRUE, #
  minkowski_p = 2,

  # stability arguments
  method  = 'hclust', # or 'kmeans' #todo hierarchical
  hclust_method      = 'ward.D', # hclust_method
  kmeans_nstart         = 10, # number of reinitializations for k-means algorithm
  n_bootstrap          = 10, # bootstrap samples
  model_based  = FALSE, # stab_predict

  # arguments distance-based k-selection methods
  n_gap = 10, # simulated datasets for Gap Statistic

  verbose=TRUE
){
  
  # Extract data
  trajectories <- extract_data(data,use)

  # Tests
  if(model_based == TRUE & method == 'hierarchical'){
    stop('Model-based instability methods are only available for k-means clustering.')
  }
  if (method=="hclust") {
    method <- "hierarchical"
  } else if (method != "kmeans") {
    stop('method must either be "hclust" or "kmeans"')
  }

  # Transform data structure for clustering input
  rearranged_trajectories = cbind(trajectories[,dimensions[1],],trajectories[,dimensions[2],])
  
  # Define containers
  kopt = list()
  seqs = list()

  
  # Stability-based k-selection methods
  if('stability' %in% compute) {
    
    if(verbose == TRUE) cat('calculating stability-based k-selection methods','\n')
    cStab_obj <- cstab::cStability(data = rearranged_trajectories,
                            kseq = kseq,
                            nB = n_bootstrap,
                            norm = TRUE,
                            predict = model_based,
                            method = method,
                            linkage = hclust_method,
                            kmIter = kmeans_nstart,
                            pbar = FALSE)
    
    kopt[['stab_kopt']] <- unlist(cStab_obj$k_instab_norm)
    seqs[['stab_seq']] <- cStab_obj$instab_path_norm
    
  }
    
  # ---- Distance-based k-selection methods
  if(any(compute %in% c('gap','jump','slope'))) {

    if(verbose == TRUE) cat('calculating distance-based k-selection methods','\n')
    cDist_obj <- cstab::cDistance(data = rearranged_trajectories,
                           kseq = kseq,
                           method = method,
                           linkage = hclust_method,
                           kmIter = kmeans_nstart,
                           gapIter = n_gap)

    if('gap' %in% compute){
      kopt[['gap_kopt']] <- unlist(cDist_obj$k_Gap)
      seqs[['gap_seq']] <- cDist_obj$Gaps
    }
    if('jump' %in% compute){
      kopt[['jump_kopt']] <- unlist(cDist_obj$k_Jump)
      seqs[['jump_seq']] <- cDist_obj$Jumps
    }
    if('slope' %in% compute){
      kopt[['slope_kopt']] <- unlist(cDist_obj$k_Slope)
      seqs[['slope_seq']] <- cDist_obj$Slopes
    }
  }

  #output
  output = list(kopt=kopt,seqs=seqs)
  return(output)

  }





