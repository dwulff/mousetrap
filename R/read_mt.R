#' Read MouseTracker files (.mt)
#' 
#' \code{read_mt} reads raw data from MouseTracker files (.mt)
#' 
#' @param File a character string specifying the filename of the .mt file. 
#' @param variables either 'all' or a characteric vector specifying the to be 
#'   extracted variables Defaults to 'all' in which case all existing variables
#'   will be extracted.
#' @param add_id a boolean specifying whether the filename and a trialnumber should
#'   be added to the data frame.
#'   
#' @return A data frame with one column per trial. Variables are ordered according
#'   to columns, x-coordinates, y-coordinates, and timestamps.
#'   array (by default called \code{aligned_trajectories}) containing the 
#' 
#' @examples
#' mt_example <- read_mt(File,columns = 'all', add_id = FALSE)
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'   
#' @export

read_mt = function(File, 
                   columns = 'all',
                   add_id  = FALSE
                   ){
  
  # ---- tests
  if(remove != 'auto' & !is.numeric(remove)) stop('remove must of value "auto" or a numeric vector')
  
  # ---- definitions
  block_names = c('RAW TRACKS (X coordinates)',
                  'RAW TRACKS (Y coordinates)',
                  'RAW TRACKS (time)')  
  variable_labels = c('X','Y','T')
  
  
  # ---- read
  # read file
  con  = file(File,'rb')
  lins = readLines(con)
  close(con)
  
  # Loop over blocks
  for(i in 1:length(block_names)){
    block      = block_names[i] 
    start      = which(lins == block)
    end        = start + which(lins[start:length(lins)] == '')[1] - 2
    nams       = strsplit(lins[start+1],',')[[1]]
    if(length(columns) == 1){
      if(columns == 'all') vars = nams
      } else {
        vars = columns
      }
    n_nams = length(nams)
    trials = (start + 2) : end
    

    # ---- process trials    
    
    ind = 0
    positions = list()
    if(i == 1) front = matrix(NA, nrow = length(trials), ncol = length(vars) + ifelse(add_filename,2,0))
    for(trial in trials){
      ind = ind + 1
      lin = strsplit(enc2utf8(lins[trial]),',')[[1]]      
      positions[[ind]] = as.numeric(lin[-c(1:n_nams)])
      if(i == 1){
        if(add_id == TRUE){
          front[ind,] = c(File,ind,lin[which(nams %in% vars)])
          } else {
          front[ind,] = lin[which(nams %in% vars)]       
          }
        }
      }
    if(i == 1){
      if(add_id == TRUE){
        colnames(front) = c('file','trial',vars)
        } else {
        colnames(front) = vars
        }
      }
    
    # ---- store trials
    
    # store in matrix
    valid_trials = which(sapply(positions,length) > 0) ; ind = 0
    pos_matrix = matrix(NA, nrow = length(valid_trials), ncol = max(sapply(positions,length)))  
    colnames(pos_matrix) = paste(variable_labels[i], 1:ncol(pos_matrix), sep='_')
    for(j in valid_trials){
      ind = ind + 1
      pos = positions[[j]]
      pos_matrix[ind,1:length(pos)] = pos
      }
    
    # combine matrices
    if(i == 1){
      data = cbind(front[valid_trials,],pos_matrix)
      } else {
      data = cbind(data,pos_matrix)
      }
    }
  
  data = data.frame(data)
  
  return(data)
  }



