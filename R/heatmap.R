#' Creates high-resolution heatmap of trajectory data
#'
#' \code{mt_heatmap_raw} creates a high-resolution heatmap of the trajectory
#'   data using gaussian smoothing.

#' @inheritParams mt_time_normalize
#' @param dimensions a character vector specifying the trajectory variables
#'   used to create the heatmap. The first two entries are used as x and
#'   y-coordinates, the third, if provided, will be added as color information.
#' @param bounds numeric vector specifying the corners (xmin, ymin, xmax, ymax) of
#'   the plot region. Note that the trajectories' start and end points are mapped to
#'   (0, 0) and (-1, 1.5), respectively. Bounds should therefore include the rectangle
#'   c(-1, 0, 1, 1.5). Defaults to c(-1.4, -.3, 1.4, 1.9).
#' @param px_size an integer specifying the size of px. A px_size of .01 implies an
#'   x-resolution of (xmax-xmin)/px_size and a y-resolution of (ymax-ymin)/px_size.
#'   Defaults to .001, resulting in an image size of 2800x2200 px if the \code{bounds} are unchanged.
#' @param upscale a numeric value by which the output resolution of the image is increased or
#'   decreased. Only applies if \code{device} is one of \code{c("tiff", "png", "pdf")}.
#' @param upsample a numeric by which the number of points used to represent individual
#'   trajectories are increased or decreased. Values of smaller than one will improve speed
#'   but also introduce a certain level of granularity.
#' @param smooth_radius a numeric specifying the standard deviation of the gaussian smoothing.
#'   If zero, smoothing is omitted.
#' @param low_pass an integer specifying allowed number of counts per pixel. This
#'   arguments limits the maximum pixel color intensity.
#' @param mean_intensity a numeric between 0 and 1 specifying the average color intensity
#'   across the entire image. Defaults to .2.
#' @param norm a logical specifying whether the data should be warped into standard space.
#'   Norm also sets standard values for bound = c(-1.4, -.3, 1.4, 1.9). 
#' @param color a character string or numeric vector of length three specifying
#'   the color used to plot the third dimension.
#' @param color_order a character string, either "increasing" or "decreasing" specifying the
#'   whether large or small values of the third dimension are assigned high color values,
#'   respectively.
#' @param n_trajectories an integer specifying the number of trajectories used to create the
#'   image. If \code{n_trajectories} is smaller than containes in the trajectorie object
#'   specified by \code{use} the \code{n_trajectories} are randomly sampled.
#' @param seed an integer specifying the seed used for the trajectory sampling.
#' @param verbose boolean specifying whether progress updates should printed.
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'
#' @return Nothing, when image is plotted using an external device. Otherwise a matrix
#'    containing the image's pixel information.
#' @export

mt_heatmap_raw = function(
  data,
  use = 'trajectories',
  dimensions = c('xpos', 'ypos'),

  # plot arguments
  bounds    = NULL,
  xres      = 1000,
  upscale   = 1,
  upsample  = 1,
  smooth_radius  = 1.5,
  low_pass  = 200,
  auto_enhance = TRUE,
  mean_intensity = .1,
  norm = TRUE,

  # color arguments
  colors      = c('black', 'blue'),
  n_shades    = c(1000, 10),

  # plot aggregate
  aggregate_lwd = 0,
  aggregate_col = 'black',

  # subsample arguments
  n_trajectories = 10000,
  seed = 1,

  # control
  verbose = TRUE
){

  # Get time
  t = proc.time()[3]

  # Data checks
  if (!length(dimensions) %in% c(2, 3)) {
    stop('Dimensions must of length 2 or 3!')
  }

  if (is.list(data)) {
    trajectories = data[[use]]
  } else {
    if (is.array(data)) {
      trajectories = data
    } else {
      stop('Data must be array or list!')
    }
  }

  if (!all(dimensions %in% dimnames(trajectories)[[2]])) {
    stop('Not all dimensions exist in data')
  }

  # Subsample trajectories -----------------------------------------------------
  # If n_trajectories is smaller than number of trajectories,
  # subsample data to n_trajectories
  if (n_trajectories < dim(trajectories)[1]) {
    if (verbose == TRUE) cat('subset trajectories','\n')
    trajectories = subsample(trajectories, n=n_trajectories, seed=seed)
  }

  # Get aggregate --------------------------------------------------------------
  # Compute aggregate x and y
  aggregate = cbind(
    colMeans(trajectories[,dimensions[1],],na.rm=T),
    colMeans(trajectories[,dimensions[2],],na.rm=T)
  )

  # Determine Dimensions ----------------------------------------------------
  # Rescale trajectories so that there are sufficient points to fill
  # the diagonal, i.e, length of diagonal divided by px

  if (verbose == TRUE) cat('spatializing trajectories','\n')

  if(is.null(bounds) & norm == FALSE){
    range_x = range(trajectories[,dimensions[1],],na.rm=T)
    range_y = range(trajectories[,dimensions[2],],na.rm=T)
    mid_x   = range_x[1] + diff(range_x) / 2
    mid_y   = range_y[1] + diff(range_y) / 2
    range_x = mid_x + (range_x - mid_x) * 1.02
    range_y = mid_y + (range_y - mid_y) * 1.02
    bounds = c(range_x[1],range_y[1],range_x[2],range_y[2])
    }
  if(norm == TRUE){
    trajectories = mt_align(trajectories,coordinates = 'mt')
    bounds = c(-1.4, -.3, 1.4, 1.9)
  }
  
  margin = ceiling(smooth_radius * 2) + ceiling(smooth_radius * 2)
  x_canv = xres - margin - 1
  
  px_size = diff(range_x) / x_canv
  
  # Determine number of resc
  xres   = x_canv
  yres   = ceiling(x_canv * diff(range_y) / diff(range_x))
  n_resc = ceiling(sqrt(xres * xres + yres * yres))
  l_diag = sqrt((bounds[3] - bounds[1])**2 + (bounds[4] - bounds[2])**2)

  
  # Determine Dimensions ----------------------------------------------------
  
  # Determine number of resc
  lengths  = mousetrap:::getLengths(
    trajectories[,dimensions[1],],
    trajectories[,dimensions[2],]
  )
  n_points = round(n_resc * (lengths / l_diag))
  n_points = n_points * upsample

  if (length(dimensions) == 2) {
    spatialized_trajectories = mousetrap:::mt_spatialize_tolong(
      trajectories,
      dimensions = dimensions,
      n_points = n_points
    )
  }

  if (length(dimensions) == 3) {
    spatialized_trajectories = mt_spatialize_tolong(
      trajectories,
      dimensions = c(dimensions),
      n_points = n_points
    )
  }

  if (aggregate_lwd > 0) {
    agg_x = aggregate[,1]
    agg_y = aggregate[,2]
    agg_l = getLength(agg_x, agg_y)
    spatialized_aggregate = mt_spatialize_tolong(
      agg_x, agg_y,
      round(upscale * n_resc * agg_l / l_diag)
    )
  }

  # Compute raw image ----------------------------------------------------------
  if (verbose == TRUE) cat('calculate image','\n')

  # retrieve image
  pts = spatialized_trajectories

  # range of pixels plus white space around image
  xs = 0 : (xres + margin) + 1
  ys = 0 : (yres + ceiling(smooth_radius * 2) + ceiling(smooth_radius * 2)) + 1

  # Remove points outside of bounds
  pts = pts[pts[,1] >= bounds[1] &
            pts[,1] <= bounds[3] &
            pts[,2] >= bounds[2] &
            pts[,2] <= bounds[4],]


  # Determine pixel locations
  x  = round(((pts[,1] - bounds[1]) / px_size) + 1)
  y  = round(((pts[,2] - bounds[2]) / px_size) + 1)

  # Determine table of pixels
  #img_df = data_frame('x'=x,'y'=y)
  #img_tb = img_df %>% group_by(x,y) %>% tally() %>% ungroup()
  img_tb = mousetrap:::tab(x, y)

  # Map pixels into matrix of zeros
  img_mat = matrix(0, ncol=length(xs), nrow=length(ys))
  img_mat[as.matrix(img_tb[,2:1]) + ceiling(smooth_radius * 2)] = img_tb[,3]

  # Store raw image, pixel locations and
  raw_img = c(t(img_mat))
  xys     = expand.grid(1:length(xs), 1:length(ys))

  # Calculate overlay information and create ovrlay image
  if (length(dimensions) == 3) {
    a = pts[,3]
    if (any(is.na(a))) {
      a[is.na(a)] = 0
      message('NAs replaced by 0')
    }

    #img_df = data_frame('x'=x,'y'=y,'a'=a)
    #img_tb = img_df %>% group_by(x,y) %>% summarize(a = mean(a)) %>% ungroup()
    img_tb = tab_sum(x,y,a)

    img_mat = matrix(0, ncol=length(xs), nrow=length(ys))
    img_mat[as.matrix(img_tb[,2:1]) + smooth_radius * 2] = img_tb[,3]
    add_img = c(t(img_mat))
  } else {
    add_img = rep(1, length(raw_img))
  }

  # get aggregate points
  if (aggregate_lwd > 0) {
    agg_x = round(((spatialized_aggregate[,1] - bounds[1]) / px_size) + 1) + smooth_radius * 2
    agg_y = round(((spatialized_aggregate[,2] - bounds[2]) / px_size) + 1) + smooth_radius * 2
    test  = (agg_x > 0 & agg_x <= max(xs) &
             agg_y > 0 & agg_y <= max(ys))
    agg_x = agg_x[test]
    agg_y = agg_y[test]
    agg   = data.frame(
      'x'=agg_x, 'y'=agg_y,
      'col'=aggregate_col, 'lwd'=aggregate_lwd
    )
  }

  # Smooth image ---------------------------------------------------------------

  smooth_img = raw_img
  if (smooth_radius > 0) {
    if(verbose == TRUE)
      cat('smooth image','\n')

    smooth_img = gaussBlur(
      smooth_img, smooth_img,
      max(xs), max(ys),
      smooth_radius
    )

    if (length(dimensions) == 3) {
      add_img = gaussBlur(
        add_img, add_img,
        max(xs), max(ys),
        smooth_radius
      )
    }
  }

  # Create, normalize, and enhance image ---------------------------------------
  # Low-pass: shave off max color intensities
  # Enhance contrast
  # Normalize image

  # create image object
  img = data.frame(xys, smooth_img, add_img)
  names(img) = c('x','y','img','a')

  # Low-pass
  img$img[img$img > low_pass * upsample] = low_pass  * upsample

  # Normalize image
  img$img = (img$img - min(img$img)) / max(img$img)

  if (length(dimensions) == 3) {
    img$a = (img$a - min(img$a)) / max(img$a)
  }

  # Enhance image
  if (auto_enhance == TRUE) {
    ms = c(); steps = c(.1, 2, 5, 10, 20, 50, 100)

    for (i in steps) {
      ms = c(ms, abs(mean(abs(img$img-1)**i - 1)));
    }

    enhance = splinefun(ms, steps)(mean_intensity)
    if (enhance > max(steps)) enhance = max(steps)
    if (enhance < 0) enhance = 1
    img$img = abs(abs(img$img - 1)**enhance - 1)

    if (verbose == TRUE) cat('enhance image by', round(enhance, 1), '\n')
    }

  # Determine colors -----------------------------------------------------------

  img$img = group(img$img, n_shades[1])
  bg = par()$bg
  if (bg == 'transparent') bg = 'white'

  if (length(dimensions) == 2) {
    img$col = colormixer(
      bg, colors[1], img$img, format='hex'
    )
  }

  if (length(dimensions) == 3) {
    if (length(colors) < 2 | length(n_shades) < 2) {
      stop('Colors and n_shades must be of length 2')
    }

    img$a = group(img$a, n_shades[2])
    color_tone = colormixer(colors[1], colors[2], img$a)
    img$col = colormixer(bg, color_tone, img$img, format='hex')
  }

  # Export raw data
  if (aggregate_lwd > 0) {
    return(list('img' = img[,c('x', 'y', 'img', 'col')], 'agg' = agg))
  } else {
    return(img[,c('x', 'y', 'img', 'col')])
  }
}

#' Plot trajectory heatmap
#' 
#' \code{mt_heatmap} plots high resolution raw trajectory maps.
#' 
#' @inheritParams mt_time_normalize
#' @inheritParams mt_heatmap_raw
#' @param file a character string giving the name of the file. If \code{NULL}
#'   the R standard device is used for plotting. Otherwise, the plotting 
#'   device is inferred from the file extension. Only supports devices 
#'   \code{tiff()}, \code{png()}, \code{pdf()}.   
#' @param plot_bounds 
#'
#' @export

mt_heatmap = function(
  data,
  use = 'trajectories',
  dimensions = c('xpos', 'ypos'),
  
  # plot arguments
  file      = 'image.tiff',
  bounds    =  NULL,
  xres      = 1000,
  upscale   = 1,
  upsample  = 1,
  smooth_radius  = 1.5,
  low_pass  = 200,
  auto_enhance = TRUE,
  mean_intensity = .1,
  plot_dims = TRUE,
  norm = TRUE,
  
  # color arguments
  colors      = c('black', 'blue'),
  n_shades    = c(1000, 10),
  
  # plot aggregate
  aggregate_lwd = 0,
  aggregate_col = 'black',
  
  # subsample arguments
  n_trajectories = 10000,
  seed = 1,
  
  # control
  verbose = TRUE  
  
){
  
  # --------- collect device
  if(!is.null(file)){
    device = substr(file,nchar(file)-3,nchar(file))
    if(!device %in% c('pdf','png','tiff')){
      stop('File != NULL requires one of .pdf, .png, or .pdf as file extension.')
      } 
    }
  
  # take time
  t = proc.time()[3]
  
  # --------- create heatmap
  heatmap = mt_heatmap_raw(
    data,
    use = 'trajectories',
    dimensions = c('xpos', 'ypos'),
    bounds    = bounds,
    xres      = xres,
    upscale   = upscale,
    upsample  = upsample,
    smooth_radius  = smooth_radius,
    low_pass  = low_pass,
    auto_enhance = auto_enhance,
    mean_intensity = mean_intensity,
    norm = norm,
    colors      = colors,
    n_shades    = n_shades,
    aggregate_lwd  = aggregate_lwd,
    aggregate_col  = aggregate_col,
    n_trajectories = n_trajectories,
    seed    = 1,
    verbose = TRUE)

  if(!is.data.frame(heatmap)){
    img = heatmap[[1]]
    agg = heatmap[[2]]
    } else {
    img = heatmap
    }
  
  
  if (verbose == TRUE) {
    cat('creating heatmap: ', max(img$x), 'x', max(img$y), 'px', '\n')
  }
  
  # --------- collect device
  if (device == 'pdf') {
    pdf(
      file,
      width=10 * (max(img$x) / max(img$y)) * upscale,
      height=10 * upscale
      )
    } else if(device == 'png') {
    tiff(
      file,
      width=max(img$x) * upscale,
      height=max(img$y) * upscale
    )
  } else if (device == 'tiff') {
    tiff(
      file,
      width=max(img$x) * upscale,
      height=max(img$y) * upscale
    )
  }
  
  # get dimensions
  range_x = range(img$x)
  range_y = range(img$y)
  
  # remove 0s
  p_img = img
  p_img = p_img[img$img>0,]
  
  # set canvas
  plot.new()
  par( mar=c(0, 0, 0, 0) )
  plot.window(xlim=range_x + c(-.5, .5),
              ylim=range_y + c(-.5, .5)
              )
  par(usr = bounds[c(1,3,2,4)])
  
  # plot points
  rect(
    p_img$x - .5, p_img$y - .5,
    p_img$x + .5, p_img$y + .5,
    col=p_img$col,
    border=NA
    )
  
  if (aggregate_lwd > 0) {
    lines(agg[,1:2], col=agg[1,3], lwd=agg[1,4])
    }
  
  # plot dimensions
  if(plot_dims){
    d_x  = diff(range_x)
    d_y  = diff(range_y)    
    xpos = range_x[1] + d_x * c(.02,.98)
    ypos = range_y[1] + d_y * c(.02,.98)
    cord = round(expand.grid(xpos,ypos))
    text(cord,labels = paste0('(',cord[,1],',',cord[,2],')'))
  }
  
  if(device %in% c('pdf','png','tiff')) {
    dev.off()
    }

  # Finalization ---------------------------------------------------------------
  # Give feedback
  if (verbose == T) {
    t = proc.time()[3] - t
    cat('heatmap created in ', round(t), 's\n', sep='')
  }
}
  


# Plot trajectory heatmap
# 
# \code{mt_heatmap} plots high resolution raw trajectory maps.
# 
# @inheritParams mt_time_normalize
# @inheritParams mt_heatmap_raw
# @param file a character string giving the name of the file. If \code{NULL}
#   the   
# 
#
# @export


mt_heatmap_plot_ggplot = function(...) {
  plot_data <- mt_heatmap(..., plot=F)
  
  return(
    ggplot2::ggplot(
      ggplot2::aes_string(x='x', y='y'),
      data=plot_data
    ) +
    ggplot2::scale_x_continuous(
      expand=c(0,0), limits=range(plot_data$x)
    ) +
    ggplot2::scale_y_continuous(
      expand=c(0,0), limits=range(plot_data$y)
    ) +
    ggplot2::geom_raster(
      fill=plot_data$col
    ) +
    ggplot2::theme(
      panel.grid=element_blank(),
      panel.border=element_blank(),
      plot.margin=unit(c(0,0,0,0), "lines"),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    ) +
    labs(x=NULL, y=NULL)
  )
}