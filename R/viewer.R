#' Illustrates individual trajectories
#'
#' \code{mt_viewer} plots individual trajectories in a specified order,
#' including trajectory meta-information, such as whether a trial is identified
#' as a change-of-mind trial.
#'
#' @inheritParams mt_time_normalize
#' @param ids a list of character or numeric vectors containing (short)
#'   identifying labels for each trajectory.
#' @param labels a list of boolean vectors describing the absence or presence of
#'   a property of interest.
#' @param sort_variables a list of character or numeric vectors used to sort the
#'   trajectories, respecting the specific order of sort_variables. If
#'   \code{NULL}, the \code{labels} argument will be used to sort the
#'   trajectories.
#' @param decreasing a boolean vector matching the length of either
#'   sort_variables or labels (depending on whether sort_variables is specified)
#'   which specifies for each variable whether the trajectories should be sorted
#'   in an increasing or decreasing order.
#' @param overlay a character string specifying the variable used to color the
#'   trajectory.
#' @param file a character string specifying the output filename of the
#'   created pdf image. If \code{NULL}, the R standard device will be used. In this
#'   case \code{directory} will be ignored. Otherwise, \code{n_pages} will be
#'   ignored.
#' @param pagewise a boolean specifying whether the program should wait for user
#'   input after filling each page of (\code{nrow} times \code{ncol}) plots.
#' @param nrow an integer specifying the number of panels per column.
#' @param ncol an integer specifying the number of panels per row.
#' @param bounds an integer vector of length four specifying the bottom, left,
#'   upper, and right bounds in terms of percent of the distance between the
#'   trajectories' start and end points.
#' @param pch plotting character passed on to \link[graphics]{points}.
#' @param cex chracter size passed on to \link[graphics]{points}.
#' @param lwd line width passed on to \link[graphics]{lines}.
#' @param plot_box boolean specifying whether a box should be drawn around each
#'   panel.
#' @param plot_index a boolean specifying whether trajectories should be colored
#'   according to the indices of the trajectory points.
#' @param spatialize a boolean specifying whether the trajectories should be
#'   spatialized (see \link{mt_spatialize}). Spatialization helps visualize the
#'   trajectory's shape, but loses time information.
#' @param n_points an integer specifying the number of points in the spatialized
#'   trajectory.
#' @param start boolean specifying wheter the trajectories' start points should
#'   be aligned.
#' @param end boolean specifying wheter the trajectories' end points should be
#'   aligned.
#'
#' @author
#'   Dirk U. Wulff <dirk.wulff@gmail.com>
#'
#' @return Nothing, if image is plotted to external device, otherwise
#'   \code{data}.

mt_viewer = function(
  data,
  use = 'trajectories',
  dimensions = c('xpos', 'ypos'),

  # additional variables
  ids = NULL,
  labels = NULL,
  sort_variables = NULL,
  decreasing = NULL,
  overlay = NULL,

  # plot arguments
  file       = NULL,
  pagewise   = TRUE,
  nrow       = 12,
  ncol       = 6,
  bounds     = c(20, 20, 20, 20),
  pch        = 16,
  cex        = .2,
  lwd        = NULL,
  plot_box   = TRUE,
  plot_index = TRUE,

  # spatialize arguments
  spatialize = TRUE,
  rescale = TRUE,
  resolution = 100,

  # align arguments
  start = TRUE,
  end   = TRUE
){

  # Tests ----------------------------------------------------------------------
  # Minimum number of dimensions
  if (!length(dimensions) %in% 2) {
    stop('Dimensions must of length 2!')
  }

  # Get trajectories
  trajectories <- mousetrap:::extract_data(data=data, use=use)

  # Required data
  if (!all(c(dimensions, overlay) %in% dimnames(trajectories)[[2]])) {
    stop('Not all dimensions exist')
  }

  # Align ----------------------------------------------------------------------
  trajectories = mt_align(
    trajectories,
    dimensions=dimensions, coordinates='norm',
    start=start, end=end
  )

  # Spatialize -----------------------------------------------------------------
  # TODO: This needs rethinking -- the scaling factor
  # should be adapted from the heatmap code.
  if (rescale == TRUE) {
    lengths  = mousetrap:::getLengths(
      trajectories[,dimensions[1],],
      trajectories[,dimensions[2],]
    )
    n_points = round(resolution * lengths / sqrt(2))

    if (is.null(overlay)) {
      trajectories = mousetrap:::mt_spatialize(
        trajectories,
        dimensions=dimensions,
        n_points=n_points
      )
    } else {
      trajectories = mousetrap:::mt_spatialize(
        trajectories,
        dimensions=c(dimensions, overlay),
        n_points=n_points
      )
    }
  }

  # Sort data ------------------------------------------------------------------
  if (!is.null(sort_variables) | !is.null(labels)) {
    if (is.null(sort_variables)) {
      sort_variables = labels
    }

    for (i in length(sort_variables):1) {
      trajectories = trajectories[order(sort_variables[[i]], decreasing=decreasing[i]),,]
      for (j in 1:length(ids)) {
        ids[[j]] = ids[[j]][order(sort_variables[[i]], decreasing=decreasing[i]),]
      }
      for(j in 1:length(categories)) {
        labels[[j]] = labels[[j]][order(sort_variables[[i]], decreasing=decreasing[i]),]
      }
    }
  }


  # Plot -----------------------------------------------------------------------

  # Determine layout
  par(mfrow = c(nrow,ncol))
  xlim = c(-1 - bounds[1] / 100, 1   + bounds[3] / 100)
  ylim = c( 0 - bounds[2] / 100, 1.5 + bounds[4] / 100 + .1)
  n_trajectories = dim(data[[use]])[1]

  # Set output device
  if (!is.null(file)) {
    if(!grepl('.pdf',file)) stop('File != NULL requires the file extension .pdf') 
    pdf(file, width=ncol, height=nrow)
  } else {
    # Limit the number of screens shown to five,
    # if the user is showing the data in her environment
    n_trajectories = (nrow * ncol) * 5
  }

  # loop over trajectories
  for (i in 1:n_trajectories) {

    # begin plot
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)

    # Get x and y points warped to mt-space
    x = -trajectories[i,dimensions[1],]
    y =  trajectories[i,dimensions[2],] * 1.5


    # Determine colors
    ind = which(!is.na(trajectories[i,'xpos',]))
    if (plot_index == TRUE) {
      time_intensity = seq(0, 1, length=length(ind))
    } else {
      time_intensity = 0
    }

    if (!is.null(overlay)) {
      overlay_values = trajectories[i,overlay,ind]
      overlay_intensity = overlay_values / max(overlay_values)
    } else {
      overlay_intensity = 0
    }
    col = rgb(time_intensity, 0, overlay_intensity)

    # plot lines
    if (!is.null(lwd)) {
      if (!is.null(overlay) | plot_time == TRUE) {
        for (j in 1:(length(x) - 1)) {
          lines(
            c(x[j], x[j+1]),
            c(y[j], y[j+1]),
            lwd=lwd, col=col[j]
          )
        }
      } else {
        lines(x, y, lwd=lwd)
      }
    }

    # plot points
    if (!is.null(pch)) {
      if (!is.null(lwd)) {
        points(x, y, col=col, cex=cex, lwd=lwd, pch=pch, bg='white')
      } else {
        points(x, y, col=col, cex=cex, lwd=lwd)
      }
    }

    # add ID
    id_pos = c(.5,.05)
    if (!is.null(ids)) {
      id = ''
      for (j in 1:length(ids)) {
        id = paste(id, ids[[j]][i], sep=ifelse(id == '', '', '  '))
      }
      text(
        xlim[1] + id_pos[1], ylim[1] + id_pos[2],
        cex=.6, font=1, adj=0,
        labels=id
      )
    }


    # Add labels
    label_pos = c(0, .1, .2, .3, .4)
    if (!is.null(labels)) {
      for (j in 1:length(labels)) {
        value = labels[[j]][i]
        rect(
          xlim[1] + label_pos[j],     ylim[1],
          xlim[1] + label_pos[j + 1], ylim[1]+.1,
          border='grey90', col=rgb(1-value, 1-value, 1-value)
        )
      }
    }

    # Draw border
    if (plot_box == TRUE) {
      box(lwd=.5, col='black')
    }

    # Wait for user
    if (pagewise == TRUE & i %% (nrow*ncol) == 0) {
      if (readline('Press [enter] to continue or [q] to abort ') == 'q') {
        break
      }
    }
  }

  # Close device if outputting to file
  if (!is.null(file)) dev.off()
}
