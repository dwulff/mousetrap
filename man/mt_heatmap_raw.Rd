% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{mt_heatmap_raw}
\alias{mt_heatmap_raw}
\title{Creates high-resolution heatmap of trajectory data}
\usage{
mt_heatmap_raw(data, use = "trajectories", dimensions = c("xpos", "ypos"),
  bounds = NULL, xres = 1000, upscale = 1, upsample = 1,
  smooth_radius = 1.5, low_pass = 200, auto_enhance = TRUE,
  mean_intensity = 0.1, device = "tiff", colors = c("black", "blue"),
  n_shades = c(1000, 10), aggregate_lwd = 0, aggregate_col = "black",
  n_trajectories = 10000, seed = 1, verbose = TRUE)
}
\arguments{
\item{data}{a mousetrap data object created using one of the mt_import
functions (see \link{mt_example} for details). Alternatively, a trajectory
array can be provided directly (in this case \code{use} will be ignored).}

\item{use}{a character string specifying which trajectory data should be
used.}

\item{dimensions}{a character vector specifying the trajectory variables
used to create the heatmap. The first two entries are used as x and
y-coordinates, the third, if provided, will be added as color information.}

\item{bounds}{numeric vector specifying the corners (xmin, ymin, xmax, ymax) of
the plot region. Note that the trajectories' start and end points are mapped to
(0, 0) and (-1, 1.5), respectively. Bounds should therefore include the rectangle
c(-1, 0, 1, 1.5). Defaults to c(-1.4, -.3, 1.4, 1.9).}

\item{upscale}{a numeric value by which the output resolution of the image is increased or
decreased. Only applies if \code{device} is one of \code{c("tiff", "png", "pdf")}.}

\item{upsample}{a numeric by which the number of points used to represent individual
trajectories are increased or decreased. Values of smaller than one will improve speed
but also introduce a certain level of granularity.}

\item{smooth_radius}{a numeric specifying the standard deviation of the gaussian smoothing.
If zero, smoothing is omitted.}

\item{low_pass}{an integer specifying allowed number of counts per pixel. This
arguments limits the maximum pixel color intensity.}

\item{mean_intensity}{a numeric between 0 and 1 specifying the average color intensity
across the entire image. Defaults to .2.}

\item{n_trajectories}{an integer specifying the number of trajectories used to create the
image. If \code{n_trajectories} is smaller than containes in the trajectorie object
specified by \code{use} the \code{n_trajectories} are randomly sampled.}

\item{seed}{an integer specifying the seed used for the trajectory sampling.}

\item{verbose}{boolean specifying whether progress updates should printed.}

\item{px_size}{an integer specifying the size of px. A px_size of .01 implies an
x-resolution of (xmax-xmin)/px_size and a y-resolution of (ymax-ymin)/px_size.
Defaults to .001, resulting in an image size of 2800x2200 px if the \code{bounds} are unchanged.}

\item{color}{a character string or numeric vector of length three specifying
the color used to plot the third dimension.}

\item{color_order}{a character string, either "increasing" or "decreasing" specifying the
whether large or small values of the third dimension are assigned high color values,
respectively.}
}
\value{
Nothing, when image is plotted using an external device. Otherwise a matrix
   containing the image's pixel information.
}
\description{
\code{mt_heatmap_raw} creates a high-resolution heatmap of the trajectory
  data using gaussian smoothing.
}
\author{
Dirk U. Wulff <dirk.wulff@gmail.com>
}

