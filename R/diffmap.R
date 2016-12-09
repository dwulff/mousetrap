#' Creates difference-heatmap of two trajectory heatmap images
#'
#' \code{mt_diffmap} creates a difference-heatmap of the trajectory
#'   data using gaussian smoothing.
#'
#' @param image1 a data frame containing the trajectory image created by
#'   \code{mt_heatmap}. The first three columns are expected to be named "x",
#'   "y", "img" and contain the respective information.
#' @param image2 see \code{image1}. Must match the dimensions of \code{image1}.
#' @param device a character specifying the plotting device. Either \code{tiff},
#'   \code{png}, or \code{pdf}. Defaults to \code{tiff}.
#' @param name a character string specifying the filename
#' @param directory a character string specifying the directory for the image
#'   file. Defaults to \code{get.wd()}.
#' @param upscale a numeric by which the output resolution of the image is
#'   increased or decreased. Only applies if \code{device} is one of
#'   \code{c("tiff", "png", "pdf")}.
#' @param smooth_radius a numeric specifying the standard deviation of the
#'   gaussian smoothing. If zero smoothing is omitted.
#' @param colors a character vector specifying the colors used to color
#'   cases where \code{image1 > image2, image1 ≈ image2, image1 < image2},
#'   respectively. Note that the colors will necessarily be used in that
#'   specific order. Defaults to c("#00863F", "#FFFFFF", "#FF1900") which
#'   specifies a green-white-red color gradient.
#' @param n_shades integer specifying the number of shades for the color
#'   gradient between the first and second, and the second and third color in
#'   \code{colors}.
#' @param plot boolean specifying whether the image should be plotted. If
#'   \code{plot} is \code{FALSE}, then the image will be returned.
#'
#' @author Dirk U. Wulff
#'
#' @return image representing the difference of image1 and image2

mt_diffmap = function(
  image1,
  image2,
  device = 'tiff',
  name = 'image',
  directory = getwd(),
  upscale = 1,
  smooth_radius = 0,
  colors = c("#00863F", "#FFFFFF", "#FF1900"),
  n_shades = 10,
  plot = TRUE
) {

  # Extract image
  if (!is.data.frame(image1)) image1 = image1[[1]]
  if (!is.data.frame(image2)) image2 = image2[[1]]

  # Compute difference
  img = image1
  img$img = image1$img - image2$img

  # Smooth image
  if (smooth_radius > 0) img$img = gaussBlur(
    img$img, img$img,
    max(img$x), max(img$y),
    smooth_radius
  )

  # Normalize data and determine colors
  v = img$img
  v = v / max(abs(v))

  if (n_shades > 0) {
    vp = v[v > 0]
    vm = v[v <= 0]
    vp = group(vp, n_shades)
    vm = group(vm, n_shades)
    v[v > 0]  = vp
    v[v <= 0] = vm
  }

  col1 = colormixer(colors[2], colors[1], abs(v), 'hex')
  col2 = colormixer(colors[2], colors[3], abs(v), 'hex')
  img$col = ifelse(v > 0, col1, col2)
  img$img = v

  # Plot image (if desired)
  if (plot == TRUE) {

    if (device == 'pdf') {
      pdf(
        paste0(directory,'/',name,'.pdf'),
        width=10 * (max(img$x)/max(img$y)) * upscale,
        height=10 * upscale
      )
    } else if (device == 'png') {
      tiff(
        paste0(directory,'/',name,'.png'),
        width=max(img$x) * upscale,
        height=max(img$y) * upscale
      )
    } else if (device == 'tiff') {
      tiff(
        paste0(directory,'/',name,'.tiff'),
        width=max(img$x) * upscale,
        height=max(img$y) * upscale
      )
    }

    plot.new()
    par(mar=c(0, 0, 0, 0))
    plot.window(
      xlim=range(img$x) + c(-.5, .5),
      ylim=range(img$y) + c(-.5, .5)
    )

    # plot points
    rect(
      img$x - .5, img$y - .5,
      img$x + .5, img$y + .5,
      col=img$col, border=NA
    )
  }

  if (device %in% c('pdf','png','tiff') & plot == TRUE) dev.off()
  if (plot == FALSE) return(img)
}
