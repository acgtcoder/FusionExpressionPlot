# This file provides functions to reduce the range of data values (as opposed
# to reducing the number of data points)

###
# Data reduction functions: Continuous to category
###

#' Maps values into labeled bins
#'
#' Basically this is a convenience wrapper around \code{\link{cut}}.
#'
#' Maps labels to a vector of numeric data based on a vector of bin ends. Each
#' label corresponds to one of n bins defined by n+1 bin boundaries (including
#' start and end boundaries, probably \code{-Inf} and \code{Inf}). It is an
#' error if the number of bins and the number of labels differ. Returns a vector
#' of the same length as data, but with each value replaced by the binLabel for
#' the bin it was sorted into. E.g if the bin labels are categories like "low",
#' "medium" and "high", this converts data to categories. Values equal to bin
#' boundaries are put in the lower bin, with the lowest bin boundary also
#' included in the lowest bin. Note, bin ends are always sorted from lowest to
#' highest. If you want the labels in the other order, just reverse the labels
#' vector.
#'
#' @param data A vector of numeric data to bin. Where values are \code{NA},
#'   above the highest bin end, or below the lowest they convert to \code{NA} in
#'   the returned label vector.
#'
#' @param binEnds The n+1 ends of the n bins into which the data is sorted. One
#'   bin might be \code{c(-Inf, Inf)}. Two bins equally weighted would be
#'   \code{c(-Inf, median(data), Inf)}. This vector will be sorted before
#'   assigning labels, so reverse the label order instead of providing bins as
#'   high to low. Data exactly matching bin ends will be in the lower bin, with
#'   the lowest bin end also part of the lowest bin. It is an error if any bin
#'   end is repeated.
#'
#' @param binLabels The vector of names corresponding to the labels of the bin a
#'   value belongs in.
#'
#' @return The vector of bin labels corresponding to the bins the data was
#'   sorted into.
#'
#' @section Errors:
#'    \describe{
#'       \item{
#'          \code{'breaks' are not unique}
#'       }{
#'          It is an error if any bin end is repeated.
#'       }
#'    \item{
#'       \code{lengths of 'breaks' and 'labels' differ}
#'    }{
#'       It is an error if the number of labels differs from the number of
#'       bins, or equivalently is not one less than the number of bin ends.
#'    }
#' }
#'
#' @seealso
#'    \code{\link{cut}}
#'
#' @examples
#' ends = c(Inf,-Inf,1,-1)    # Sorted to c(-Inf, -1, 1, Inf)
#' labels = c('-', '0', '+' )
#' mapLabels( c(-Inf,-5,-1,-0.5,0,0.5,1,5,Inf), ends, labels)
#' #=> [1] "-", "-", "-", "0", "0", "0", "0", "+", "+"
#'
#' ends = c(-10, 0, 10)
#' labels = c('-', '+')
#' mapLabels( c(-20, -10, -5, 0, NA, 5.2, 10, Inf), ends, labels)
#' #=> [1] NA  "-" "-" "-" NA  "+" "+" NA
#'
#' @export
mapLabels <- function(data, binEnds, binLabels) {
   if (length(binLabels) != length(binEnds) -1 ) {
      stop("lengths of 'breaks' and 'labels' differ")
   }
   bins = cut( data, binEnds, labels = 1:( length(binEnds) - 1 ), include.lowest = TRUE );
   return( binLabels[ as.numeric( levels(bins)[bins] )]);
}

#' Map colors to values
#'
#' Maps colors to a vector of numeric data based on a palette of colors and a
#' vector of bin ends. Can use either a named palette from the
#' \pkg{RColorBrewer} package or a list of explicit colors. Data is binned and
#' labeled with the associated color as per \code{\link{mapLabels}}. This is not
#' intended to produce a smooth gradient, but instead compresses the information
#' in the data into a discrete set of colors.
#'
#' This is really just a convenience function to allow using color brewer
#' palettes easily when binning data. Colors are no different than other labels
#' and can be applied generically using \code{\link{mapLabels}}. Note, this
#' requires the \pkg{RColorBrewer} package if using \code{brewerPaletteName}
#'
#' @section Color Brewer Palettes:
#'
#'   If no colors are specified, a vector of colors will be assigned using the
#'   brewer palette named, \code{'RdBu'} if not specified. Each palette comes
#'   in multiple variants with a different number of colors. The number of
#'   colors will be determined automatically from the \code{binEnds} provided.
#'   Using color brewer palettes is just a convenient shortcut to specifying a
#'   list of colors and is ignored if a vector of colors is specified. If no
#'   version of the given brewer palette has the number of levels required by
#'   the bin ends, an error will occur.
#'
#'   For diverging values, can have 3 - 11 bins, allowed values are:
#'
#'   \preformatted{BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral}
#'
#'   For sequential values, can have 3 - 8 bins, allowed values are:
#'
#'   \preformatted{Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd
#' Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd}
#'
#'        (Note: it is Greys, not Grays)
#'
#'   For qualitative palettes, can have 3 to (n) bins where n varies by palette:
#'
#'   \preformatted{Accent (8), Dark2 (8), Paired (12), Pastel1 (9), Pastel2 (8),
#' Set1 (9), Set2 (8), Set3 (12)}
#'
#'   Note: palettes are only specified in one direction, and binEnds are always
#'   sorted from lowest to highest before using, so set reverse=TRUE if you want
#'   to reverse the color order.
#'
#' @param x The vector of numeric data to map to a color
#'
#' @param binEnds The sequence of cut-points defining the bins into which the
#'   data will be split. Each bin will correspond to a color in the color vector
#'   or palette and provides for mapping each value to a color. Values exactly
#'   matching bin ends will be in the lower bin. \code{NA}s or values outside
#'   the range will map to \code{NA}. Normally the first and last bin ends are
#'   the special values \code{Inf} and \code{-Inf}. Note: bin ends will be
#'   sorted from smallest to largest before mapping, so make sure labels are
#'   ordered to match.
#'
#' @param brewerPaletteName The name of the \pkg{RColorBrewer} palette to use
#'   for bin labels, by default \code{'RdBu'}. Setting labels causes this to be
#'   ignored.
#'
#' @param colors A vector of colors to use as bin names. The number of colors
#'   will have to be one less than the number of bin ends. If this is specified,
#'   any \code{brewerPaletteName} is ignored and the \pkg{RColorBrewer} package
#'   is not needed
#'
#' @param reverse Reverse the order of the colors. If explicitly specifying the
#'   colors though \code{colors}, one could instead just reverse them. However,
#'   when using \code{brewerPaletteName}, this is the  only way to reverse the
#'   palette order. The default value is \pkg{FALSE}.
#'
#' @return The vector of colors the data maps to.
#'
#' @seealso
#'    \code{\link{mapLabels}}
#'
#' @examples
#' library(RColorBrewer);
#'
#' ### Review from RColorBrewer
#' # Display available palettes
#' display.brewer.all()
#'
#' # Display colors in a palette
#' paletteName <- "RdBu"
#' levels <- 4
#' RColorBrewer::brewer.pal(levels, paletteName)
#'
#' ### Using mapColors
#' x <- c(-100, -1, 0, 1, 100)
#' ends <- c(-Inf,-10,10,Inf)
#'
#' # Map data to 3 colors using the default (RdBu) palette
#' mapColors( x, binEnds= ends )
#' #=> [1] "#EF8A62" "#F7F7F7" "#F7F7F7" "#F7F7F7" "#67A9CF"
#'
#' # Map data to 3 colors using the reverse of the default (RdBu) palette
#' mapColors( x, binEnds= ends, reverse= TRUE)
#' #=> [1] "#67A9CF" "#F7F7F7" "#F7F7F7" "#F7F7F7" "#EF8A62"
#'
#' # Map to a specified palette name
#' mapColors( x, binEnds= ends, brewerPaletteName= 'Set1' )
#' #=> [1] "#E41A1C" "#377EB8" "#377EB8" "#377EB8" "#4DAF4A"
#'
#' # Map data using an explicit color vector
#' colorMap <- c('red', 'violet', 'blue')
#' mapColors( x, binEnds= ends, colors= colorMap)
#' #=> [1] "red"    "violet" "violet" "violet" "blue"
#'
#' mapColors( x, binEnds= ends, colors= colorMap, reverse= TRUE)
#' #=> [1] "blue"   "violet" "violet" "violet" "red"
#'
#' @export
mapColors <- function(
   x, binEnds, brewerPaletteName= 'RdBu', reverse=FALSE, colors=NULL
) {
   if (! is.null(colors)) {
      labels <- colors
   } else {
      labels= RColorBrewer::brewer.pal( length(binEnds) - 1, brewerPaletteName)
   }

   if (reverse) {
      labels <- rev(labels)
   }

   labelCol <- mapLabels( x, binEnds, labels);
   return(labelCol);
}
