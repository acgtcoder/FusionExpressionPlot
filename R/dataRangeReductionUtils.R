# This file provides functions to reduce the range of data values (as opposed
# to reducing the number of data points)

###
# Data reduction functions: Continuous to category
###

#' Map labels to values
#'
#' Maps labels to a vector of numeric data based on a
#' vector of bin ends. Each label corresponds to one of n bins defined by n+1
#' bin boundaries (including start and end boudaries, probably -Inf and Inf).
#' Returns a vector of the same length as data, but with each value replaced by
#' the binLabel for the bin it was sorted into. E.g if the bin labels are
#' categories like low, medium and high, this converts data to categories. Values equal to bin boundaries
#' are put in the lower bin, with the lowest bin boundary also included in the
#' lowest bin.
#'
#' Basically this is a convienience wrapper around cut().
#'
#' @param x A vector of numeric data to bin. NA values and values above or
#'   below the lowest bin end are allowed, they convert to NA's in the returned
#'   label vector.
#'
#' @param binEnds The n+1 ends of the n bins into which the data is sorted. One
#'   bin would be c(-Inf, Inf). Two bins ~ equally weighted would be c(-Inf,
#'   median(data), Inf). This vector will be sorted before assigning labels, so
#'   reverse the label order instead of providing bins as high to low. Bin ends
#'   belong to the lower bin, with the lowest bin end also part of the lowest
#'   bin.
#'
#' @param binLabels The vector of names corresponding to the labels of the bin a
#'   value belongs in.
#'
#' @return The vector of bin labels corresponding to the bins the data was
#'   sorted into.
#'
#' @export
mapLabels <- function(data, binEnds, binLabels) {
   bins = cut( data, binEnds, labels = 1:( length(binEnds) - 1 ), include.lowest = TRUE );
   return( binLabels[ as.numeric( levels(bins)[bins] )]);
}

#' Map colors to values
#'
#' Maps colors to a vector of numeric data based on a pallette of colors and a
#' vector of bin ends. Can use either a named pallete from the color brewer
#' package or a list of explicit colors. Data is binned and labeled with the
#' associated color as per \code{\link{sortDataIntoBins}}. This does not
#' produce smooth gradient, but instead compresses the information in the data
#' into a discrete set of colors.
#'
#' This is really just a convienience function to allow using color brewer
#' pallettes easily when binning data. Colors are no different than other
#' lables and can be applied generically using mapLabels
#'
#' @section Color Brewer Palletes:
#'
#'   If no colors are specified, a vector of colors will be assigned using the
#'   brewer palette named (RdBu if none specified). The number of colors needed
#'   will be determined automatically from the binEnds provided. The color brewer
#'   palette is just a convienient shortcut to specifying a list of colors and
#'   is ignored if a vector of colors is specified. If no version of the given
#'   palette has the number of levels required by the bin ends, an error will
#'   occur.
#'
#'   For diverging values, can have 3 - 11 bins, allowed values are:
#'
#'   \preformatted{BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral}
#'
#'   For sequential values, can have 3 - 8 value, allowed values are:
#'
#'   \preformatted{Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples
#'   RdPu Reds YlGn YlGnBu YlOrBr YlOrRd}
#'
#'   For qualitative palletes, can have 3 to (n) where n varies by palette:
#'
#'   \preformatted{Accent (8), Dark2 (8), Paired (12), Pastel1 (9), Pastel2 (8), Set1 (9),
#'   Set2 (8), Set3 (12)}
#'
#'   Note: palletes are only specified in one directon, and binEnds are always
#'   sorted from lowest to highest before using, so set reverse=TRUE if you want
#'   to reverse the color order.
#'
#' @examples
#' # Display available palletes.
#' library(RColorBrewer);
#' display.brewer.all()
#'
#' # Display colors in a pallette
#' palletteName <- "RdBu"
#' levels <- 4
#' brewer.pal(levels, palletteName)
#'
#' x <- c(-100, -1, 0, 1, 100)
#'
#' # Map data to 3 colors using the default (RdBu) pallette
#' mapColors( x, binEnds= c(-Inf,-10,10,Inf))
#'
#' # Map data to 3 colors using the reverse of the default (RdBu) pallette
#' mapColors( x, binEnds= c(-Inf,-10,10,Inf), reverse=TRUE)
#'
#' # Map data using an explicit color vector
#' colorMap <- ('red', 'violet', 'blue')
#' mapColors( x, binEnds= c(-Inf,-10,10,Inf), colors= colorMap)
#'
#' @param x The vector of numeric data to map to a color
#'
#' @param binEnds The sequence of cut-points defining the bins into which the
#'   data will be split. Each bin will correspond to a color in the color
#'   vector or palette and provides for mapping each value to a
#'   color. Values exactly matching bin ends will be in the lower bin.
#'   NAs or values outside the range will map to NAs. Normally the first and
#'   last bin ends are the special values Inf and -Inf. Note: bin ends will be
#'   sorted from smallest to largest before mapping, so make sure lables are
#'   ordered to match.
#'
#' @param brewerPaletteName The name of the RColorBrewer palatte to use for bin
#' labels, by default 'RdBu'. Setting labels causes this to be ignored.
#'
#' @param colors A vector of colors to use as bin names. The number of colors will
#' have to be one less than the number of bin ends. If this is specified, any
#' brewerPaletteName is ignored
#'
#' @param reverse Reverse the order of the colors If specifying the colors
#' directly, could just reverse them before specifying them. However, when using
#' the color brewer pallette, might want to reverse the pallette order. The
#' default value is FALSE.
#'
#' @return the vector of colors the data maps to.
#'
#' @export
mapColors <- function(
   x, binEnds, brewerPaletteName= 'RdBu', reverse=FALSE, colors=NULL
) {
   if (! is.null(colors)) {
      labels <- colors
   } else {
      labels= brewer.pal( length(binEnds) - 1, brewerPaletteName)
   }
   if (reverse) {
      labels <- rev(labels)
   }
   labelCol <- mapLabels( x, binEnds, labels);
   return(labelCol);
}
