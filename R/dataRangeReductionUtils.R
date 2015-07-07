# This file provides functions to reduce the range of data values (as opposed
# to reducing the number of data points)

###
# Data reduction functions: Continuous to category
###

#' Sort values into bins
#'
#' Given a vector of numeric values returns a vector of bin names, where the n
#' bins to sort into are defined by n+1 bin boundaries (including start and end
#' boudaries, probably -Inf and Inf). Basically a convienience wrapper around
#' cut().
#'
#' Bins the data values into n bins as s specified by the n+1 binEnds, then
#' returns a vector of the same length as data, but with each value replaced by
#' the binLabel for the bin it was sorted into. E.g if the bin labels are
#' colors, this converts data into a color scale. Values equal to bin boundaries
#' are put in the lower bin, with the lowest bin boundary also included in the
#' lowest bin.
#'
#' @param data A vector of numeric data to bin. NA values and values above or
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
#' @param binLabels The vector of names corresponding to the labesl of the bin a
#'   value belongs in.
#'
#' @return The vector of bin labels corresponding to the bins the data was
#'   sorted into.
#'
#' @export
sortDataIntoBins <- function(data, binEnds, binLabels) {
   bins = cut( data, binEnds, labels = 1:( length(binEnds) - 1 ), include.lowest = TRUE );
   return( binLabels[ as.numeric( levels(bins)[bins] )]);
}
