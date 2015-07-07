# This file provides utility functions that add to or wrap the behavior of
# GenomicRanges package GRange objects.

#' Create a simple GRanges object.
#'
#' Create a simple Granges object based on a sequence of exons on one strand of
#' one chromosome. This could be a gene, a transcript, or one side of a fusion.
#' Defaults to '*' if strand is not specified. Coordinates are given as if the
#' exons were on the positive strand, in 1-based coordinates. This will be
#' backawards from the actual sequence if the gene is on the negative strand.
#'
#' @param start A vector giving the starting (lower) genomic coordinates of the
#'   exons. Must be specified.
#' @param end A vector giving the ending (higher) genomic coordinates of the
#'   exons. Must be specified.
#' @param chr The chromosome (singular) on which these exons occur, like 'chr1'
#'   or 'chrX'. Must be specified.
#' @param strand The strand the exons are actually on, as '+', '-', or '*'.
#'   Defaults to '*' if not specified.
#' @return The genomicRanges object built from the exons.
#'
#' @export
grNew <- function( start, end, chr, strand= '*' ) {

   gr <- GRanges(
      ranges= IRanges( start=start, end=end ),
      seqnames= chr,
      strand= strand
   );
   return(gr)
}

#' Add a data column to a GRanges object.
#'
#' Adds a vector as a column in a GRanges object, assuming the data is in the
#' correct order. This data can be anything describing the range elements, for
#' instance a vector of colors. This is essentially a convienience wrapper
#' around the elementMetadata() function from the GRanges package.
#'
#' @param gr The granges object to add to add a column to.
#'
#' @param name The name of the column when added
#'
#' @param vec The data. If not the same length as the GRanges object, will
#'   attempt to wrap values. If not an even multiple, will error.
#'
#' @return The granges object, with the appended column.
#'
#' @export
grAddColumn <- function( gr, name, vec ) {
   elementMetadata(gr)[[name]] <- vec;
   return(gr);
}

###
# GRanges object functions
###

#' Add a category column to a GRange object
#'
#' Adds a column of labels as metadata to a GRanges object, mapped from an
#' existing column of numeric metadata. Which column to map can be specified by
#' column name or index. The data from that column will be binned based on the
#' given binEnds and the matching bin label will be used. By default assumes you
#' are assinging a color map using a colorBrewere pallete, but any category
#' labels can be used. To add a category column based on an external vector,
#' just use sortDataIntoBins and grAddColumn.
#'
#' @section Color Brewer Palletes:
#'
#'   If no labels are specified, a vector of colors will be assigned using the
#'   brewer palette named (RdBu if none specified). The number of colors needed
#'   will be determined automatically from the binEnds provided. To use
#'   different colors, just specify them manually as bin labels. Doing so is no
#'   different than assigning any lables, colors or otherwise. The color brewer
#'   palette is just a convienient shortcut and is ignored if labels are given.
#'   If no version of the given palette has the number of levels required by the
#'   bin ends an error will occur.
#'
#'   For diverging values, can have 3 - 11 bins, allowed values are
#'
#'   BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral.
#'
#'   For sequential values, can have 3 - 8 value, allowed values are:
#'
#'   Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples
#'   RdPu Reds YlGn YlGnBu YlOrBr YlOrRd.
#'
#'   For qualitative palletes, can have 3 to (n) where n varies by palette:
#'
#'   Accent (8), Dark2 (8), Paired (12), Pastel1 (9), Pastel2 (8), Set1 (9),
#'   Set2 (8), Set3 (12)
#'
#' @examples
#' # Display available palletes.
#' library(RColorBrewer);
#' display.brewer.all()
#'
#' @param gr The GRanges object to annotate.
#'
#' @param column The metadata column in the GRange object to determine labels
#' from. Can specify the column by name or by (metadata) column index.
#'  map. Defaults to 1.
#'
#' @param binEnds The sequence of cut-points defining the bins into which the
#'   data will be split. Each bin will correspond to a label in the provided
#'   labels vector or color palette, providing a mapping of each value to a
#'   label or color. Values exactly matching bin ends will be in the lower bin.
#'   NAs or values outside the range will map to NAs. Normally the first and
#'   last bin ends are the special values Inf and -Inf. Note: bin ends will be
#'   sorted from smallest to largest before mapping, so make sure lables are
#'   ordered to match
#'
#' @param brewerPaletteName The name of the RColorBrewer palatte to use for bin
#' labels, by default 'RdBu'. Setting labels causes this to be ignored.
#'
#' @param labels The vector of labels to use as bin names. The number of colors will
#' have to be one less than the number of bin ends. If specified, brewerPaletteName is ignored
#'
#' @param as.column The name of the metadata column added to the GRanges object.
#'
#' @return the GRanges object with the new column added, or if annotate is set
#'   FALSE, just the label (color) vector.
#'
#' @export
grLabelMap <- function(
   gr, binEnds, column=1, brewerPaletteName= 'RdBu',
   labels= rev(brewer.pal( length(binEnds) - 1, brewerPaletteName)),
   as.column = 'exonColor', annotate=TRUE
) {
   colorCol <- sortDataIntoBins(elementMetadata(gr)[[column]],binEnds,colors);
   if (annotate) {
      return(grAddColumn(gr, as.column, colorCol));
   }
   else {
      return(colorCol);
   }
}
