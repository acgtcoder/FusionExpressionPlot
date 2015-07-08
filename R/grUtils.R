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

#' Convert a genomic elements string into an GRanges object.
#'
#' Parses a string into a sequence of genomic elements, adding specified
#' chromosome and strand. Each genomic elements is assumed to be two genomic
#' positions separated by a delimiter (-), with the low position first and the
#' high position second, in 1 based coordinates. The sequence of elements are
#' separated by a delimiter (,). All elements must be on the same sequence block
#' (chromosome) and must be on the same strand (by default '*' for unspecified).
#' Positions are given as if on the positive strand, actual positions would be
#' reversed if on the negative strand. Whitespace around delimiters is ignored.
#' If alternate delimiters that are whitespace are used, then they must be
#' present but additional whitespace is ignored.
#'
#' @param elementString The string to parse into genomic ranges.
#' @param chr The chromosome for the genomic ranges.
#' @param strand The strand for the chromosome, must be '*', '+' or '-'.
#'   Defaults to '*' if not specified.
#' @param seqDelim The sequence delimiter, usually use the default '\\s*,\\s*'. Can be
#'   any reqular expression.
#' @param posDelim The start-end position separator, usually use the default
#'   '\\s*-\\s*'. Can be any regular expression.
#'
#' @return A GRanges object based on the input elements.
#'
#' @export
grFromElementString <- function(
   elementString, chr, strand='*',
   seqDelim='\\s*,\\s*', posDelim='\\s*-\\s*'
) {
   elements <- strsplit(elementString, seqDelim)[[1]]
   starts = as.numeric(sapply( sapply( elements,function(x) strsplit(x, posDelim)), "[[", 1))
   ends = as.numeric(sapply( sapply(elements, function(x) strsplit(x, posDelim)), "[[", 2))
   gr <- grNew( start=starts, end=ends, chr=chr,strand=strand)
   return(gr)
}

#' Convert a genomic location string into a GRange object.
#'
#' Given a genomic location string, converts it to a GRanges object. The
#' location string is parsed based on a regular expression with up to three
#' capture groups, capturing the chromosome, the strand, and the element
#' location sub-strings. The element string is then parsed into a sequence of
#' start and stop ranges, and the this data is used to create and return a
#' GRange object. The regexp used to match with must contain at least three
#' named capture groups: \code{(?<chr>..)}, \code{(?<elements>...)}, and
#' \code{(?<strand>..)}. Matching is done using \code{perl=TRUE}).
#'
#' The default location string parsed is assumed to look something like:
#' \preformatted{    "chr1:1-100,200-300,400-500:-"}
#'
#' White-space is ignored around the delimiters ':', ',', and '-'.
#'
#' @param locationString The string to parse into a GRange object
#'
#' @param captureRE The regular expression used to extract the substrings for
#' the chromosome (chr), the element list (elements), and strand (strand).
#'
#' @return A GRange object corresponding to the location string.
#'
#' @export
grFromLocationString <- function(
   locationString, captureRE='(?<chr>.*?)\\s*:\\s*(?<elements>.*?)\\s*:\\s*(?<strand>[-+*]?)'
) {
   matchResults <- regexpr(captureRE, locationString, perl= TRUE)
   parsed <- regexprNamedMatches(matchResults, locationString)
   gr <- grFromElementString( parsed[1,'elements'],
                            chr=    parsed[1,'chr'],
                            strand= parsed[1, 'strand']
   )
   return(gr)
}