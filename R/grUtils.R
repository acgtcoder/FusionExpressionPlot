# This file provides utility functions that add to or wrap the behavior of
# GenomicRanges package GRange objects.

#' Create a simple GRanges object.
#'
#' Create a simple Granges object based on a sequence of exons on one strand of
#' one chromosome. This could be a gene, a transcript, or one side of a fusion.
#' Defaults to '*' if strand is not specified. Coordinates are given as if the
#' exons were on the positive strand, in 1-based coordinates. This will be
#' backawards from the actual sequence if the gene is on the negative strand.
#' Ranges will be sorted by start and then by end.
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
   return(sort(gr))
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

#' Convert a genomic element string into an GRanges object.
#'
#' Parses a genomic elements string into ranges, then adds a specified
#' chromosome and strand and returns a GRanges object. The elements string parsed
#' by default looks like:
#' \preformatted{    "1-100,3-500,2-3" }
#'
#' Spaces around the between delimiter ',' and the from delimiter '-' are ignored;
#' other delimiters (regular expressions) can be specified.  All elements must be on the same sequence block
#' (chromosome) and must be on the same strand (by default '*' for unspecified).

#' Positions are given in 1 based coordinates as if on the positive strand with the low position first and the
#' high position second. Actual positions would be
#' reversed if on the negative strand.
#'
#' @param elementString The string to parse into genomic ranges.
#' @param chr The chromosome for the genomic ranges.
#' @param strand The strand for the chromosome, must be '*', '+' or '-'.
#'   Defaults to '*' if not specified.
#' @param betweenDelim The sequence delimiter, usually use the default '\\s*,\\s*'. Can be
#'   any regular expression.
#' @param fromDelim The start-end position separator, usually use the default
#'   '\\s*-\\s*'. Can be any regular expression.
#'
#' @return A GRanges object based on the input elements, with ranges sorted by start and then
#' end position
#'
#' @export
grFromElementString <- function(
   elementString, chr, strand='*',
   betweenDelim='\\s*,\\s*', fromDelim='\\s*-\\s*'
) {
   elements <- strsplit(elementString, betweenDelim)[[1]]
   starts = as.numeric(sapply( sapply( elements,function(x) strsplit(x, fromDelim)), "[[", 1))
   ends = as.numeric(sapply( sapply(elements, function(x) strsplit(x, fromDelim)), "[[", 2))
   gr <- grNew( start=starts, end=ends, chr=chr,strand=strand)
   return(sort(gr))
}

.grGetChr <- function( gr ) {
   chrFactor <- unique(seqnames(gr))
   if (length(chrFactor) > 1) {
      stop("GRange object has more than one chromosome.")
   } else if (length(chrFactor) == 0) {
      return('')
   } else {
      return(as.character(chrFactor[1]))
   }
}

.grGetStrand <- function( gr ) {
   strandFactor <- unique(strand(gr))
   if (length(strandFactor) > 1) {
      stop("GRange object has more than one strand.")
   } else if (length(strandFactor) == 0) {
      return('')
   } else {
      return(as.character(strandFactor[1]))
   }
}

#' Extract a genomic element string from a GRange object
#'
#' Returns a string representation of the ranges in a GRange object if they are
#' on only one chromosome of one strand. The two ends of each range are
#' separated by a '-' by default but any string can be specified for this 'from'
#' delimiter. Additionally, each range element is separated by a ',' by default,
#' but any string can be specified for this 'between' delimiter. If there are no
#' ranges in the GRange object, then an empty string is returned.
#'
#' It is an error to try and generate an element string from a GRanges object
#' representing ranges on more than one chromosome or on more than one strand.
#'
#' This is the reciprocal of \code{\link{grFromElementString}}.
#'
#' @examples
#' gr <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chr1', strand='-')
#' grToElementString(gr)
#' #=> [1] "1-100,200-300,400-500"
#' grToElementString(gr, fromDelim=' to ', betweenDelim='; ')
#' #=> [1] "1 to 100; 200 to 300; 400 to 500"
#' gr1 <- grFromElementString(grToElementString(gr), chr='ANY')
#' grToElementString(gr1) == grToElementString(gr)
#' #=> [1] TRUE
#'
#' @param gr The GRange object to extract an element string from.
#'
#' @param fromDelim The string to print between the start and end position of a range. Defaults to '-'
#'
#' @param betweenDelim The string to print between the range elements. Defaults to ','
#'
#' @return A string equivalent of the ranges in the GRange object provided.
#'
#' @export
grToElementString <- function( gr, fromDelim='-', betweenDelim=',' ) {
   # These called only for error handling
   .grGetChr(gr)
   .grGetStrand(gr)
   if (length(gr) == 0) {
      return('')
   }
   starts <- start(gr)
   ends <- end(gr)
   return( paste( starts, ends, sep=fromDelim, collapse=betweenDelim ));
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
#' @param ... Allows passing alternate fromDelim and betweenDelim regular
#' expressions through to grFromElementString, which parses the element part of
#' the locationstring
#'
#' @return A GRange object corresponding to the location string.
#'
#' @export
grFromLocationString <- function(
   locationString, captureRE='(?<chr>.*?)\\s*:\\s*(?<elements>.*?)\\s*:\\s*(?<strand>[-+*]?)',
   ...
) {
   matchResults <- regexpr(captureRE, locationString, perl= TRUE)
   parsed <- regexprNamedMatches(matchResults, locationString)
   gr <- grFromElementString( parsed[1,'elements'],
                            chr=    parsed[1,'chr'],
                            strand= parsed[1, 'strand'],
                            ...
   )
   return(gr)
}

#' Extract a genomic location string from a GRange object
#'
#' Returns a string representation of a GRange object if it contains only ranges
#' on just one chromosome and on just one strand. The string returned is
#' specified using the format characters \{\{c\}\} for the chromosme, \{\{S\}\}
#' for the strand and \{\{e\}\} for the element string (formatted as with
#' grToElementString()). The default format returned is
#' "\{\{c\}\}:\{\{e\}\}:\{\{s\}\}"
#'
#' It is an error to try and generate a location string from a GRanges object
#' representing ranges on more than one chromosome or on more than one strand.
#'
#' This is the reciprocal of /code{/link{grFromLocationString}}.
#'
#' @examples
#' gr <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chr1', strand='-')
#' grToLocationString(gr)
#' #=> [1] "chr1:1-100,200-300,400-500:-"
#' grToLocationString(gr, format="{{c}}({{s}}){{e}}")
#' #=> [1] "chr1(-)1-100,200-300,400-500"
#' grToLocationString(gr, format="On {{c}}: {{e}}.", fromDelim=' to ', betweenDelim=', ')
#' #=> [1] "On chr1: 1 to 100, 200 to 300, 400 to 500."
#' gr1 <- grFromLocationString(grToLocationString(gr))
#' grToLocationString(gr1) == grToLocationString(gr)
#' #=> [1] TRUE
#'
#' @param gr The GRange object to extract a location string from.
#'
#' @param format The string to print, translating \{\{e\}\} to an element
#'   location string, \{\{c\}\} as the chromosme, \{\{s\}\} as the strand.
#'   Defaults to "\{\{c\}\}:\{\{e\}\}:\{\{s\}\}",
#'
#' @param ... Allows passing parameters to the \code{\link{grToElementString}}
#'   used to generate the elements string.
#'
#' @return A string equivalent to the GRange object provided, formatted as
#'   specified.
#'
#' @export
grToLocationString <- function( gr, format= '{{c}}:{{e}}:{{s}}', ... ) {
   if (length(gr) == 0) {
      return('')
   }
   c <- .grGetChr(gr)
   s <- .grGetStrand(gr)
   e <- grToElementString( gr, ... )
   return( templateFill(format) )
}
