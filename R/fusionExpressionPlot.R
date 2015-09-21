# library(GenomicRanges);
# library(RColorBrewer);

#' Plot Expression Profiles for Fused Genes
#'
#' Plots exons from two genes as boxes scaled by exon size and colored based on
#' expression in a sample. Fusion points between the two genes are represented
#' as triangles at the appropriate position. This allows visualizing the drop in
#' relative expression of the exons that are part of a putative expressed fusion
#' transcript relative to the exons only in the unfused gene transcripts.
#' Normalization is performed across an entire cohort of samples and then again
#' within each gene, so color differences show exons whose expression is
#' unusually variant within the genes against a background of similar samples.
#'
#' @docType package
#' @name FusionExpressionPlot
#'
#' @import IRanges
#' @import GenomicRanges
#' @import RColorBrewer
NULL

###
# Utility functions: testing.
###

# Die with "message" if "test" is false.
test_ok <- function( test, message="Test failed!." ) {
   if (! test) { stop( message ) }
}

# Die with "message" if "got" and "want" are not identical.
test_is <- function( got, want, message="" ) {
   test_ok(
      (! identical(got, want)),
      message= paste( c("\n", "Failed test: ", message, "\n", "Got: \"", got, "\" but wanted: \"", want, "\""), sep= "", collapse="" )
   )
}


###
# GRanges object factories
###

#' Extract GRange objects from cohort exon expression data
#'
#' Given a dataframe describing exon expression for a sample cohort, generate a
#' list of GRange objects for one sample, by selected genes. Each GRange object
#' is identified by gene name and includes a meta-data column giving the
#' expression data.
#'
#' @section Cohort exon expression data:
#'
#'   The input dataframe must have the following columns:
#'
#'   \describe{
#'       \item{chr}{
#'   The chromosome this exon is on, as 'chr1'...'chr22', chrX', 'chrY', or
#'   'chrM'.}
#'       \item{start}{
#'   The position of the first base of the exon.}
#'       \item{end}{
#'   The position of the last base of the exon.}
#'       \item{strand}{
#'   The strand of this exon (+,-, or *). All exons in a gene should have the
#'   same strand.}
#'       \item{gene}{
#'   The name of this gene. All exons in a gene should have the same name, and
#'   no exons from different genes should have the same name.}
#'       \item{***}{
#'   The normalized expression value for this exon in sample named "***". Should
#'   have one column for every sample in the cohort. Sample names may not be
#'   "chr", "start", "end", "strand", "gene", and must be unique.}
#'   }
#'
#'   Start and end are assumed to be 1 based, but this does not really matter as
#'   long as all positions here and in any position data (i.e. fusions) used
#'   with this data all have the same base.
#'
#' @param df A dataframe with gene model and exon expression. Required to have
#'   the following columns: 'chr', 'start', 'end', 'strand', 'gene', and one
#'   column for each sample.
#'
#' @param genes The names of the genes to extract from the cohort data frame.
#'
#' @param sample The name of the sample to extract from the df. Required
#'
#' @param column The metadata column name under which the expression data is
#'   added to the GRange objects. By default this is 'exonData'
#'
#' @return A List of GRanges object for the gene, with start and end positions
#'   of the exon specified, and with a metadata column "relativeExpression"
#'   giving the relative expression for the selected gene in the selected
#'   sample.
#'
#' @export
grFromCohortDf <- function( df, genes, sample, column='exonData' ) {
   grList = list()
   for (gene in genes ) {
      geneDf <- df[df$gene == gene,c('chr', 'start', 'end', 'strand', 'gene', sample)];

      gr <- GRanges(
         ranges= IRanges( start= geneDf$start, end= geneDf$end),
         seqnames= geneDf$chr,
         strand= geneDf$strand
      );
      gr <- grAddColumn(gr, column, geneDf[,sample]);
      grList[[gene]] = gr
   }
   return(grList);
}

###
# Specialized functions for use in relative expression graphs
###

# Add a column of colors to gr based on relativeExpression column
expressionGrLabelMap <- function( gr ) {
   exonFillColor <- mapColors( gr$relativeExpression,
                               binEnds= c(Inf,0.3,0.1,-0.1,-0.3,-Inf),
                               brewerPaletteName= 'RdBu', reverse=TRUE
   )
   gr <- grAddColumn(gr, 'exonFillColor', exonFillColor)
   return(gr)
}

#' Convert from genomic feature to model coordinates.
#'
#' Converts a vector of (+ strand) genomic positions to a vector of coordinates
#' in a model as described by a GRanges object. The model describes a sequence
#' of non-ovelapping elements separated by gaps of varying sizes. Genomic
#' positions are maped 1-1 within elements, but scale to a fixed gap width when
#' they fall between elements. The GRanges object describing the element model
#' must be all on one strand of one chromosome. '-' strand models are reveresed
#' so all coordinates are 5' - 3' transformed.
#'
#' @section How to convert to model coordinates:
#'
#'   Model coordinates start with '1', corresponding to the genomic position of
#'   the first element in the GRanges object. Model coordinate end with a value
#'   that is the sum of all the element widths plus the number of gaps * the
#'   gapthWidth. Genomic positions before the model (less than any element
#'   start) are converted to model coordinate 0; genomic positions after the
#'   mode (greater than any element end are converted to the last model
#'   coordinate + 1. This can be controlled by the outOfBounds parameter (Not
#'   implemented).
#'
#'   If a genomic position falls within a gap, it is converted to a gap-relative
#'   position, then shifted 1/2 a position in the 5' direction. This shifted,
#'   gap-relative position [0.5 ... actualGapWidth - 0.5] is then scaled to fit
#'   (multiplied by gapWidth/actualGapWidth). The next largest integer (ceiling)
#'   is taken as the offset into the gap in the model, with the returned model
#'   coordinate being "the sum of the widths of the exons before" + (gapWidth *
#'   "the number of gaps before") + "this offset". If the strand is negative,
#'   the coordinate is relative to the last exon base of the mode instead of the
#'   first.
#'
#' @section The conversion formula:
#'
#'   For positive strand coordinate values:
#'
#'   \deqn{ coordinate =   sum( widths(elementsBefore)) + count(gapsBefore) *
#'   gapWidth + ceiling( (p - endPos(elementJustBefore) - 0.5) *
#'   (gapWidth/thisGapWidth)) }
#'
#'   For negative strand coordinate values, an additional transformation is
#'   applied:
#'
#'   \deqn{ coordinate =  sum( widths( allElements )) + count(allGapsB) *
#'   gapWidth + 1 - posCoordinate }
#'
#' @param pos    Vector of integers (required).
#'
#'   A vector of chromosome genomic positions to convert to model coordinates.
#'   Assumed to be on the same chromosome as the provided gene model, although
#'   not required to be within the gene model.
#'
#' @param gene  GRanges object (required).
#'
#'   Object giving the model, i.e. gene exons. All elements must be on the same
#'   chromosome on the same strand in sequence order with no overlaps and no 0
#'   length introns. Assumes '-' strand genes are described as if their shadow
#'   on the '+' strand was the real gene, with ony the '-' indicateing anything
#'   aobut a reverse order.
#'
#' @param gapWidth  Integer (default = 100)
#'
#'   Fixed gap width. Equivalent is size to an exon of the same length.
#'
#' @param outOfBounds 'warning' (default) | 'error' | 'ok'
#'
#'   By default, genomic positions can be supplied that are outside the model,
#'   but a warning will be given, once. "Position outside the gene model found:
#'   #". Setting outOfBounds='error' will cause the first such position found to
#'   abort. Setting outOfBounds='ok' means no warning or error will be
#'   generated.
#'
#' @return Integer vector of genomic positions converted to model coordinates.
#'
#' @examples
#' library(GenomicRanges);
#' grPos <- GRanges(
#'    ranges= IRanges(
#'       start= c( 101, 301, 511, 513, 813 ),
#'       end= c( 200, 500, 511, 612, 822 )
#'    ),
#'    seqnames= c( "chr1" ),
#'    strand= c( "+" )
#' );
#'
#' @export
genomicToModelCoordinates <- function(
   pos, gene, gapWidth= 100, outOfBounds='warning'
) {
   #
   # TODO: Valdate gene
   # Check that gene is GRanges object, single chromosome, single stranded,
   # reduced, and ordered.

   # No overlaps, sorted.
   gene = reduce(gene);

   # Initialize return coordinates
   modelCoordinate = rep(NA, length(pos));

   eStarts = start(gene);
   eEnds = end(gene);
   eWidths = width(gene)
   iWidths = width(gaps(gene))[-1] # Drop gap" before first exon.
   lastModelPos = sum(eWidths) + length(iWidths) * gapWidth
   strand = as.character(strand(gene)[1]);

   # loop over each pos in input. TODO: convert to apply?
   for (i in 1:length(pos)){
      p = pos[i]

      # Handle point before gene model, or after, if '-' strand model
      if ( p < eStarts[1] ) {

         # User is to be notified when position is outside gene model
         if (outOfBounds != 'ok') {
            message = paste0("Position outside the gene model found: ", p);
            if (outOfBounds == 'warn' || outOfBounds == 'warning') {
               warning( message );
            }
            if (outOfBounds == 'die' || outOfBounds == 'error') {
               stop( message );
            }
         }

         # Depending on strand, return 1 past end or one before start.
         modelCoordinate[i] <- if (strand == '-') lastModelPos + 1 else 0;
      } # END choice = before model (after if '-')


      # Handle point after gene model, or before, if '-' strand model
      else if (p > eEnds[length(gene)]) {

         # User is to be notified when position is outside gene model
         if (outOfBounds != 'ok') {
            message = paste0("Position outside the gene model found: ", p);
            if (outOfBounds == 'warn' || outOfBounds == 'warning') {
               warning( message );
            }
            if (outOfBounds == 'die' || outOfBounds == 'error') {
               stop( message );
            }
         }

         # Depending on strand, return  0 or 1 past end.
         modelCoordinate[i] = if (strand == '-') 0 else lastModelPos + 1;
      } # END choice = after model (before if '-')

      # Handle point within model
      else {

         # Contribution to model coordinate of prior complete exons and introns
         # TODO: inefficient as sorted
         eBefore= max( c( 0, which( p > eEnds )));
         iBefore= max(c(which( p >= eStarts ))) - 1;
         wholePart = sum( eWidths[0:eBefore]) + iBefore * gapWidth;

         # Contribution to model coordinate of partial exon or intron pos p is in.
         fractionPart = 0
         # In exon eBefore + 1
         if (eBefore == iBefore) {
            fractionPart = p - eStarts[eBefore + 1] + 1;

            # In intron iBefore + 1
         } else if (eBefore == iBefore + 1) {
            fractionPart = ceiling((p - eEnds[eBefore] - 0.5) * gapWidth/(iWidths[iBefore + 1]));

         } else {
            stop( paste0("Bad program!, eBefore= ", eBefore, "; iBefore= ", iBefore, "; p= ", p, "."));
         }

         # Coordinate
         modelCoordinate[i] <- wholePart + fractionPart;
         if (strand == '-') {
            modelCoordinate[i] <- lastModelPos - modelCoordinate[i] + 1
         }
      } # END choice = within model
   } # END loop = input genomic positions vector
   return( modelCoordinate );
}

#' Plotable rectangles from GRanges Object
#'
#' Given a genomicRanges Object, return the data needed to plot it as a sequence
#' of rectangles based on widths. Gaps are all fixed at a single intronWidth.
#' Height, fill color, and border color for the exon rectangles can be specified
#' on a per-exon basis (vector-wrapped), or given in the genomicRanges object
#' (as metadata colums, by  default exonHeight, exonFillColor, and
#' exonBorderColor). Returns the rectangle as a list of data, as specified
#' below. Uses the genomicToModelCoordinates function, passing throught the
#' parameters needed
#'
#' @param gr The GenomicRanges object to plot.
#'
#'   Object giving the gene to calculate the parameters for. Passed through to
#'   genomicToModelCoordinates as its gene parameter. Besides the criteria
#'   specified there, may also contain data columns 'exonHeight',
#'   'exonFillColor' and 'exonBorderColor'. If present, the values in these
#'   columns will be used instead on the equivalently names parameters to this
#'   function.
#'
#' @param intronWidth The fixed spacing to use for gaps between rectangles.
#'
#'   The width of the itron, the x spacing between the exon rectangles. Passed
#'   through to genomicToModelCoordinates as its gapWidth parameter. Defaults to
#'   100.
#'
#' @param exonHeights The vector of heights to use for the rectangles.
#'
#'   If fewer than the number of elements are givien, will be wrapped. If not
#'   specified exonHeightsCol must contain the column name in the GRanges object
#'   to use instead, by default it looks for column 'exonHeight'.
#'
#' @param exonFillColors The vector of fill colors to use for the rectangles.
#'
#'   If fewer than the number of elements are givien, will be wrapped. If not
#'   specified exonFillColors must contain the column name in the GRanges object
#'   to use instead, by default it looks for column 'exonFillColor'.
#'
#' @param exonBorderColors The vector of border colors to use for the
#'   rectangles.
#'
#'   If fewer than the number of elements are givien, will be wrapped. If not
#'   specified exonBorderColors must contain the column name in the GRanges
#'   object to use instead, by default it looks for column 'exonBorderColor'.
#'
#' @param exonHeightsCol The name of GRanges column to use for exonHeights, if
#'   not specified. By default 'exonHeight'
#'
#' @param exonFillColorsCol The name of GRanges column to use for
#'   exonFillColors, if not specified. By default 'exonFillColor'
#'
#' @param exonBorderColorsCol The name of GRanges column to use for
#'   exonBorderColors, if not specified. By default 'exonBorderColor'
#'
#' @return A list with the following elements:
#'
#'   \describe{
#'       \item{xStarts}{The low x coordinate; the left side of each rectangle.}
#'       \item{xEnds}{The high x coordinate; the right side of each rectangle.}
#'       \item{yBottoms}{The low y coordinate; the bottom of each rectangle.}
#'       \item{yTops}{The high y coordinate; the top of each rectangle.}
#'       \item{fillColors}{The fill color for each rectangle.}
#'       \item{fillColors}{The border color for each rectangle.}
#'       \item{xRange}{c(<min x value>, <max x value).}
#'       \item{yRange}{c(<min y value>, <max y value).}
#'   }
#'
#' @export
grToRect <- function (
   gr, intronWidth= 100,
   exonHeights= 100, exonFillColors= 'blue', exonFillColorsCol= 'exonFillColors',
   exonBorderColors= 'black', exonBorderColorsCol= 'exonBorderColors',
   exonHeightsCol='exonHeights'
) {

   rect <- list();
   rect$onReverseStrand = (as.character(strand(gr)[1]) == '-');

   if (rect$onReverseStrand) {
      rect$xStarts  <- rev(genomicToModelCoordinates( end(  gr), gr, intronWidth ));
      rect$xEnds    <- rev(genomicToModelCoordinates( start(gr), gr, intronWidth ));
      rect$yBottoms <- rep( 0, length( gr ));
      rect$yTops        <- rev(elementMetadata(gr)[['exonHeight']]);
      rect$fillColors   <- rev(elementMetadata(gr)[['exonFillColor']]);
      rect$borderColors <- rev(elementMetadata(gr)[['exonBorderColor']]);
   }
   else {
      rect$xStarts  <- genomicToModelCoordinates( start(gr), gr, intronWidth );
      rect$xEnds    <- genomicToModelCoordinates( end(  gr), gr, intronWidth );
      rect$yBottoms <- rep( 0, length( gr ));
      rect$yTops        <- elementMetadata(gr)[['exonHeight']];
      rect$fillColors   <- elementMetadata(gr)[['exonFillColor']];
      rect$borderColors <- elementMetadata(gr)[['exonBorderColor']];
   }

   if ( is.null( rect$yTops        )) { rect$yTops        <- rep_len( exonHeights,      length( gr )); }
   if ( is.null( rect$fillColors   )) { rect$fillColors   <- rep_len( exonFillColors,   length( gr )); }
   if ( is.null( rect$borderColors )) { rect$borderColors <- rep_len( exonBorderColors, length( gr )); }

   rect$xRange <- c( min(rect$xStarts),  max(rect$xEnds));
   rect$yRange <- c( min(rect$yBottoms), max(rect$yTops));

   return(rect);
}

# Given a gene-exon plotable object, shift it by x and y bp
offsetRect <- function ( rect, xAdd=0, yAdd=0 ) {
   rect$xStarts  <- rect$xStarts + xAdd;
   rect$xEnds    <- rect$xEnds   + xAdd;
   rect$xRange   <- rect$xRange  + xAdd;

   rect$yBottoms <- rect$yBottom + yAdd;
   rect$yTops    <- rect$yTops   + yAdd;
   rect$yRange   <- rect$yRange  + yAdd;

   return(rect)
}

demo.grToRect <- function() {

   main = 'Defaults';
   sub = 'Gene model for TMPRSS2 (-)';
   gr <- demo.makeGrForTmprss2();
   rect <- grToRect( gr )
   plot( rect$xRange, rect$yRange, type="n", sub=sub, main=main);
   rect( rect$xStarts, rect$yBottoms, rect$xEnds, rect$yTops,
         col= rect$fillColors, border= rect$borderColors
   );

   main = 'Internal exon plot parameters';
   sub = 'Gene model for Test (+)';
   gr <- demo.makeGrForTest();
   rect <- grToRect( gr )
   plot( rect$xRange, rect$yRange, type="n", sub=sub, main=main);
   rect( rect$xStarts, rect$yBottoms, rect$xEnds, rect$yTops,
         col= rect$fillColors, border= rect$borderColors
   );

   main = 'Pair of genes';
   sub = 'Test( + ) and TMPRSS2 (-)';
   gr1 <- demo.makeGrForTest();
   gr2 <- demo.makeGrForTmprss2();
   rect1 <- grToRect( gr1 )
   rect2 <- offsetRect( grToRect( gr2 ), 500 + rect1$xRange[2] );
   plot( c( rect1$xRange[1], rect2$xRange[2] ),
         c( min(rect1$yRange[1], rect2$yRange[1] ), max(rect1$yRange[2], rect2$yRange[2])),
         type="n", sub=sub, main=main);
   rect( rect1$xStarts, rect1$yBottoms, rect1$xEnds, rect1$yTops,
         col= rect1$fillColors, border= rect1$borderColors
   );
   rect( rect2$xStarts, rect2$yBottoms, rect2$xEnds, rect2$yTops,
         col= rect2$fillColors, border= rect2$borderColors
   );

}

#' Convert points to polygons
#'
#' Converts a vector of coordinate points to a data object that can be drawn
#' with \code{\link{polygon}}. It calculate the values need to plot the points
#' as polygons specified by shape. Currently only triangles pointing up and
#' triangles pointing down can be specified. Each point can be described
#' separately by passing in a vector for the \code{shapes} parameter. Other than
#' x and y coordinates of the points to plot and the two scale paramters,
#' provided values will be wrapped as needed.
#'
#' @section Shape scaling:
#'
#'   The shape of a pixel is not necessarily square, so a square of 10 x 10
#'   pixels may not be square but instead rectangular. To make the polygon
#'   shapes come out symetrical, the relative heights and widths of the pixels
#'   must be specified. This is the purpose of the \code{scale} and
#'   \code{resScale} parameters. Two parameters are provided to allow specifying
#'   the plot scaling ratio and the screen resolution scaling ratio separately.
#'
#' @param x The x plot coordinates where the the poly-points should be drawn
#'
#' @param y The y plot coordinates where the the poly-points should be drawn
#'
#' @param widths The widths in plot coordinates of the poly-points. Default =
#'   5.
#'
#' @param heights The height in plot coordinates of the poly-points. Default =
#'   5.
#'
#' @param sides The position of the poly-points over the coordinate point. Can
#'   have a value from 1-9, oriented like a 3 x 3 touch pad, 1 has uper left
#'   corner of poly as the point, 5 has the middle over the point, and 9 has the
#'   lower right as the point. Default = 5, centered.
#'
#' @param shapes The shape of the poly-point. Allowed shapes are
#'   \code{'triangleDown'} (the default) and \code{'triangleUp'}.
#'
#' @param fillColors The colors filling the poly-point. Default = \code{'red'}.
#'
#' @param borderColors The colors of the poly-point borders. Default =
#'   \code{'black'}.
#'
#' @param scale The scale factor of the x and y axis of the plot, i.e if plot
#'   is 100x500, scale should be \code{c(1,5)} or \code{c(0.2,1)}. By default
#'   this is \code{c(1,1)}.
#'
#' @param resScale The scale factor of the device, i.e. if the screen
#'   resolution is 1024x768, this should be \code{c(1,1024/768)} or
#'   \code{c(768/1024,1)}. By default this is \code{c(1,1)}.
#'
#' @return Returns a list with the data needed to plot the points as polygons.
#'   The elements of the list are:
#'
#' \tabular{ll}{
#'    \code{count} \tab number of points plotted as polygons\cr
#'    \code{xRadii} \tab 1/2 widths of points, after wrapping\cr
#'    \code{yRadii} \tab 1/2 heights of points, after wrapping\cr
#'    \code{shapes} \tab The shapes parameter, after wrapping\cr
#'    \code{sides} \tab The sides parameter, after wrapping\cr
#'    \code{fillColors} \tab The fillColors parameter, after wrapping\cr
#'    \code{borderColors} \tab The borderColors parameter, after wrapping\cr
#'    \code{x} \tab X coordinates of polygons to draw, \code{NA} separated\cr
#'    \code{y} \tab Y coordinates of polygons to draw, \code{NA} separated\cr
#'    \code{xRange} \tab lowest and highest x value drawn\cr
#'    \code{yRange} \tab lowest and highest y values drawn\cr
#' }
#'
#' @export
pointToPoly <- function (
   x, y, widths=20, heights=20, shapes= 'triangleDown', sides= 5,
   fillColors='red', borderColors='black',
   scale=c(4,1), resScale=c(1,1)
) {
   if (length(x) != length(y)) {
      stop( "x and y lengths differ" );
   }

   poly <- list()
   numOfPoints <- length(x);
   poly$count        <- numOfPoints
   poly$xRadii       <- rep_len( widths,       length.out = numOfPoints ) * scale[1] * resScale[1] / 2;
   poly$yRadii       <- rep_len( heights,      length.out = numOfPoints ) * scale[2] * resScale[2] / 2;
   poly$shapes       <- rep_len( shapes,       length.out = numOfPoints );
   poly$sides        <- rep_len( sides,        length.out = numOfPoints );
   poly$fillColors   <- rep_len( fillColors,   length.out = numOfPoints );
   poly$borderColors <- rep_len( borderColors, length.out = numOfPoints );

   px = numeric(0);
   py = numeric(0);
   for (i in 1:numOfPoints) {
      sideXchange <- c(      poly$xRadii[i], -1 * poly$xRadii[i],              0 )[ (poly$sides[i] %%  3) +1 ];
      sideYchange <- c( -1 * poly$yRadii[i],                   0, poly$yRadii[i] )[ (poly$sides[i] + 2) %/% 3];

      if (poly$shapes[i] == 'triangleDown') {
         px = c( px, c( x[i],                  x[i] - poly$xRadii[i], x[i] + poly$xRadii[i] ) + sideXchange, NA);
         py = c( py, c( y[i] - poly$yRadii[i], y[i] + poly$yRadii[i], y[i] + poly$yRadii[i] ) + sideYchange, NA);
      }
      else if (poly$shapes[i] == 'triangleUp') {
         px = c( px, c( x[i],                  x[i] - poly$xRadii[i], x[i] + poly$xRadii[i] ) + sideXchange, NA);
         py = c( py, c( y[i] + poly$yRadii[i], y[i] - poly$yRadii[i], y[i] - poly$yRadii[i] ) + sideYchange, NA);
      }
      else {
         stop(paste0( "shape ", poly$shapes[i], " is unknown to pointToPoly." ));
      }
   }

   poly$x = px
   poly$y = py
   poly$xRange=c( min(px, na.rm = TRUE), max(px, na.rm = TRUE) );
   poly$yRange=c( min(py, na.rm = TRUE), max(py, na.rm = TRUE) );

   return( poly );
}

demo.pointToPoly <- function() {
   main = 'Defaults';
   sub = 'Plot: 300 x 300; Points: (100,200) - 5x5, (200,100) - 5x5.';
   pol <- pointToPoly(c(100,200), c(200,100));
   plot( c(0,300), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors);

   main = 'Resolution Scaling 1680x1050';
   sub = 'Plot: 300 x 300; Points: (100,200) - 20x20, (200,100) - 20x20.';
   pol <- pointToPoly(c(100,200), c(200,100), resScale=c(1,1680/1050));
   plot( c(0,300), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors);

   main = 'Options: shapes, fillColors, and borderColors';
   sub = 'Plot: 300 x 300; Points: (100,200) - 5x5, (200,100) - 5x5.';
   pol <- pointToPoly(c(100,200), c(200,100), shapes=c('triangleDown', 'triangleUp'),
                      fillColors=c('blue','green'), borderColors= c('green', 'red') );
   plot( c(0,300), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors);

   main = 'Asymetric points, symetric plot (not normally what is wanted)';
   sub = 'Plot: 300 x 300; Points: (100,200) - 20x5, (200,100) - 5x20.';
   pol <- pointToPoly( c(100,200), c(200,100), widths=c(20, 5), heights=c(5, 20) );
   plot( c(0,300), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors);

   main = 'Asymetric plot, unscaled points (not normally what is wanted)';
   sub = 'Plot: 1200 x 300; Points: (100,200) - 5x5, (200,100) - 5x5.';
   pol <- pointToPoly( c(100,200), c(200,100) );
   plot( c(0,1200), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors);

   main = 'Asymetric plot, points scaled by c(4,1)';
   sub = 'Plot: 1200 x 300; Points: (100,200) - 5x5, (200,100) - 5x5.';
   pol <- pointToPoly( c(100,200), c(200,100), scale=c(4,1));
   plot( c(0,1200), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors );

   main = 'Equivalent but not conceptually clear - asymetric shapes';
   sub = 'Plot: 1200 x 300; Points: (100,200) - 20x5, (200,100) - 20x5.';
   pol <- pointToPoly( c(100,200), c(200,100), widths=c(20,20), heights=c(5,5));
   plot( c(0,1200), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors );

   main = 'Y positioning: sides = 1:9';
   sub = 'Plot: 300 x 300; Points every 20 from x=100 at y = 100.';
   pol <- pointToPoly(c(100,120,140,160,180,200,220,240,260),
                      c(100,100,100,100,100,100,100,100,100),
                      sides=c(1,2,3,4,5,6,7,8,9), resScale=c(1,1680/1050),
                      widths=10, heights=10);
   plot( c(0,300), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors);

   main = 'X positioning: sides = 1:9';
   sub = 'Plot: 300 x 300; Points every 20 from y=100 at x = 100.';
   pol <- pointToPoly(c(100,100,100,100,100,100,100,100,100),
                      c(100,120,140,160,180,200,220,240,260),
                      sides=c(1,2,3,4,5,6,7,8,9), resScale=c(1,1680/1050),
                      widths=10, heights=10);
   plot( c(0,300), c(0,300), type="n", sub=sub, main=main);
   polygon( pol$x, pol$y, col=pol$fillColors, border=pol$borderColors);

}

###
# Plotting
###

#' Generate the plot data object for a fusion expression plot.
#'
#' Merge data from multiple genes and multiple point sets into a list for
#' plotting as a fusion expression graph. Uses \code{\link{grToRect}} and
#' \code{\link{pointToPoly}}. This is a lower-level function called by
#' plotting functions like \code{\link{plotFusionExpressionOne}} or
#' \code{\link{plotFusionExpressionPair}}
#'
#' @seealso \code{\link{grToRect}}, \code{\link{pointToPoly}},
#' \code{\link{plotFusionExpressionOne}}, \code{\link{plotFusionExpressionPair}}
#' \code{\link{plotFusionExpressionPlotData}}
#'
#' @param genes List of genomic ranges objects to plot, by gene name.
#'
#' @param fusions List of vectors of genomic fusion points to plot, by gene
#'   name.
#'
#' @param geneNames Vector of gene names.
#'
#' @param geneSep Distance between genes in the plot, in bases. Default = 500.
#'
#' @param fusionSep Distance between exon tops and tops of plotted fusion
#'   points, in bases. Default = 40
#'
#' @param nameSep Distance between the plot and the gene name annotation, in
#'   bases. Default = 20.
#'
#' @return A "fusion expression plot data" object, a list with the data needed
#'   to plot the gene models, expression levels, and fusion points. This list
#'   has the folowing elements:
#'
#' \tabular{ll}{
#'    geneRects   \tab The coordinates of the rectangles to plot as exons\cr
#'    fusionPolys \tab The coordinates of the polygons (triangles) used as fusion point markers\cr
#'    geneNames   \tab The names of the genes being plotted\cr
#'    xRange      \tab The min and max x values in the plot\cr
#'    yRange      \tab The min and max y values in the plot\cr
#' }
#'
#' @export
getFusionExpressionPlotData <- function( genes, fusions, geneNames,
                                         geneSep= 500, fusionSep= 40, nameSep= 20) {
   ###
   #   genes = list of GRanges objects, by gene name
   #   fusions = list of vectors of fusion genomic positions, by gene name
   #   geneNames = list of gene names, must match those in genes and fusions
   ###

   geneRects = list();
   fusionPolys = list();
   xOffset = 0;
   xMin = c();
   xMax = c();
   yMin = c();
   yMax = c();
   for (gene in geneNames) {
      gr <- genes[[gene]];
      geneRect <- grToRect( gr );
      geneRects[[gene]] <- offsetRect( geneRect, xAdd= xOffset, yAdd= nameSep );
      xMin <- c( xMin, geneRects[[gene]]$xRange[1] );
      xMax <- c( xMax, geneRects[[gene]]$xRange[2] );
      yMin <- c( yMin, geneRects[[gene]]$yRange[1] - nameSep );  # Leaving nameSep space on bottom
      if (! is.null(fusions)) {
         genomePos <- fusions[[gene]];
         xFusion <- genomicToModelCoordinates( genomePos, gr );
         xFusion <- xFusion + xOffset;
         yFusion <- rep(geneRect$yRange[2] + fusionSep, length(xFusion));
         fusionPolys[[gene]] <- pointToPoly( xFusion, yFusion );
         yMax <- c( yMax, fusionPolys[[gene]]$yRange[2] );
      } else {
         fusionPolys <- NULL
         yMax <- c( yMax, geneRects[[gene]]$yRange[2] );
      }
      xOffset <- xOffset + geneRect$xRange[2] + geneSep;
   }
   xRange = c(min(xMin), max(xMax));
   yRange = c(min(yMin), max(yMax));


   return( list(geneRects=geneRects,fusionPolys=fusionPolys, geneNames=geneNames, xRange=xRange, yRange=yRange) );
}


demo.getFusionExpressionPlotData <- function() {

   geneNames <- c('clown', 'TMPRSS2');
   genes <- list();
   genes[[geneNames[1]]] <- demo.makeGrForTest();
   genes[[geneNames[2]]] <- demo.makeGrForTmprss2();
   fusions <- list()
   fusions[[geneNames[1]]] <- c(101, 400, 822);
   fusions[[geneNames[2]]] <- c(42836479, 42837000, 42879992);
   return( getFusionExpressionPlotData( genes, fusions, geneNames ));

}

demo.getFusionExpressionPlotData4 <- function() {

   geneNames <- c('clownA', 'TMPRSS2', 'clownB', 'clownC' );
   genes <- list();
   genes[[geneNames[1]]] <- demo.makeGrForTest();
   genes[[geneNames[2]]] <- demo.makeGrForTmprss2();
   genes[[geneNames[3]]] <- demo.makeGrForTest();
   genes[[geneNames[4]]] <- demo.makeGrForTest();

   fusions <- list()
   fusions[[geneNames[1]]] <- c(101, 400, 822);
   fusions[[geneNames[2]]] <- c(42836479, 42837000, 42879992);
   fusions[[geneNames[3]]] <- c(101, 400);
   fusions[[geneNames[4]]] <- c(400, 822);
   return( getFusionExpressionPlotData( genes, fusions, geneNames ));

}

demo.getFusionExpressionPlotData2 <- function() {

   geneNames <- c('TMPRSS2a', 'TMPRSS2b');
   genes <- list();
   genes[[geneNames[1]]] <- demo.makeGrForTmprss2();
   genes[[geneNames[2]]] <- demo.makeGrForTmprss2();
   fusions <- list()
   fusions[[geneNames[1]]] <- c(42837000, 42879992);
   fusions[[geneNames[2]]] <- c(42836479, 42837000);
   return( getFusionExpressionPlotData( genes, fusions, geneNames ));

}

demo.getFusionExpressionPlotData1 <- function() {

   geneNames <- c( 'TMPRSS2' );
   genes <- list();
   genes[[geneNames[1]]] <- demo.makeGrForTmprss2();
   fusions <- list()
   fusions[[geneNames[1]]] <- c(42837000, 42879992);
   return( getFusionExpressionPlotData( genes, fusions, geneNames ));

}


#' Plot a fusionExpressionPlotData object
#'
#' Plots the data object generated by \code{\link{getFusionExpressionPlotData}}
#' and saves it as a pdf file. This is a lower-level function called by plotting
#' functions like \code{\link{plotFusionExpressionOne}} or
#' \code{\link{plotFusionExpressionPair}}
#'
#' @param data A data object generated by
#'   \code{\link{getFusionExpressionPlotData}}.
#'
#' @param title The title to put on the graph.
#'
#' @param file The filename (pdf) to save the plot as. By default this is
#'   \code{temp.pdf}
#'
#' @param scale The number of "bases" per inch in the plot, default = 1000
#'
#' @param ratio The height/width ratio of the plot.
#'
#' @param scaleOffset The distance in bases that the scale annotation is set
#'   below the gene exon graph.
#'
#' @return nothing
#'
#' @seealso \code{\link{getFusionExpressionPlotData}}
#'   \code{\link{plotFusionExpressionOne}}
#'   \code{\link{plotFusionExpressionPair}}
#'
#' @export
plotFusionExpressionPlotData <- function ( data, title, file= 'temp.pdf',
                                            scale=1000, ratio =4, scaleOffset = 15
) {
   geneNames <- data$geneNames;
   geneRects <- data$geneRects;
   fusionPolys <- data$fusionPolys;
   numGenes = length(geneNames);
   xRange = data$xRange;
   yRange = data$yRange;
   yRange[1] = yRange[1] - scaleOffset
   if (length(geneRects) != numGenes ) { stop( "gene names don't match geneRects data" ); }
   if (! is.null(fusionPolys) && length(fusionPolys) != numGenes ) { stop( "gene names don't match fusionPolys data" ); }

   figWidth = (data$xRange[2] - data$xRange[1] + 1)/scale;
   figHeight = ((data$yRange[2] - data$yRange[1] + 1)/scale) * ratio;
   pdf( file=file, width=   figWidth + 1,
        height=  figHeight + 1
   );
   par( xpd=NA, mai= c(0.5, 0.5, 0.5, 0.5), omi= c(0,0,0,0) );

   plot( data$xRange, data$yRange, type="n", axes=F, xlab="", ylab="", main=title );

   for (i in 1:numGenes) {
      rect( geneRects[[i]]$xStarts, geneRects[[i]]$yBottoms, geneRects[[i]]$xEnds, geneRects[[i]]$yTops,
            col= geneRects[[i]]$fillColors, border= geneRects[[i]]$borderColors
      );
      if (! is.null(fusionPolys)) {
         polygon( fusionPolys[[i]]$x, fusionPolys[[i]]$y, col=fusionPolys[[i]]$fillColors, border=fusionPolys[[i]]$borderColors);
      }
      text( geneRects[[i]]$xStarts[1], min(geneRects[[i]]$yBottoms - scaleOffset/2), labels=geneNames[i], adj=c(0,1) );
   }
   lines( c(data$xRange[2] - 500, data$xRange[2]), c(yRange[1] + scaleOffset, yRange[1] + scaleOffset), lwd=4);
   text( data$xRange[2] - 250, yRange[1], labels='500 bp', adj=c(0.5,1))
   dev.off();

}

#' Plot a two-gene fusion expression plot.
#'
#' This generates one fusion expression plot file showing the expression of two
#' genes and associated (possibly multiple, possibly none) fusion points between
#' them. For plotting fusions to non-gene regions or just one gene, see
#' \code{\link{plotFusionExpressionOne}}. Essentially a wrapper around two
#' low-level functions: \code{\link{getFusionExpressionPlotData}} which
#' generates the plotting data given expression and fusion information, and
#' \code{\link{plotFusionExpressionPlotData}} which handles the actual plot
#' file creation.
#'
#' @param geneName1 The name of the first gene to plot.
#'
#' @param fusions1 A vector of genomic fusion ends within the first gene \code{geneName1}
#'
#' @param geneName2 The name of the second gene to plot.
#'
#' @param fusions2 A vector of genomic fusion ends within the second gene \code{geneName1}
#'
#' @param sample The sample name, for the plot and for lookup of expression in
#'   the \code{cohortExpressionDF}.
#'
#' @param cohortExpressionDF Provides the gene models for the specified gene and
#'   the normalized expression for the specified sample. See
#'   \code{\link{normExpressionData}} for format.
#'
#' @seealso \code{\link{normExpressionData}},
#'   \code{\link{plotFusionExpressionPlotData}},
#'   \code{\link{getFusionExpressionPlotData}},
#'   \code{\link{getCohortExonExpressionData}}
#'   \code{\link{plotFusionExpressionOne}}
#'
#' @return Nothing
#'
#' @export
plotFusionExpressionPair <- function( geneName1, geneName2, fusions1, fusions2, sample, cohortExpressionDF ) {

   gr1 <- grFromCohortDf(cohortExpressionDF, geneName1, sample, column='relativeExpression')[[geneName1]];
   gr2 <- grFromCohortDf(cohortExpressionDF, geneName2, sample, column='relativeExpression')[[geneName2]];

   gr1 <- expressionGrLabelMap( gr1 );
   gr2 <- expressionGrLabelMap( gr2 );

   geneNames <- c(geneName1, geneName2);
   genes <- list();
   genes[[geneNames[1]]] <- gr1
   genes[[geneNames[2]]] <- gr2
   fusions <- list()
   fusions[[geneNames[1]]] <- fusions1;
   fusions[[geneNames[2]]] <- fusions2;
   data <- getFusionExpressionPlotData( genes, fusions, geneNames );


   plotFusionExpressionPlotData( data, title=paste0( sample, " ", geneName1, "~", geneName2 ),
                                  file= paste0( sample, ".", geneName1, ".", geneName2, ".pdf" ) );
}

#' Plot a one-gene fusion expression plot.
#'
#' This generates one fusion expression plot file showing the expression of one
#' gene and associated (possibly multiple, possibly none) fusion points to
#' non-gene regions For plotting two genes, see
#' \code{\link{plotFusionExpressionPair}}. Essentially a wrapper around two
#' low-level functions: \code{\link{getFusionExpressionPlotData}} which
#' generates the plotting data given expression and fusion information, and
#' \code{\link{plotFusionExpressionPlotData}} which handles the actual plot
#' file creation.
#'
#' @param geneName1 The name of the gene being plotted.
#'
#' @param fusions1 A vector of genomic fusion ends within the \code{geneName1}
#'
#' @param sample The sample name, for the plot and for lookup of expression in
#'   the \code{cohortExpressionDF}.
#'
#' @param cohortExpressionDF Provides the gene models for the specified gene and
#'   the normalized expression for the specified sample. See
#'   \code{\link{normExpressionData}} for format.
#'
#' @seealso \code{\link{normExpressionData}},
#'   \code{\link{plotFusionExpressionPlotData}},
#'   \code{\link{getFusionExpressionPlotData}},
#'   \code{\link{getCohortExonExpressionData}}
#'   \code{\link{plotFusionExpressionPair}}
#'
#' @return Nothing
#'
#' @export
plotFusionExpressionOne <- function( geneName1, fusions1, sample, cohortExpressionDF ) {

   gr <- grFromCohortDf(cohortExpressionDF, geneName1, sample, column='relativeExpression')[[geneName1]];

   gr <- expressionGrLabelMap( gr );

   genes <- list();
   genes[[geneName1]] <- gr
   fusions <- NULL
   if (! is.null(fusions1)) {
      fusions <- list()
      fusions[[geneName1]] <- fusions1;
   }
   data <- getFusionExpressionPlotData( genes, fusions, geneName1 );


   plotFusionExpressionPlotData( data, title=paste0( sample, "   ", geneName1 ),
                                  file= paste0( sample, ".", geneName1, ".pdf" ) );
}


demo.plotFusionExpressionPlotData <- function() {
   data <- demo.getFusionExpressionPlotData();
   plotFusionExpressionPlotData( data, title="Testing: clown + TMPRSS2", file = 'clown_TMPRSS2.pdf' );
}

demo.plotFusionExpressionPlotData1 <- function() {
   data <- demo.getFusionExpressionPlotData1();
   plotFusionExpressionPlotData( data, title="TMPRSS2", file = 'TMPRSS2.pdf' );
}

demo.plotFusionExpressionPlotData4 <- function() {
   data <- demo.getFusionExpressionPlotData4();
   plotFusionExpressionPlotData( data, title="clownA TMPRSS2 clownB clownC", file = 'clownA_TMPRSS2_clownB_clownC.pdf' );
}

demo.plotFusionExpressionPlotData2 <- function() {
   data <- demo.getFusionExpressionPlotData2();
   plotFusionExpressionPlotData( data, title="Testing: TMPRSS2 + TMPRSS2", file = 'TMPRSS2_TMPRSS2.pdf' );
}


# plot.geneModelFusionOLD <- function(
#    gene,       color= NULL,  height= NULL,
#    gene2= NULL, color2= NULL, height2= NULL,
#    fusions= NULL,    fusions2= NULL,
#    fusionNames= NULL,   fusionNames2= NULL,
#    geneName= "GENE_1", geneName2= "GENE_2",
#    layout= c(left=0.01, gene=0.45, sep=0.08, gene2=0.45, right=0.01),
#    intronWidth= 100, scale= 'ABSOLUTE', exonBorder= 'black',
#    fusionOffset= 10, plot= TRUE, title=NULL, size= c(1000,100)
# ) {
#    # Plots one or two gene models (GenomicRanges objects) and annotates them
#    # with fusion points and gene names. Exons are rectangles whose width is
#    # proportional to the size of the exon. Exons are separated by white-
#    # space of fixed width corresponding to introns. This does not vary with
#    # intron size.
#    #
#    # Exons will be displayed with the color specified by the (wrapped) vector
#    # 'color'. If no color vector is provided, one will be read from the
#    # gene metadata column 'color'. If that is not found, the default is black.
#    #
#    # Exons will be displayed with the height specified by the (wrapped) vector
#    # 'height'. If no height vector is provided, one will be read from the
#    # gene metadata column 'height'. If that is not found, the default is 100.
#    #
#    #      A line underneath each gene will represent that gene's exon scale.
#    # Genes will by default be scaled to equal sizes. This fixes the visual width
#    # of a gene regardless of the size of the other displayed gene. Very long and
#    # very short genes will normally display as the same total width, but with
#    # different scale lines. The same gene will have the same visual width
#    # across different plots. To display both genes to the same scale, set
#    # scale= 'TRUE'RELATIVE'. This means both genes will have a visual width
#    # corresponding to their total length, but will have different visual widths
#    # in different plots. E.g. A short gene plotted with another short gene will
#    # be wider when compared to the same short gene paired with a long gene.
#    # More control over plot layout is available using the layout parameter as
#    # described below.
#    #
#    #     PARAM: gene: A GRanges object giving the exon model for the gene.
#    # Assumed to be non-overlapping, in-order, same strand, same chromosome
#    # Described as if on the positive strand (5' to 3' order). 'color' and
#    # 'height' metadata columns will be used as the color and height parameters.
#    # if not provided, defaulting to 'black' and 100 respectively.
#    #     PARAM: color= NULL: The color of the exon blocks.
#    # May be specified in gene as metadata column 'color', but over-ridden by
#    # this vector if specified. Wrapped if needed.
#    #     PARAM: height= NULL: The height of the exon blocks.
#    # May be specified in gene as metadata column 'height', but over-ridden by
#    # this vector if specified. Wrapped if needed.
#    #     PARAM: gene2= NULL: A second GRanges object, plotted after the first.
#    # If missing, associated arguments are ignored. 'color' and 'height' metadata
#    # columns can be used to provide the color2 and height2 parameters,
#    # defaulting to 'black' and 100 respectively.
#    #     PARAM: color2= NULL: as color but applied to gene2
#    #     PARAM: height2= NULL: as height but applied to gene2
#    #     PARAM: fusions: A vector of (+) genomic positions delimiting the fusion
#    # points to plot in gene 1. Points are displayed as colored triangles
#    # cycling through the colors 1:8.
#    #     PARAM: fusions2= NULL: As fusions, but applied to gene2
#    #     PARAM: fusionNames= NULL: Text centered above the fusions on gene.
#    # Often a name like "chr3:34,567-34,999". Ignored if no fusions specified.
#    # Note, names will overlap if positions are close.
#    #     PARAM: fusionNames2= NULL: As fusionNames, but applied to gene2
#    #     PARAM: geneName= NULL: Text below gene on the left side.
#    # Often a name like "TP53".
#    #     PARAM: geneName2= NULL: As geneName, but applied to gene2
#    #     PARAM: intronWidth= 100:  Absolute intron width.
#    # All introns are scaled to be this width before the whole gene is scaled
#    #     PARAM: exonBorder= 'black' The color of the exon box outline.
#    #    PARAM: scale= 'FIXED'|'RELATIVE'|'ABDSOLUTE': How to scale the plot of
#    # the two genes. If 'FIXED', plots the two genes at whatever scale is needed
#    # to fit into the layout windows specified at the size specified If
#    # 'RELATIVE' plots the genes so combined they fit into the combined gene +
#    # gene2 layout windows, but with relative window sizes based on the
#    # displayed length, not the layout. If 'ABSOLUTE' is given, displays as
#    # 'RELATIVE' but ignores size and uses one pixel per position for both genes,
#    # reverse calculating the size from the layout. 'FIXED' and 'RELATIVE' will
#    # look the same for one gene plots; 'RELATIVE' and 'ABSOLUTE' will look the
#    # same if zoomed, but will have different plot sizes at the same resolution.
#    # Case sensitive! NOTE: ABSOLUTE is not implemented!
#    #     PARAM: layout= c(left=0.01,gene=0.45,sep=0.08,gene2=0.45,right=0.01):
#    # Plots gene models into windows with these relative spaces. Doesn't have to
#    # add to one. If plotting only one gene, only left, gene, and right are used.
#    # If plotting two genes with scale='relative'RELATIVE', then the two genes
#    # will share gene + gene2, each taking the fraction based on its total exonic
#    # length. Note: probably should not have gene != gene2 when plotting with
#    # scale= 'FIXED' (the default)
#    #    PARAM: size= c(1000,100): sets the default total scaling for the plot,
#    # in pixels. This will determine the size when printed, depending on the
#    # pixel to inch setting of the printing device. size[1] (x width) is ignored
#    # if scale is 'ABSOLUTE'
#    #    PARAM: plot= TRUE: Draw plot and return data (list) invisibly.
#    # Set false not to plot and to autoprint returned data.
#    #    PARAM: title= NULL: Plot a title
#    #
#    #    RETURNS: List with all calculated data for plot.
#    #
#    #    TODO: Validate gene and gene2, if provided.
#    #    TODO: Consider y position of fusion end markers where height varies.
#    # Includes figuring out how to position fusionNames.
#    ###
#
#    # Guard
#    withGene2 = if ( is.null(gene2) ) FALSE else TRUE;
#
#    # Set up data structures for genes and the grapah, attach them to one list
#    # (data) for returning.
#    data  <- list();
#    graph <- list();
#    data$graph <- graph;
#    g1 <- list();
#    data$g1 <- g1;
#    if (withGene2) {
#       g2 <- list();
#       data$g2 <- g2;
#    }
#
#    #
#    # Validation
#    #
#
#    # PARAM scale validation
#    if ( ! scale %in% c('FIXED', 'RELATIVE', 'ABSOLUTE') ) {
#       stop("Scale must be one of 'FIXED', 'RELATIVE', or 'ABSOLUTE'.");
#    }
#
#    #
#    # PARAM layout validation
#    # Ensures that all named layout fields exist, are greater than 0 when used
#    # and are zero when not used.
#    #
#
#    layoutNames <- c('left', 'gene', 'sep', 'gene2', 'right');
#    if (any((layout < 0))) {
#       stop("Layout may not contain negative values.");
#    }
#
#    if (is.na(layout['gene']) || layout['gene'] == 0) {
#       stop("Layout must specify a non-zero gene fraction.");
#    }
#
#    # Assume if have gene2 and no layout fraction for gene2 that the gene
#    # layout fraction contains the sum of the layout for both genes, but
#    # must define a separator value
#    # Convienience for default "FIXED" scale.
#    if ( withGene2 ) {
#       if ( is.na( layout['sep'] ) || layout['sep'] == 0 ) {
#          stop("With gene2, layout must specify a non-zero sep fraction.");
#       }
#       if ( is.na( layout['gene2'] ) || layout['gene2'] == 0 ) {
#          layout['gene'] <- layout['gene']/2;
#          layout['gene2'] <- layout['gene']
#       }
#    }
#    else {
#       # no gene2 being plotted
#       layout['gene2'] <- 0
#       layout['sep'] <- 0
#    }
#
#    if (is.na(layout['left']))  { layout['left']  <- 0; }
#    if (is.na(layout['right'])) { layout['right'] <- 0; }
#
#    layout <- layout[layoutNames];
#    if (any(is.na(layout))) {
#       stop( paste0( "Bad Program. Layout should not have NA's: ", layout ));
#    }
#
#    #
#    # END: PARAM layout validation
#    #
#
#    #
#    # END: validation
#    #
#
#    #
#    # Transform layout to pixels, accounting for scale and size.
#    #
#
#    gene1width = sum( width(gene) ) + intronWidth * ( length( gaps(gene) ) - 1);
#    gene2width = 0;
#    if (withGene2) {
#       gene2width = sum( width(gene2) ) + intronWidth * ( length( gaps(gene2) ) - 1);
#    }
#
#    # Correct layout units where gene and gene 2 regions are relative to lengths
#    if (scale == 'RELATIVE' || scale == 'ABSOLUTE') {
#       # Keeping all 5 regions, but re-allocating the gene space for gene and
#       # gene2 based on the sum of the space allcated for gene and gene2
#       # split according to relative length
#
#       layoutUnit = sum( layout[c('gene', 'gene2')] ) / (gene1width + gene2width);
#       layout['gene']  = gene1width * layoutUnit;
#       layout['gene2'] = gene2width * layoutUnit;
#    }
#    graph$layout = layout;
#
#    # Convert from weighting to regions of whole graph (normalize so summ is 1)
#    regions = layout / sum(layout);
#    graph$regions = regions;
#
#    # Convert from fractional regions to pixel widths.
#    if (scale == 'ABSOLUTE') {
#       # Already have correct ration of regions for gene widths and normalized
#       # total width, so just need to figure out what the total size should be
#       # multiply by that, will get back sans fp error) the absolute sizes.
#       # 1 / fraction(gene + gene2) = new size[1] / absoulte_size(gene + gene2);
#       size[1] <- (gene1width + gene2width) / ( regions['gene'] + regions['gene2'] )
#    }
#    regionSizes = regions * size[1]
#    graph$regionSizes = regionSizes;
#
#    g1$xStarts <- genomicToModelCoordinates( start(gene), gene, intronWidth, regionSizes['gene'] ) + regionSizes['left'];
#    g1$xEnds   <- genomicToModelCoordinates( end(  gene), gene, intronWidth, regionSizes['gene'] ) + regionSizes['left'];
#    g1$yBottoms <- rep( 0, length(gene) );
#    g1$yTops    <- if( missing(height) ) elementMetadata(gene)[['height']] else rep( height, length(gene) );
#    g1$yMax     <- max(g1$yTops);
#    g1$colors <- if( missing(color) ) elementMetadata(gene)[['color']]  else rep(  color, length(gene) );
#    g1$xFusion <- if( ! is.null(fusions) ) genomicToModelCoordinates(fusions, gene, intronWidth, regionSizes['gene']) + regionSizes['left'] else NULL;
#    g1$yFusion <- if( ! is.null(fusions) ) rep( g1$yMax + fusionOffset, length(fusions) ) else NULL;
#
#    if (withGene2) {
#       g2 = list();
#       g2$offset = regionSizes['left'] + regionSizes['gene'] + regionSizes['sep']
#       g2$xStarts = genomicToModelCoordinates( start(gene2), gene2, intronWidth, regionSizes['gene2'] ) + g2$offset;
#       g2$xEnds   = genomicToModelCoordinates( end(  gene2), gene2, intronWidth, regionSizes['gene2'] ) + g2$offset;
#       g2$yBottoms = rep( 0,      length(gene2) );
#       g2$yTops  <- if( missing( height2 )) elementMetadata(gene2)[['height']] else rep( height2, length( gene2 ));
#       g2$yMax     <- max(g2$yTops);
#       g2$colors <- if( missing(  color2 )) elementMetadata(gene2)[['color']]  else rep(  color2, length( gene2 ));
#       g2$xFusion <- if( ! is.null( fusions2 )) genomicToModelCoordinates( fusions2, gene2, intronWidth, regionSizes['gene2'] ) + g2$offset else NULL;
#       g2$yFusion <- if( ! is.null( fusions2 )) rep(g2$yMax + fusionOffset, length(fusions2)) else NULL;
#    }
#
#    graph$xRange = c(0, size[1]);
#    graph$yRange = c(0, g1$yMax + 1);
#    if (! is.null(geneName) ) {
#       graph$gene1NamePos = c( regionSizes['left'], min(g1$yBottoms) );
#    }
#    if (! is.null(geneName2)  && withGene2) {
#       graph$gene2NamePos = c( g2$offset, min(g2$yBottoms) );
#    }
#
#    if (plot) {
#       # Allow plotting outside graph (in margins)
#       oldClipping <- par()$xpd;
#
#       scale=1000;
#       ratio=4
#       par(xpd=NA);
#       pdf( file="temp.pdf", width=graph$xRange[2]/scale + 2.2, height=ratio*(graph$yRange[2]/scale) + 2.2);
# #     pdf( file="temp.pdf");
#
#
#       # Set up plotting page and plot coordinates
#       plot(graph$xRange, graph$yRange, type="n", axes=F, xlab="", ylab="", mar=c(0.5, 0, 0.5, 0) + 0.1);
#       # Add exon boxes
#       if (withGene2) {
#          rect(
#             c(g1$xStarts, g2$xStarts), c(g1$yBottoms, g2$yBottoms),
#             c(g1$xEnds,   g2$xEnds),   c(g1$yTops,    g2$yTops),
#             col = c(g1$colors, g2$colors), border= exonBorder
#          );
#       }
#       else {
#          rect(
#             g1$xStarts, g1$yBottoms,
#             g1$xEnds,   g1$yTops,
#             col = g1$colors, border= exonBorder
#          );
#       }
#
#       # Add fusions
#       if (! is.null(fusions)) {
#          plot.polyPoints( g1$xFusion, g1$yFusion, shape='triangleDown', col=1:8, ratio=0.1);
#          if (! is.null(fusionNames) ) {
#             text( g1$xFusion, g1$yFusion, fusionNames, pos=3)
#          }
#       }
#       if (! is.null(fusions2)  && withGene2) {
#          plot.polyPoints( g2$xFusion, g2$yFusion, shape='triangleUp', col=1:8, ratio=0.1);
#          if (! is.null(fusionNames2) ) {
#             text( g2$xFusion, g2$yFusion, fusionNames2, pos=3)
#          }
#       }
#
#       # Add gene names
#       if (! is.null(geneName)) {
#          strand = paste0( "(", strand(gene)[1], ")" );
#          text( graph$gene1NamePos[1], graph$gene1NamePos[2], paste(geneName, strand), pos=1);
#       }
#       if (! is.null(geneName2)  && withGene2) {
#          strand = paste0( "(", strand(gene2)[1], ")" );
#          text( graph$gene2NamePos[1], graph$gene2NamePos[2], paste(geneName2, strand), pos=1);
#       }
#
#       if (! is.null(title)) {
#          title(title);
#       }
#       # Restore par
#
# #      dev.off();
#       par(xpd=oldClipping);
#
#       return( invisible(data) );
#    }
#    else {
#       return( data );
#    }
# }

# plot.fusionExpressionOLD <- function (
#    fusionDF, expresionDF,
#    sample, gene= NULL, gene2= NULL,
#    is.normalized = TRUE, is.skipNoFusions= TRUE,
#    is.pairedFusionEnds= TRUE, is.orderedFusionEnds= FALSE,
#    baseDir = "./", plotType= "pdf", borderCol= 'black'
# ) {
#    ###
#    #     Plots one or more gene expression plots, returning a dataframe
#    # describing the plots, including the filename. By default plots will be
#    # generated into the current working directory as pdfs. See the baseDir
#    # and plotType options to change defaults.
#    #     Expression data and gene model information is input as the expressionDF
#    # data frame, and fusion points are input as the fusionDF data frame. These
#    # data frames can be cohort-wide, or limited only to what is to be ploted.
#    # It is assumed the expression data has been cohort-wide normalized. Which
#    # sample and genes to include in a plot are specified by the (wrapped) input
#    # vectors sample, gene, and gene2. One plot is generated for each
#    # supplied triplet as long as there is at least one fusion displayed. See the
#    # is.normalized to change defaults.
#    #     Each plot describes the relative exon expression pattern for two genes
#    # in a sample and the fusion points within those genes in that sample. By
#    # default only fusions with an end in each gene are shown, but the order in
#    # which the fusion ends are specified in the fusionDF is ignored See the
#    # is.pairedFusionEnds and is.orderedFusionEnds options to change defaults.
#    ###
#    #     PARAM: fusionDF - The fusion positions data frame. Fusion end order is
#    # ignored unless is.orderedFusionEnds is set.
#    #     PARAM: expressionDF - exon expresion data frame, normalized unless
#    # is.normalized is set false
#    #     PARAM: sample - name of sample to plot
#    #     PARAM: gene - name of first gene to plot
#    #     PARAM: gene2 - name of second gene to plot
#    #     PARAM: is.normalized= TRUE - It is assumed the provided expression data
#    # is normalized. If not, set this to FALSE and the expression data frame will
#    # be normalized prior to processing. It is probabbly better to do this before
#    # providing the data, as this is a long-running operation for big data frames
#    # and is best saved for reuse. If using this option to normalize, the
#    # provided expression data frame MUST contain expression data for EVERY
#    # sample in the cohort (for the genes being ploted), including data for these
#    # genes from samples without fusions.
#    #     PARAM: borderCol - The color of the exon borders.
#    #     PARAM: is.skipNoFusions= TRUE - If neither gene plotted for a sample
#    # has a fusion, the plot will not be generated unliss this is set FALSE.
#    # If any fusion is plotted in either gene, the plot will always be generated.
#    #     PARAM: is.pairedFusionEnds= TRUE - Only fusions with one end in each
#    # gene will be included in a plot. If fusions to different genes exist in
#    # the same sample, these fusion points will be plotted only if this is set
#    # FALSE. TODO: This is unimplemented.
#    #     PARAM: is.orderedFusionEnds= FALSE - Genes are plotted in order as
#    # gene and then gene2. Fusions with an end in both genes are normally
#    # included regardless of being specifeid in the fusion data frame as gene -
#    # gene2 or gene2 - gene. To only show fusions in the specified order, set
#    # this to TRUE.  TODO: This is unimplemented.
#    #     PARAM: baseDir= "./" - The directory to generate plots into
#    #     PARAM: plotType= "pdf" - The type of plot to generate, currently only
#    # pdf is supported.
#    ###
#    #     RETURNS: df - Dataframe describing the plots made or considered based
#    # on input data. If expression data exists for the given sample and both
#    # genes, a plot will be considered. Various options will determine if a plot
#    # is made. If a plot is not made, the filename will be NA. If is.onlyPlotted
#    # is set, the returned data frame will contain only those records for which
#    # plots were made (those where filename is not NA).
#    #     df$sample - The name of the sample considered for plotting.
#    #     df$gene - The name of the gene to plot on the left.
#    #     df$gene2 - The name of the gene to plot on the right.
#    #     df$file - The name of the plot file, NA if this was not plotted.
#    #     df$baseDir - The base dir, included to allow access by absolute path.
#    ###
#
#    # Wrap samples, gene and gene2 lists to get list of plots. Makes it easy to
#    # plot a pair of genes for many samples, or many genes for 1 sample, or
#    # several genes for several samples. Probably not good if have same gene name
#    # in both gene and gene2 lists, unless only have one gene in one of those
#    # lists.
#    plotDF <- data.frame(sample,gene,gene2,file=c(NA_character_),baseDir,stringsAsFactors=FALSE);
#    plotNames <- character(0);
#
#    # Optimization for large fusion candidate data frames and many plots. Little
#    # additional overhead if small list.
#    fusionCandidatesDF <- getCandidateFusions(fusionDF, sample, gene, gene2);
#
#    # Plot each line in data frame and
#    for (row in 1:nrow(plotDF)) {
#       plot = list();
#       plot$sample <- plotDF$sample[row];
#       plot$gene1 <- plotDF$gene[row];
#       plot$gene2 <- plotDF$gene2[row];
#       fdf <- getCandidateFusions(fusionCandidatesDF, plot$sample, plot$gene1, plot$gene2);
#       if (is.pairedFusionEnds && is.orderedFusionEnds) {
#          selectVec <- fdf$gene1 == plot$gene1 & fdf$gene2 == plot$gene2;
#       }
#       else if (is.pairedFusionEnds && ! is.orderedFusionEnds) {
#          selectVec <- (fdf$gene1 == plot$gene1 & fdf$gene2 == plot$gene2) |
#                          (fdf$gene1 == plot$gene2 & fdf$gene2 == plot$gene1)
#       }
#       else if (! is.pairedFusionEnds && is.orderedFusionEnds) {
#          selectVec <- fdf$gene1 == plot$gene1 | fdf$gene2 = plot$gene2;
#       }
#       else if (! is.pairedFusionEnds && ! is.orderedFusionEnds) {
#          selectVec <- fdf$gene1 == plot$gene1 | fdf$gene1 == plot$gene2 |
#                         fdf$gene2 == plot$gene1 | fdf$gene2 == plot$gene2;
#       }
#       else {
#          stop("Bad program! plot.fusionExpression() logic not correct!");
#       }
#       fdf <- fdf[selectVec,];
#       selectVec <- fdf$gene1 == plot$gene1 | fdf$gene2 == plot$gene2;
#       plot$fusions1 <- fdf[ selectVec, c('position1')];
#       plot$fusions2 <- fdf[ selectVec, c('position2')];
#       if ( length(plot$fusions1) != 0 || length(plot$fusions2) != 0 ) {
#          plot$name = paste0(plot$sample, '.', plot$gene1, '.', plot$gene2, '.', plotType );
#          plotDF$file[row] <- plot$name;
#
#          expressionGr1 <- expressionGrFactory( expresionDF, gene= plot$gene1, sample= plot$sample );
#          expressionGr2 <- expressionGrFactory( expresionDF, gene= plot$gene2, sample= plot$sample );
#          cutPoints <- c(Inf, 0.3, 0.1, -0.1, -0.3, -Inf);
#          expressionGr1 <- expressionGrLabelMap(expressionGr1);
#          expressionGr2 <- expressionGrLabelMap(expressionGr2);
#
#          pdf(file=paste0(baseDir,'/',plot$name), width=1000/96, height=(100+220)/96 );
#          plot.geneModelFusion(
#             gene=expressionGr1, height= 100, gene2= expressionGr2, height2= 100,
#             fusions= plot$fusions1, fusions2= plot$fusions2,
#             geneName= plot$gene1, geneName2= plot$gene2,
#             title=plot$name, size= c(1000,100), exonBorder= borderCol
#          );
#          dev.off();
#
#       }
#
#    }
#    return(plotDF);
# }


demo.makeGrForTmprss2 <- function() {
   exons <-
      "42836479-42838080,42839661-42839813,42840323-42840465,42842575-42842670,42843733-42843908,42845252-42845423,42848504-42848547,42851099-42851209,42852403-42852529,42860321-42860440,42861434-42861520,42866283-42866505,42870046-42870116,42879877-42879992";
   iRanges <- IRanges(start = as.numeric(sub("-.+", "", unlist(
      strsplit(exons, ",")
   ))),
   end =   as.numeric(sub(".+-", "", unlist(
      strsplit(exons, ",")
   ))));
   gr <- GRanges(seqnames <- as.factor("chr21"),
                 ranges   <- iRanges,
                 strand   <- '-');
   return(gr)
}

demo.makeGrForTest <- function() {
   gr <- GRanges(
      ranges= IRanges( start= c( 101, 301, 511, 513, 813 ), end= c( 200, 500, 511, 612, 822 )),
      seqnames= c( "chr1" ),
      strand= c( "+" )
   );
   gr <- grAddColumn(gr, 'exonFillColor', c('red', 'orange', 'yellow', 'green', 'blue'));
   gr <- grAddColumn(gr, 'exonBorderColor', c('blue', 'green', 'black', 'orange', 'red'));
   gr <- grAddColumn(gr, 'exonHeight', c(50, 100, 200, 100, 50));
   return(gr);
}


# runGenePlotDemo <- function(doPlot=TRUE) {
#    TMPRSS2$GRanges <- demo.makeGrForTmprss2();
#    colors = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020");
#    binEnds = c(0,0.25,0.5,1,2,Inf);
#    TMPRSS2$GRanges <- addColumn(
#       TMPRSS2$GRanges, 'color', sortDataIntoBins(
#          c(0.01,0.125,0.25,0.25,0.5,0.5,1,1,2,2,4,4,8,16), binEnds, colors
#       )
#    );
#    # 108370   uc002yxa.2	chr21:39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39772496-39772567,39774479-39774559,39775428-39775631,39795332-39795483,39817327-39817544,39870287-39870428:-	ERG|2078	CDSstart=125;CDSstop=1564
#    # 108372	uc011aek.1	chr21:39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39772496-39772567,39774479-39774559,39775428-39775631,39795332-39795483,39870287-39870428:-	ERG|2078	CDSstart=183;CDSstop=1346
#    # 108374	uc010gnv.2	chr21:39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39774479-39774559,39775428-39775631,39795332-39795483,39947586-39947671,39956768-39956869:-	ERG|2078	CDSstart=229;CDSstop=1320
#    # 108376	uc010gnw.2	chr21:39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39772496-39772567,39774479-39774559,39775428-39775631,39795332-39795483,39817327-39817544,39947586-39947671,39956768-39956869,40032446-40032591:-	ERG|2078	CDSstart=296;CDSstop=1756
#    # 108378	uc010gnx.2	chr21:39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39774479-39774559,39775428-39775631,39795332-39795483,39817327-39817544,39947586-39947671,39956768-39956869,40032446-40032591:-	ERG|2078	CDSstart=296;CDSstop=1684
#    # 108380	uc011ael.1	chr21:39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39772496-39772567,39774479-39774559,39775428-39775631,39795332-39795483,39817327-39817544,39947586-39947671,39956768-39956869,40033582-40033704:-	ERG|2078	CDSstart=273;CDSstop=1733
#    # 108382	uc002yxb.2	chr21:39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39774479-39774559,39775428-39775631,39795332-39795483,39817327-39817544,39947586-39947671,39956768-39956869,40033582-40033704:-	ERG|2078	CDSstart=273;CDSstop=1661
#    # 108384	uc011aem.1	chr21:39763872-39764366,39795332-39795483,39817327-39817544,39870287-39870428:-	ERG|2078	CDSstart=125;CDSstop=640
#    # 108388	uc002yxc.3	chr21:39771999-39772567,39774479-39774559,39775428-39775631,39795332-39795483,39817327-39817544,39947586-39947671,39956768-39956869,40033582-40033704:-	ERG|2078	CDSstart=273;CDSstop=1250
#    ERG = list();
#    ERG$chr <- as.factor( "chr21" );
#    ERG$exons = "39751952-39755845,39762917-39762964,39763581-39763637,39764298-39764366,39774479-39774559,39775428-39775631,39795332-39795483,39817327-39817544,39947586-39947671,39956768-39956869,40033582-40033704"
#    ERG$iranges <- IRanges(
#       start= as.numeric(sub("-.+", "", unlist(strsplit(ERG$exons, ",")))),
#       end=   as.numeric(sub(".+-", "", unlist(strsplit(ERG$exons, ","))))
#    );
#    ERG$strand = "-"
#    ERG$GRanges <- GRanges(
#       seqnames= ERG$chr,
#       ranges=ERG$iranges,
#       strand=ERG$strand
#    );
#    ERG$GRanges <- addColumn(
#       ERG$GRanges, 'color', sortDataIntoBins(
#          c(18,12,6,3,1,0.9,0.5,0.3,0.2,0.1,0.05), binEnds, colors
#       )
#    );
#
#    gafFile = "./data/TCGA.hg19.June2011.gaf.geneModels.txt"
#    pradExonFileFile = "./data/pradExonQuantFiles.txt"
#    files <- c(
#       "/datastore/rclbg/nextgenout3/seqware-analysis/illumina/120425_UNC11-SN627_0224_BD0VAKACXX/seqware-0.7.0_Mapsplice-0.7.4/120425_UNC11-SN627_0224_BD0VAKACXX_5_GATCAG/exon_quantification.txt",
#       "/datastore/rclbg/nextgenout3/seqware-analysis/illumina/120425_UNC11-SN627_0224_BD0VAKACXX/seqware-0.7.0_Mapsplice-0.7.4/120425_UNC11-SN627_0224_BD0VAKACXX_6_GATCAG/exon_quantification.txt",
#       "/datastore/rclbg/nextgenout3/seqware-analysis/illumina/120425_UNC11-SN627_0224_BD0VAKACXX/seqware-0.7.0_Mapsplice-0.7.4/120425_UNC11-SN627_0224_BD0VAKACXX_7_TAGCTT/exon_quantification.txt"
#    );
#
#
#    plot.geneModelFusion(gene2=ERG$GRanges, geneName2="ERG", height2=100,
#              gene=TMPRSS2$GRanges, geneName="TMPRSS2", height=100,
#              fusions2= c(39882952, 39752952, 39754952), fusionNames2= c("39,782,952","39,752,952","39,754,952"),
#              fusions= c(42852530), fusionNames= "42,852,530",
#              plot=TRUE, title="Its a graph!", fusionOffset=3, scale= 'FIXED'
#    );
#
# }
