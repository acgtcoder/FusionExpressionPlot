# library(GenomicRanges);
# library(RColorBrewer);

#' @import IRanges
#' @import GenomicRanges
#' @import RColorBrewer
NULL

DEBUG <- 1;

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
# Data reduction functions: Continuous to category
###

# Get bin name for the range (bin) data points fall into. Good for colors.
# Note, binEnds will be sorted before assigning labels.
dataToBinLabels <- function( data, binEnds, binLabels) {
   ###
   #     Bins the data as specified by binEnds, then returns a vector of the
   # same length as data, but with each value replaced by the binLabel. If the
   # bin labels are colors, this converts data into a color scale. If the bin
   # labels are the upper (or lower) bin ends, in converts values into their
   # category maxima (minima).
   #
   #     PARAM: data - A vector of numeric data to bin.
   #     PARAM: binEnds - The end of the bins into which the data is sorted.
   # One bin would be c(-Inf, Inf). Two bins equally weighted would be
   # c(-Inf, median(data), Inf). This vector will be sorted before assigning
   # labels, so reverse the label order instead of providing bins as high to low.
   #     PARAM: binLabels - The vector of values corresponding to the bins to
   # convert the data into. For the two bin example above, might use
   # c('low','high') or c('blue','red');
   bins = cut(data, binEnds, labels = 1:(length(binEnds) - 1));
   return(binLabels[as.numeric(levels(bins)[bins])]);
}

###
# GRanges object functions
###

# Add a data column to a GRanges object.
grAddColumn <- function( gr, name, vec ) {
   elementMetadata(gr)[[name]] <- vec;
   return(gr);
}

# Add a column of color values to a GRange object based on a numeric data column.
# Must supply the ranges and the colors to go with those ranges. Usess colorBrewer
grColorMap <- function(
   gr, binEnds, column=1, brewerPaletteName= 'RdBu',
   colors= rev(brewer.pal( length(binEnds) - 1, brewerPaletteName)),
   as.column = 'exonColor', annotate=TRUE
) {
   ###
   #      Add a column of color values as metadata to a GRanges object, mapped
   # from an existing column of numeric metadata. Which column to map can be
   # specified by column name or index. The data from that column will be binned
   # based on a specified sequence of binEnds.
   #
   #      The color assigned will be those from the brewer palette named (with
   # the number of levels auto selected based on the bin ends). It is possible
   # to assigne a vector of colors directly, in which case the brewer palette
   # name is ignored. Note that the order of the bin ends matters, as probably
   # want high to low if representing as red-blue, but may want low to high if
   # using a sequential pallete.
   #
   #     Returns the gRanges object with the color metadata column added, or if
   # annotate is set false, returns just the color vector.
   ###
   #     PARAM gr - The GRanges object to annotate.
   #     PARAM column= 1 - The column index or name of the numeric data in
   # gr to color map.
   #     PARAM binEnds= The sequence of cut-points defining the bins into which
   # the data will be split. Each bin will correspond to a color in the provided
   # color palette, providing a mapping of each value to a color. Bins are open
   # on the left and closed on the right (meaning values exactly matching bin
   # ends will be in the bin to the left.) No data value should be less than or
   # equal to the lowest bin end, or greater than the upper bin end, otherwise
   # NA values will result. Normally one would have the first and last bin end
   # as the special values Inf and -Inf.
   #     PARAM brewerPaletteName= 'RdBu' - The name of the RColorBrewer palatte.
   # For example palattes: library(RColorBrewer); display.brewer.all(). Setting
   # colors causes this to be ignored.
   #     For diverging values, can have 3 - 11 bins, allowed values are:
   # BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral.
   #     For sequential values, can have 3 - 8 value, allowed values are:
   # Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples
   # RdPu Reds YlGn YlGnBu YlOrBr YlOrRd.
   #     For qualitative palletes, can have 3 to (n) where n varies by palette:
   # Accent (8), Dark2 (8), Paired (12), Pastel1 (9), Pastel2 (8), Set1 (9),
   # Set2 (8), Set3 (12)
   #     PARAM colors= * - Can specify an explicit vector of colors to
   # override the color-brewer palatte selection. The number of colors will
   # have to be one less than the number of bin ends.
   #     PARAM as.column= 'color' - The name of the metadata column added to
   # the GRanges object.
   ###
   #     RETURNS gr | vector - The GRanges object passed in, but with the color
   # columnn added, or, if annotate is set FALSE, returns just the color vector.
   ###

   colorCol <- dataToBinLabels(elementMetadata(gr)[[column]],binEnds,colors);
   if (annotate) {
      return(grAddColumn(gr, as.column, colorCol));
   }
   else {
      return(colorCol);
   }
}

###
# GRanges object factories
###

# Get a list of GRange objects by gene for a column of sample data
# extracted from a data frame with gene + exon rows and sample columns.

#' @export
grFromCohortDf <- function( df, genes, sample, column='exonData' ) {
   ###
   #     Given a dataframe describing exon expression (for all samples in the
   # cohort to consider), a vectr df gene names, and a sample name, returns a
   # list of GRanges object, by gene, with the sample data column as an extra
   # data column.
   #     The input dataframe must have the following columns: chr, start, end,
   # strand, gene, and then one column for each sample. The data in the sample
   # columns are assumed to be the same measurement for each exon in every gene.
   # Start and end are assumed to be 1 based, but this does not really matter
   # as long as all positions here and in any position data (i.e. fusions)
   # used with this data all have the same base.
   ###
   #     PARAM: df - A dataframe with gene model and exon expression.
   #          df$chr    - The chromosome this exon is on, as 'chr1'...'chr22',
   # 'chrX', 'chrY', or 'chrM'.
   #          df$start  - The position of the first base of the exon.
   #          df$end    - The position of the last base of the exon.
   #          df$strand - The strand of this exon (+,-, or *). All exons in a
   # gene should have the same strand.
   #          df$gene   - The name of this gene. All exons in a gene should have
   # the same name, and no exons from different genes should have the same name.
   #          df$***    - The normalized expression value for this exon in sample
   # named "***". Should have one column for every sample in the cohort. Sample
   # names may not be "chr", "start", "end", "strand", "gene", or be repeated.
   #      PARAM: genes - The gene names to extract from the df. Must exist.
   #      PARAM: sample - The sample to extract from the df. Must exist.
   #      PARAM: column - The metadata column name under which the expression
   # data is added to the granges object.
   ###
   #     RETURN: A List of GRanges object for the gene, with start and end positions of
   # the exon specified, and with a metadata column "relativeExpression"
   # giving the relative expression for the selected gene in the selected sample.
   ###
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

# The grColorMap as used with exon expression
expressionGrColorMap <- function( gr ) {
   return( grColorMap( gr, column='relativeExpression', as.column = 'exonFillColor',
                       binEnds= c(Inf,0.3,0.1,-0.1,-0.3,-Inf) ))
}

#
# Data for plotting objects
#

# Convert genomic positions to plot positions based on a gene model
# with fixed exon size.
genomicToModelCoordinates <- function(
   pos, gene, gapWidth= 100, outOfBounds='warning'
) {
   ###
   #      Converts pos, a vector of (+ strand) genomic positions to a vector of
   # coordinates in a model as described by the GRanges object gene. The
   # coordinates are 1-1 with those in the gene elements, but with gaps
   # scaled to a fixed gapWidth. The GRanges gene must be all on one strand of
   # one chromosome. '-' strand gene models are reveresed so all coordinates are
   # 5' - 3' transformed (based on first element in granges object).
   #      The first genomic start position for any element in the GRanges object
   # is model coordinate '1'. The last model coordinate is the sum of all the
   # element widths plus the number of gaps * the gapthWidth. Genomic positions
   # before the model (less than any element start) are converted to model
   # coordinate 0; genomic positions after the mode (greater than any element
   # end are converted to the last model coordinate + 1. This can be controlled
   # by the outOfBounds parameter (Not implemented).
   #
   #      Due to the required scaling and rounding, there is no way to exactly
   # define the model coordinate of genomic positions in gaps such that both
   # gaps larger than gapWidth and gaps smaller than gapWidth are mapped to
   # intuitively and identically in both positive and negative strand genes.
   # The chosen method is as follows:
   #
   #     If a genomic position falls within a gap, it is converted to a
   # gap-relative position, then shifted 1/2 a position in the 5' direction.
   # This shifted, relative position [0.5 ... thisGapWidth - 0.5] is then
   # scaled by the gapthWidth (gapWidth/thisGapWidth). The next largest integer
   # (ceiling) is taken as the offset into the gap in the model, with the
   # returned model coordinate being "the sum of the widths of the exons before"
   #  + (gapWidth * "the number of gaps before") + "this offset". If the strand
   # is negative, the coordinate is relative to the last exon base of the mode.
   # More formally:
   #
   #    coordinate = sum( widths(elementsBefore)) + count(gapsBefore) * gapWidth
   #       + ceiling( (p - endPos(elementJustBefore) - 0.5) *
   #                   (gapWidth/thisGapWidth))
   #      if (strand == '--') { coordinate = sum( widths( allElements ))
   #                                       + count(allGapsB) * gapWidth + 1
   #                                       - coordinate
   ###
   #     PARAM: pos    (REQ)  IntVec
   # A vector of chromosome genomic positions to convert to model coordinates.
   # Assumed to be on the same chromosome as the provided gene model, although
   # not required to be within the gene model.
   #      PARAM: gene  (REQ)  GRanges
   # Object giving the model, i.e. gene exons.
   # All elements must be on the same chromosome on the same strand in sequence
   # order with no overlaps and no 0 length introns. Assumes '-' strand genes
   # are described as if their shadow on the '+' strand was the real gene, with
   # ony the '-' indicateing anything aobut a reverse order.
   #      PARAM: gapWidth=100  INT
   # Fixed gap width. Equivalent is size to an exon of the same length.
   #      PARAM: outOfBounds='warning'
   # By default, genomic positions can be supplied that are outside the model,
   # but a warning will be given, once. "Position outside the gene model found: #".
   # Setting outOfBounds='error' will cause the first such position found to
   # abort. Setting outOfBounds='ok' means no warning or error will be generated.
   #      RETURNS: IntVec
   # The genomic positions converted to model coordinates.
   #     TODO: Valdate gene
   # Check that gene is GRanges object, single chromosome, single stranded,
   # reduced, and ordered.
   ###

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

test_genomicToModelCoordinates <- function() {

   # Exons: 101-200,301-500,511-511,513-612,813-822
   #          100    200       1      100     10
   # Introns:   201-300,501-510,512-512,613-812
   #              100     10       1      200
   geneModel <- demo.makeGrForTest();

   # Before
   got = genomicToModelCoordinates( c(1, 100, 101), geneModel, 100 );
   want = c(0, 0, 1);
   stopifnot( got == want );
   # First exon
   got = genomicToModelCoordinates( c(101, 150, 200), geneModel, 100 );
   want = c(1, 50, 100);
   stopifnot( got == want );
   # First intron
   got = genomicToModelCoordinates( c(201, 250, 300), geneModel, 100 );
   want = c(101, 150, 200);
   stopifnot( got == want );
   got = genomicToModelCoordinates( c(301, 400, 500), geneModel, 100 );
   want = c(201, 300, 400);
   stopifnot( got == want );
   got = genomicToModelCoordinates( c(501, 502, 510), geneModel, 100 );
   want = c(405, 415, 495);
   stopifnot( got == want );
   got = genomicToModelCoordinates( c(511, 512, 513), geneModel, 100 );
   want = c(501, 551, 602);
   stopifnot( got == want );
   got = genomicToModelCoordinates( c(613, 712, 812), geneModel, 100 );
   want = c(702, 751, 801);
   stopifnot( got == want );
   # After
   got = genomicToModelCoordinates( c(822, 823, 10000), geneModel, 100 );
   want = c(811, 812, 812);
   stopifnot( got == want );

   geneModel <- GRanges(
      ranges= IRanges( start= c( 101, 301, 511, 513, 813 ), end= c( 200, 500, 511, 612, 822 )),
      seqnames= c( "chr1" ),
      strand= c( "-" )
   );

   # Before
   got = genomicToModelCoordinates( c(822, 823, 10000), geneModel, 100 );
   want = c(1, 0, 0);
   stopifnot( got == want );
   # First exon
   got = genomicToModelCoordinates( c(822, 818, 813), geneModel, 100 );
   want = c(1, 5, 10);
   stopifnot( got == want );

   return(TRUE);
}

# Given a gene model, return a data object that can be used to draw a gene-exon
# plot or rectangles
grToRect <- function (
   gene, intronWidth= 100, exonHeights= 100,
   exonFillColors= 'blue', exonBorderColors= 'black'
) {
   # Given a genomicRanges Object, return the data needed to plot it as a
   # sequence of rectangles based on the exon widths. Introns are all fixed at
   # the single intronWidth. Height, fill color, and border color for the exon
   # rectangles can be specified on a per-exon basis (vector-wrapped), or given in
   # the genomicRanges object (as metadata colums" exonHeight, exonFillColor,
   # and exonBorderColor). Returns a list of data as specified below. Uses the
   # genomicToModelCoordinates function, passing throught the parameters
   # needed
   ###
   #      PARAM: gene  (REQ)  GRanges
   # Object giving the gene to calculate the parameters for. Passed through to
   # genomicToModelCoordinates as its gene parameter. Besides the criteria
   # specified there, may also contain data columns 'exonHeight', 'exonFillColor'
   # and 'exonBorderColor'. If present, the values in these columns will be used
   # instead on the equivalently names parameters to this function.
   #     PARAM: intronWidth=100  Int
   # The width of the itron, the x spacing between the exon rectangles. Passed
   # through to genomicToModelCoordinates as its gapWidth parameter.
   ###

   rect <- list();
   rect$onReverseStrand = (as.character(strand(gene)[1]) == '-');

   if (rect$onReverseStrand) {
      rect$xStarts  <- rev(genomicToModelCoordinates( end(  gene), gene, intronWidth ));
      rect$xEnds    <- rev(genomicToModelCoordinates( start(gene), gene, intronWidth ));
      rect$yBottoms <- rep( 0, length( gene ));
      rect$yTops        <- rev(elementMetadata(gene)[['exonHeight']]);
      rect$fillColors   <- rev(elementMetadata(gene)[['exonFillColor']]);
      rect$borderColors <- rev(elementMetadata(gene)[['exonBorderColor']]);
   }
   else {
      rect$xStarts  <- genomicToModelCoordinates( start(gene), gene, intronWidth );
      rect$xEnds    <- genomicToModelCoordinates( end(  gene), gene, intronWidth );
      rect$yBottoms <- rep( 0, length( gene ));
      rect$yTops        <- elementMetadata(gene)[['exonHeight']];
      rect$fillColors   <- elementMetadata(gene)[['exonFillColor']];
      rect$borderColors <- elementMetadata(gene)[['exonBorderColor']];
   }

   if ( is.null( rect$yTops        )) { rect$yTops        <- rep_len( exonHeights,      length( gene )); }
   if ( is.null( rect$fillColors   )) { rect$fillColors   <- rep_len( exonFillColors,   length( gene )); }
   if ( is.null( rect$borderColors )) { rect$borderColors <- rep_len( exonBorderColors, length( gene )); }

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

# Convert a vector of coordinate points to a data object that can be drawn
# with polygon to plot a shape at each of those coordinates.
pointToPoly <- function (
   x, y, widths=20, heights=20, shapes= 'triangleDown', sides= 5,
   fillColors='red', borderColors='black',
   scale=c(4,1), resScale=c(1,1)
) {

   ###
   #     Calculate the values need to plot a sequence of points as polygons
   # specified by shape. Currently only triangles pointing up and triangles
   # pointing down can be specified. Each point can be described separately by
   # passing in a vector for the defining parameter. Other than x and y
   # coordinates of the points to plot and the two scale paramters, provided
   # values will be wrapped as needed. The returned x and y coordinates are
   # in block separated by NA values, which will be ploted by the "polygon"
   # command as separated polygons.
   #     The scale and resScale values affect only the shape and size of the
   # polygons. It does not affect the base positions drawn at. i.e. They allows
   # for symetrical triangles to be drawn symetrically on graphs with different
   # x and y scales, and with different device resolutions. Conceptually this
   # could be handled by pre-scaling the widths and heights, but it is easier to
   # separate this out.
   #     PARAM: x - The x position in plot coordinates of the poly-points
   #     PARAM: y - The y position in plot coordinates of the poly-points
   #     PARAM: widths= 5 - The widths in plot coordinates of the poly-points
   #     PARAM: heights= 5 - The height in plot coordinates of the poly-points
   #     PARAM: sides= 5 - The position of the (larger) poly over the point. Can
   # have a value from 1-9, oriented like a touch pad, 1 has uper left corner
   # of poly as the point, 5 has the middle over the point, and 9 has the
   # lower right as the point.
   #     PARAM: shapes= 'triangleDown' - Shape of 'polygons'point'. Allowed
   # shapes are triangleDown and triangleUp.
   #     PARAM: fillColors= 'red' - Color filling the 'point'.
   #     PARAM: borderColors= 'black' - Color of the 'point' borders.
   #     PARAM: scale= c(1,1) - The scale factor of the x and y axis of the
   # plot, i.e if plot is 100x500, scale should be c(1,5) or c(0.2,1)
   #     PARAM: resScale= c(1,1) -The scale of the device, i.e. if the screen
   # resolution is 1024x768, this should be c(1,1024/768) or c(768/1024,1).
   ###
   #     RETURNS: List with data needed to plot points as polygons, with
   #          $count = number of points plotted as polygons.
   #          $xRadii = 1/2 widths of points, after wrapping.
   #          $yRadii = 1/2 heights of points, after wrapping.
   #          $shapes = The shapes parameter, after wrapping.
   #          $sides = The sides parameter, after wrapping.
   #          $fillColors = The fillColors parameter, after wrapping.
   #          $borderColors = The borderColors parameter, after wrapping.
   #          $x = X coordinates of polygons to draw, NA separated
   #          $y = Y coordinates of polygons to draw, NA separated
   #          $xRange = lowest and highest x value drawn
   #          $yRange = lowest and highest y values drawn
   ###
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

# Merge data from multiple genes and multiple point sets into a data structure
# for plotting as a fusion expression graph. Uses the gene to exon rect data
# function above
getFusionExpressionPlotData <- function( genes, fusions, geneNames,
                                         geneSep= 500, fusionSep= 40, nameSep= 20) {
   ###
   #   genes = list of GRanges objects, by gene name
   #   fusions = list of vectors of fusion genomic positions, by gene name
   #   geneNames = list of gene names, must match those in genes and fusions
   ###

   if (DEBUG) {
      print( "In getFusionExpressionPlotData() with:");
      print(genes);
      if (is.null(fusions)) {
         print(NULL)
      } else {
         print(fusions)
      }
      print(geneNames);
   }

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


# Given the data from a merged fusion plot, plot it.
plot.fusionExpressionPlotData <- function ( data, title, file= 'temp.pdf',
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

   plot( data$xRange, data$yRange, type="n", axes=F, xlab="", ylab="", main=title, );

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

plot.fusionExpressionPair <- function( geneName1, geneName2, fusions1, fusions2, sample, cohortExpressionDF ) {

   if (DEBUG) {
      print("Called plot.fusionExpressionPair");
      print("geneName1:"); print(geneName1);
      print("geneName2:"); print(geneName2);
      print("fusions1:"); print(fusions1);
      print("fusions2:"); print(fusions2);
      print("sample:"); print(sample);
      print("cohortExpressionDF:"); print(cohortExpressionDF);

   }
   gr1 <- grFromCohortDf(cohortExpressionDF, geneName1, sample, column='relativeExpression')[[geneName1]];
   if (DEBUG) { print("Gene 1, before color mapping"); print( gr1 ); }
   gr2 <- grFromCohortDf(cohortExpressionDF, geneName2, sample, column='relativeExpression')[[geneName2]];
   if (DEBUG) { print("Gene 2, before color mapping"); print( gr2 ); }

   gr1 <- expressionGrColorMap( gr1 );
   if (DEBUG) { print("Gene 1, after color mapping"); print( gr1 ); }
   gr2 <- expressionGrColorMap( gr2 );
   if (DEBUG) { print("Gene 2, after color mapping"); print( gr2 ); }

   geneNames <- c(geneName1, geneName2);
   genes <- list();
   genes[[geneNames[1]]] <- gr1
   genes[[geneNames[2]]] <- gr2
   fusions <- list()
   fusions[[geneNames[1]]] <- fusions1;
   fusions[[geneNames[2]]] <- fusions2;
   data <- getFusionExpressionPlotData( genes, fusions, geneNames );


   plot.fusionExpressionPlotData( data, title=paste0( sample, " ", geneName1, "~", geneName2 ),
                                  file= paste0( sample, ".", geneName1, ".", geneName2, ".pdf" ) );
}

plot.fusionExpressionOne <- function( geneName1, fusions1, sample, cohortExpressionDF ) {

   if (DEBUG) {
      print(geneName1)
      if (is.null(fusions1)) {
         print(NULL)
      } else {
         print(fusions1)
      }
      print(sample)
   }

   gr <- grFromCohortDf(cohortExpressionDF, geneName1, sample, column='relativeExpression')[[geneName1]];

   gr <- expressionGrColorMap( gr );

   genes <- list();
   genes[[geneName1]] <- gr
   fusions <- NULL
   if (! is.null(fusions1)) {
      fusions <- list()
      fusions[[geneName1]] <- fusions1;
   }
   data <- getFusionExpressionPlotData( genes, fusions, geneName1 );


   plot.fusionExpressionPlotData( data, title=paste0( sample, "   ", geneName1 ),
                                  file= paste0( sample, ".", geneName1, ".pdf" ) );
}


demo.plot.fusionExpressionPlotData <- function() {
   data <- demo.getFusionExpressionPlotData();
   plot.fusionExpressionPlotData( data, title="Testing: clown + TMPRSS2", file = 'clown_TMPRSS2.pdf' );
}

demo.plot.fusionExpressionPlotData1 <- function() {
   data <- demo.getFusionExpressionPlotData1();
   plot.fusionExpressionPlotData( data, title="TMPRSS2", file = 'TMPRSS2.pdf' );
}

demo.plot.fusionExpressionPlotData4 <- function() {
   data <- demo.getFusionExpressionPlotData4();
   plot.fusionExpressionPlotData( data, title="clownA TMPRSS2 clownB clownC", file = 'clownA_TMPRSS2_clownB_clownC.pdf' );
}

demo.plot.fusionExpressionPlotData2 <- function() {
   data <- demo.getFusionExpressionPlotData2();
   plot.fusionExpressionPlotData( data, title="Testing: TMPRSS2 + TMPRSS2", file = 'TMPRSS2_TMPRSS2.pdf' );
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
#          expressionGr1 <- grColorMap(expressionGr1, cutPoints, column='relativeExpression' );
#          expressionGr2 <- grColorMap(expressionGr2, cutPoints, column='relativeExpression' );
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
   exons <- "42836479-42838080,42839661-42839813,42840323-42840465,42842575-42842670,42843733-42843908,42845252-42845423,42848504-42848547,42851099-42851209,42852403-42852529,42860321-42860440,42861434-42861520,42866283-42866505,42870046-42870116,42879877-42879992";
   iRanges <- IRanges(
      start= as.numeric(sub("-.+", "", unlist(strsplit(exons, ",")))),
      end=   as.numeric(sub(".+-", "", unlist(strsplit(exons, ","))))
   );
   gr <- GRanges(
      seqnames <- as.factor( "chr21" ),
      ranges   <- iRanges,
      strand   <- '-'
   );
   return( gr )
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
#       TMPRSS2$GRanges, 'color', dataToBinLabels(
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
#       ERG$GRanges, 'color', dataToBinLabels(
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
