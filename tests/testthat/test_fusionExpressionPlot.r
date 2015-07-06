context("Functions in the fusionExpressionPlot.R file.")

#
# grNew
#

test_that( 'grNew creates correct genomicRanges objects, everything specified', {
   start <- c( 101, 301, 511, 513, 813 )
   end <- c( 200, 500, 511, 612, 822 )
   chr <- 'chr1'
   strand <- '+'
   got <- grNew( start= start, end= end, chr=chr, strand= strand )
   want <- GRanges(
      ranges= IRanges( start= c( 101, 301, 511, 513, 813 ), end= c( 200, 500, 511, 612, 822 )),
      seqnames= c( "chr1" ),
      strand= c( "+" )
   );
   expect_is(got, 'GRanges')
   expect_equal(got, want)
   expect_equal( start(got), start )
   expect_equal( end(got), end )
   expect_equal( as.character( seqnames(got) ), rep( chr, length(start) ))
   expect_equal( as.character( strand(got) ), rep( strand, length(start) ))
})

test_that( 'grNew, one exon, strand by default', {
   start <- 101
   end <- 200
   chr <- 'chrX'
   strand <- '*'

   got <- grNew( start= start, end= end, chr=chr )

   want <- GRanges(
      ranges= IRanges( start= c( 101 ), end= c( 200 )),
      seqnames= c( "chrX" ),
      strand= c( "*" )
   );
   expect_is(got, 'GRanges')
   expect_equal(got, want)
   expect_equal( start(got), start )
   expect_equal( end(got), end )
   expect_equal( as.character( seqnames(got) ), rep( chr, length(start) ))
   expect_equal( as.character( strand(got) ), rep( strand, length(start) ))
})

#
# genomicToModelCoordinates()
#

# Exons: 101-200,301-500,511-511,513-612,813-822
#          100     200      1      100     10
# Introns:   201-300,501-510,512-512,613-812
#              100     10       1      200
# Models:  1-100,  201-400,    501-501,    602-701,    802-811
#          1-100,10101-10300,20301-20301,30302-30401,40402-40411
#          1-100,102-301,303-303,305-404,406-415

makeGrPos <- function() {
   gr <- grNew(
      start=  c( 101, 301, 511, 513, 813 ),
      end=    c( 200, 500, 511, 612, 822 ),
      chr=    "chr1",
      strand= "+"
   );
   return(gr);
}

test_that( 'genomicToModelCoordinates; coordinates <= first +exon are ok.', {
   got = genomicToModelCoordinates( c(1, 100, 101), makeGrPos(), 100 );
   want = c(0, 0, 1)
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(1, 100, 101), makeGrPos(), 10000 );
   want = c(0, 0, 1)
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(1, 100, 101), makeGrPos(), 1 );
   want = c(0, 0, 1)
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinates >= last +exon are ok.', {
   got = genomicToModelCoordinates( c(822, 823, 10000), makeGrPos(), 100 );
   want = c(811, 812, 812);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(822, 823, 10000), makeGrPos(), 10000 );
   want = c(40411, 40412, 40412);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(822, 823, 10000), makeGrPos(), 1 );
   want = c(415, 416, 416);
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinates in +exons are ok.', {
   # first exon/exon = gap
   got = genomicToModelCoordinates( c(101, 150, 200), makeGrPos(), 100 );
   want = c(1, 50, 100);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(101, 150, 200), makeGrPos(), 10000 );
   want = c(1, 50, 100);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(101, 150, 200), makeGrPos(), 1 );
   want = c(1, 50, 100);
   expect_equal(want, got)

   # internal exon/exon larger than gap
   got = genomicToModelCoordinates( c(301, 400, 500), makeGrPos(), 100 );
   want = c(201, 300, 400);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(301, 400, 500), makeGrPos(), 10000 );
   want = c(10101, 10200, 10300);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(301, 400, 500), makeGrPos(), 1 );
   want = c(102, 201, 301);
   expect_equal(want, got)

   # last exon/exon smaller than gap
   got = genomicToModelCoordinates( c(813, 818, 822), makeGrPos(), 100 );
   want = c(802, 807, 811);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(813, 818, 822), makeGrPos(), 10000 );
   want = c(40402, 40407, 40411);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(813, 818, 822), makeGrPos(), 1 );
   want = c(406, 411, 415);
   expect_equal(want, got)

   # exon of size 1/single value works
   got = genomicToModelCoordinates( 511, makeGrPos(), 100 );
   want = 501;
   expect_equal(want, got)

   got = genomicToModelCoordinates( 511, makeGrPos(), 10000 );
   want = 20301;
   expect_equal(want, got)

   got = genomicToModelCoordinates( 511, makeGrPos(), 1 );
   want = 303;
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinates in +gaps are ok.', {
   # first gap ( = gapWidth)
   got = genomicToModelCoordinates( c(201, 250, 300), makeGrPos(), 100 );
   want = c(101, 150, 200);
   expect_equal(want, got)

   # internal gap
   got = genomicToModelCoordinates( c(501, 502, 505, 510), makeGrPos(), 100 );
   want = c(405, 415, 445, 495);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(501, 502, 505, 510), makeGrPos(), 1 );
   want = c(302, 302, 302, 302);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(501, 502, 505, 510), makeGrPos(), 10000 );
   want = c(10800, 11800, 14800, 19800);
   expect_equal(want, got)

   # last gap ( > gapWidth)
   got = genomicToModelCoordinates( c(613, 712, 812), makeGrPos(), 100 );
   want = c(702, 751, 801);
   expect_equal(want, got)

   got = genomicToModelCoordinates( c(613, 712, 812), makeGrPos(), 1 );
   want = c(405, 405, 405);
   expect_equal(want, got)

   # borders, size 1 intron
   got = genomicToModelCoordinates( c(511, 512, 513), makeGrPos(), 100 );
   want = c(501, 551, 602);
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinate extremes ok, +model?', {

   # Exons: 101-100100,200101-200200
   #          100000        100
   # Introns:   100101-200100
   #                100000

   # Model:   1-100000,100101-100200
   gr <- grNew(
      start=  c(    101, 200101 ),
      end=    c( 100100, 200200 ),
      chr=    "chr1",
      strand= "+"
   );

   # Before
   got = genomicToModelCoordinates(
      c(1, 100, 101), gr, 100 );
   want = c(0,   0, 1);
   expect_equal(want, got)

   # Exon 1
   got = genomicToModelCoordinates(
      c(100, 101, 102, 5100, 100099, 100100, 100101), gr, 100 );
   want = c(  0,   1,  2,  5000,  99999, 100000, 100001);
   expect_equal(want, got)

   # Intron 1
   got = genomicToModelCoordinates(
      c(100100, 100101, 150101, 200100, 200101), gr, 100 );
   want = c(100000, 100001, 100051, 100100, 100101);
   expect_equal(want, got)

   # Exon 2
   got = genomicToModelCoordinates(
      c(200100, 200101, 200102, 200199, 200200, 200201), gr, 100 );
   want = c(100100, 100101, 100102, 100199, 100200, 100201);
   expect_equal(want, got)

   # After
   got = genomicToModelCoordinates(
      c(200200, 200201, 200202), gr, 100 );
   want = c(100200, 100201, 100201);
   expect_equal(want, got)

})

# Exons: 101-200,301-500,511-511,513-612,813-822
#          100     200      1      100     10
# Introns:   201-300,501-510,512-512,613-812
#              100     10       1      200

# Model:   811-712,611-412,311-311,210-111,10-1

makeGrNeg <- function() {
   gr <- grNew(
      start=  c( 101, 301, 511, 513, 813 ),
      end=    c( 200, 500, 511, 612, 822 ),
      chr=    "chr1",
      strand= "-"
   );
   return(gr);
}

test_that( 'genomicToModelCoordinates; coordinates <= first -exon are ok.', {

   got = genomicToModelCoordinates( c(822, 823, 10000), makeGrNeg(), 100 );
   want = c(1, 0, 0);
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinates >= last -exon are ok.', {
   got = genomicToModelCoordinates( c(1, 100, 101), makeGrNeg(), 100 );
   want = c(812, 812, 811);
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinates in -exons are ok.', {
   # last exon/exon = gap
   got = genomicToModelCoordinates( c(101, 150, 200), makeGrNeg(), 100 );
   want = c(811, 762, 712);
   expect_equal(want, got)

   # internal exon/exon larger than gap
   got = genomicToModelCoordinates( c(301, 400, 500), makeGrNeg(), 100 );
   want = c(611, 512, 412);
   expect_equal(want, got)

   # first exon/exon smaller than gap
   got = genomicToModelCoordinates( c(813, 818, 822), makeGrNeg(), 100 );
   want = c(10, 5, 1);
   expect_equal(want, got)

   # exon of size 1/single value works
   got = genomicToModelCoordinates( 511, makeGrNeg(), 100 );
   want = 311;
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinates in -gaps are ok.', {
   # first gap ( = gapWidth)
   got = genomicToModelCoordinates( c(201, 250, 300), makeGrNeg(), 100 );
   want = c(711, 662, 612);
   expect_equal(want, got)

   # internal gap ( < gapWidth)
   got = genomicToModelCoordinates( c(501, 502, 505, 510), makeGrNeg(), 100 );
   want = c(407, 397, 367, 317);
   expect_equal(want, got)

   # last gap ( > gapWidth)
   got = genomicToModelCoordinates( c(613, 712, 812), makeGrNeg(), 100 );
   want = c(110, 61, 11);
   expect_equal(want, got)

   # borders, size 1 intron
   got = genomicToModelCoordinates( c(511, 512, 513), makeGrNeg(), 100 );
   want = c(311, 261, 210);
   expect_equal(want, got)
})

test_that( 'genomicToModelCoordinates; coordinate extremes ok, -model?', {

   # Exons: 101-100100,200101-200200
   #          100000        100
   # Introns:   100101-200100
   #                100000

   # Model:   1-100,201-100200
   gr <- grNew(
      start=  c(    101, 200101 ),
      end=    c( 100100, 200200 ),
      chr=    "chr1",
      strand= "-"
   );

   # Before
   got = genomicToModelCoordinates(
          c(     1,    100,    101), gr, 100 );
   want = c(100201, 100201, 100200);
   expect_equal(want, got)

   # Exon 1
   got = genomicToModelCoordinates(
          c(   100,    101,     102,  50100, 100099, 100100, 100101), gr, 100 );
   want = c(100201, 100200,  100199,  50201,  202, 201, 200);
   expect_equal(want, got)

   # Intron 1
   got = genomicToModelCoordinates(
          c(100100, 100101, 150101, 200100, 200101), gr, 100 );
   want = c(201, 200, 150, 101, 100);
   expect_equal(want, got)

   # Exon 2
   got = genomicToModelCoordinates(
       c(200100, 200101, 200102, 200199, 200200, 200201), gr, 100 );
   want = c(101,    100,     99,      2,      1,      0);
   expect_equal(want, got)

   # After
   got = genomicToModelCoordinates(
          c(200200, 200201, 200202), gr, 100 );
   want = c(1,           0,      0);
   expect_equal(want, got)

})

#
# sortDataIntoBins()
#

test_that( 'sortDataIntoBins; normal use ok?', {
   ends = c(-Inf, -1, -0.01, 0, 0.01, 1, Inf)
   labels = c('---', '--', '-', '+', '++', '+++')

   data = c(-Inf, -3e5, -1.001, -1, -0.5, -0.01, -0.001, -0, 0)
   got = sortDataIntoBins(data, ends, labels)
   want = c('---', '---', '---', '---', '--', '--', '-', '-', '-')
   expect_equal(got, want)

   data = c(Inf, 3e5, 1.001, 1, 0.5, 0.01, 0.001, +0, 0)
   got = sortDataIntoBins(data, ends, labels)
   want = c('+++', '+++', '+++', '++', '++', '+', '+', '-', '-')
   expect_equal(got, want)

   data = 7
   got = sortDataIntoBins(data, ends, labels)
   want = '+++'
   expect_equal(got, want)
})

test_that( 'sortDataIntoBins; missing/out of range values?', {
   ends = c(-10, -1, 0, 1, 10)
   labels = c('--', '-', '+', '++')
   data = c(NA_integer_, 2, 23)

   got = sortDataIntoBins(data, ends, labels)
   want = c(NA_character_, '++', NA_character_)
   expect_equal(got, want)
})

test_that( 'sortDataIntoBins; sorts input bin values?', {
   ends = c(Inf,-Inf,1,-1)
   labels = c('-', '0', '+' )
   data = c(-Inf,-5,-1,-0.5,0,0.5,1,5,Inf)

   got = sortDataIntoBins(data, ends, labels)
   want = c('-','-','-','0','0','0','0','+','+')
   expect_equal(got, want)
})

#
# grAddColumn()
#

test_that( 'grAddColumn; normal use ok?', {
   gr <- grNew(
      start=  c( 101, 301, 511, 513, 813 ),
      end=    c( 200, 500, 511, 612, 822 ),
      chr=    "chr1",
      strand= "+"
   );
   exonFillColor <- c('red', 'orange', 'yellow', 'green', 'blue')
   want <- GRanges(
      ranges= IRanges( start= c( 101, 301, 511, 513, 813 ), end= c( 200, 500, 511, 612, 822 )),
      seqnames= c( "chr1" ),
      strand= c( "+" ),
      exonFillColor= c('red', 'orange', 'yellow', 'green', 'blue')
   );
   got <- grAddColumn( gr, 'exonFillColor', exonFillColor);

   # manual creation is the same as adding column
   expect_identical(got, want);

   # can recover the vector from the object
   expect_equal(got$exonFillColor, exonFillColor);

   # Not affecting the original, but creating a new copy.
   expect_null(gr$exonFillColor);

})

test_that( 'grAddColumn; wrapping data values works?', {
   gr <- grNew(
      start=  c( 101, 301, 511, 513, 813, 1000 ),
      end=    c( 200, 500, 511, 612, 822, 1100 ),
      chr=    "chr1",
      strand= "+"
   );

   oneWrapped <- grAddColumn(gr, 'col_1', 1);
   expect_equal( oneWrapped$col_1, rep(1,6) );

   threeWrapped <- grAddColumn( gr, 'col_3', c('a', 'b', 'c') );
   expect_equal( threeWrapped$col_3, rep( c('a', 'b', 'c'), 2 ));
})

#
# grLabelMap()
#

test_that( 'grLabelMap; adds or returns correct column for default pallette', {
   fail( 'Test column returned not implemented' )
   fail( 'Test column added not implemented' )
   fail( 'Test original gr unaffected not implemented' )
})

test_that( 'grLabelMap; adds or returns correct column for given pallette', {
   fail( 'Test column returned not implemented' )
   fail( 'Test column added not implemented' )
   fail( 'Test original gr unaffected not implemented' )
})

test_that( 'grLabelMap; adds or returns correct column for given labels', {
   fail( 'Test column returned not implemented' )
   fail( 'Test column added not implemented' )
   fail( 'Test original gr unaffected not implemented' )
})

#
# grFromCohortDf()
#

test_that( 'grLabelMap from cohort with one sample, one gene', {
   fail( 'Returns all genes if unfiltered, with correct granges objects' )
   fail( 'Returns all and only selected genes if filtered, with correct granges objects' )
})

test_that( 'grLabelMap from cohort with two samples, two gene', {
   fail( 'Returns all genes for 1st sample if unfiltered, with correct granges objects' )
   fail( 'Returns all and only selected genes for 1st sample if filtered, with correct granges objects' )
   fail( 'Returns all genes for 2nd sample if unfiltered, with correct granges objects' )
   fail( 'Returns all and only selected genes for 2nd sample if filtered, with correct granges objects' )
})

#
# grToRect()
#

test_that( 'grToRect generates correct rectangle, all defaults', {
   fail( 'Correct xStarts' )
   fail( 'Correct xEnds' )
   fail( 'Correct yBottoms' )
   fail( 'Correct yTops' )
   fail( 'Correct fill colors' )
   fail( 'Correct border colors' )
   fail( 'Correct xRanges' )
   fail( 'Correct yRanges' )
})

test_that( 'grToRect generates correct rectangle, just one exon', {
   fail( 'Correct xStarts' )
   fail( 'Correct xEnds' )
   fail( 'Correct yBottoms' )
   fail( 'Correct yTops' )
   fail( 'Correct fill colors' )
   fail( 'Correct border colors' )
   fail( 'Correct xRanges' )
   fail( 'Correct yRanges' )
})

test_that( 'grToRect generates correct rectangle, new column names', {
   fail( 'Correct xStarts' )
   fail( 'Correct xEnds' )
   fail( 'Correct yBottoms' )
   fail( 'Correct yTops' )
   fail( 'Correct fill colors' )
   fail( 'Correct border colors' )
   fail( 'Correct xRanges' )
   fail( 'Correct yRanges' )
})

test_that( 'grToRect generates correct rectangle, manual values', {
   fail( 'Correct xStarts' )
   fail( 'Correct xEnds' )
   fail( 'Correct yBottoms' )
   fail( 'Correct yTops' )
   fail( 'Correct fill colors' )
   fail( 'Correct border colors' )
   fail( 'Correct xRanges' )
   fail( 'Correct yRanges' )
})

test_that( 'grToRect generates correct rectangle, manual values wrap', {
   fail( 'Correct xStarts' )
   fail( 'Correct xEnds' )
   fail( 'Correct yBottoms' )
   fail( 'Correct yTops' )
   fail( 'Correct fill colors' )
   fail( 'Correct border colors' )
   fail( 'Correct xRanges' )
   fail( 'Correct yRanges' )
})

test_that( 'grToRect generates correct rectangle, manual values and new column names', {
   fail( 'Correct xStarts' )
   fail( 'Correct xEnds' )
   fail( 'Correct yBottoms' )
   fail( 'Correct yTops' )
   fail( 'Correct fill colors' )
   fail( 'Correct border colors' )
   fail( 'Correct xRanges' )
   fail( 'Correct yRanges' )
})

makeGrFull <- function() {
   gr <- grNew(
      start=  c( 101, 301, 511, 513, 813 ),
      end=    c( 200, 500, 511, 612, 822 ),
      chr=    "chr1",
      strand= "+"
   );
   gr <- grAddColumn(gr, 'exonFillColor', c('red', 'orange', 'yellow', 'green', 'blue'));
   gr <- grAddColumn(gr, 'exonBorderColor', c('blue', 'green', 'black', 'orange', 'red'));
   gr <- grAddColumn(gr, 'exonHeight', c(50, 100, 200, 100, 50));
   return(gr);
}

