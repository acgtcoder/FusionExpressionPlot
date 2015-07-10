context("Unit tests for functions in the grUtils.R file.")

describe( 'grNew()', {
   it('creates a genomicRanges object if everything specified', {
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
   it( 'Can creates a one-exon genomicRanges object with correct default strand', {
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
})

describe( 'grAddColumn()', {
   it( 'Adds a columns as specified to a GRanges object', {
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
   it( 'can wrapp data columns if needed', {
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
})

describe( 'grFromElementString()', {
   it( 'Creates correct gr from string using default delimiters', {
      got <- grFromElementString( "1-100,200-300,400-500", 'chrX')
      want <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chrX', strand = '*');
      expect_equal(got,want);
   })
   it( 'Creates correct gr from string using default delimiters when string has spaces in it around the delimiters', {
      got <- grFromElementString( "1 - 100 , 200 -300 ,400- 500, 600\t - \t\t 700", 'chrX')
      want <- grNew( start=c(1,200,400,600), end=c(100,300,500,700), chr='chrX', strand = '*');
      expect_equal(got,want);
   })
   it( 'Creates correct gr from string with 1 exon only', {
      got <- grFromElementString( "1-100", 'chr1', strand='-')
      want <- grNew( start=c(1), end=c(100), chr='chr1', strand = '-');
      expect_equal(got,want);
   })

   it( 'Creates correct gr from string using alternate delimiters', {
      got <- grFromElementString( "1..100; 200..300; 400..500", 'chr12', strand='+', fromDelim='\\.\\.', betweenDelim='; ')
      want <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chr12', strand = '+');
      expect_equal(got,want);
   })
})

describe( 'grFromLocationString()', {
   it('creates correct gr from string using default captureRE', {
      got <- grFromLocationString( "chr1:1-100,200-300,400-500:*")
      want <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chr1', strand = '*');
      expect_equal(got,want);
   })

   it('creates correct gr from string using new captureRE', {
      captureRE <- '(?<chr>chr.+?)\\s*\\(\\s*(?<strand>[-+*])\\s*\\)\\s*(?<elements>.+)'
      got <- grFromLocationString( captureRE=captureRE, "chrM_rcs(+) 1-100, 200-300, 400-500" )
      want <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chrM_rcs', strand = '+');
      expect_equal(got,want);
   })

   it('creates correct gr from string if specifying new elementString parameters', {
      got <- grFromLocationString( "chr1:1..100;200..300;400..500:-", betweenDelim=";", fromDelim='\\.\\.')
      want <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chr1', strand = '-');
      expect_equal(got,want);
   })
})

describe("grToElementString()", {
   it("it outputs an element string describing the GRanges object using default between delimiter ',' and from delimiter '-'", {
      gr <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chrM_rcs', strand = '+');
      got <- grToElementString( gr )
      expect_is( got, "character")
      want <- "1-100,200-300,400-500"
      expect_equal(got, want)
   })
   it("it dies if the GRanges object has more than one chromosome", {
      gr1 <- grNew( start=c(1,200), end=c(100,300), chr='chr1', strand = '+');
      gr2 <- grNew( start=c(1,200), end=c(100,300), chr='chr2', strand = '+');
      gr <- suppressWarnings(c(gr1,gr2))
      # ensure actually has two chromosomes - testing the test
      expect_equal(length(unique(seqnames(gr))), 2)

      wantErrorRE <- "GRange object has more than one chromosome."
      expect_error( grToElementString( gr ), wantErrorRE )
   })
   it("it dies if the GRanges object has more than one strand", {
      gr1 <- grNew( start=c(1,200), end=c(100,300), chr='chr1', strand = '+');
      gr2 <- grNew( start=c(1,200), end=c(100,300), chr='chr1', strand = '-');
      gr <- suppressWarnings(c(gr1,gr2))
      # ensure actually has two strands - testing the test
      expect_equal(length(unique(strand(gr))), 2)

      wantErrorRE <- "GRange object has more than one strand."
      expect_error( grToElementString( gr ), wantErrorRE)
   })
   it("it returns empty if the GRanges object has no elements", {
      gr <- GRanges()
      expect_equal( length(gr), 0)
      want <- ""
      got <- grToElementString( gr )
      expect_equal(got, want)
   })
   it("allows between and from delimiters to be specified, including blanks", {
      gr <- grNew( start=c(1,200,400), end=c(100,300,500), chr='chrM_rcs', strand = '+');
      got <- grToElementString( gr, betweenDelim='; ', fromDelim=' ' )
      want <- "1 100; 200 300; 400 500"
      expect_equal(got, want)
   })
})