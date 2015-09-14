context("gafGeneModels.R")

describe( "extractGeneModels()", {
   # This is difficult to test. Going with fragile solution of tiny input file
   # (testGaf) and expected output file (testGafGeneModels). Testing by exact
   # match of output file.
   #
   # Maintainance Note:
   #      If errors are discovered in the selection of records to
   # extract, representative records that reporduce the error should be added
   # to the input file, and when skipped by mistake, also to the output file

   testGaf <- 'data/mock.gaf'
   testGafContentProxy <- tools::md5sum( testGaf )
   testGafLineCount <- 15

   testGafGeneModels <- 'data/mock.gaf.extract.geneModels'
   testGafGeneModelsContentProxy <- tools::md5sum( testGafGeneModels )
   outputGeneModelsFile <- paste0( testGaf, ".geneModels" )

   describe( "Behavior with default parameters", {

      it( "outputs correctly named file in same directory as the GAF?", {
         extractGeneModels( testGaf )
         expect_true( file.exists( outputGeneModelsFile ))

         # Clean up after ourselves
         if ( file.exists( outputGeneModelsFile )) {
            unlink( outputGeneModelsFile );
         }
      })

      it( "selects all gene records and only genes records from the GAF?", {
         extractGeneModels( testGaf )
         gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
         # Need equivalent as md5 value named with filename.
         expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )

         # Clean up after ourselves
         if ( file.exists( outputGeneModelsFile )) {
            unlink( outputGeneModelsFile );
         }
      })

      it("returns the correct stats", {
         got <- extractGeneModels( testGaf )
         expect_is( got, "list")
         expect_equal( length(got), 11)

         expect_equal( got$gaf,             testGaf                 )
         expect_equal( got$gaf_real,        normalizePath( testGaf ))
         expect_equal( got$gaf_md5,         testGafContentProxy     )
         expect_equal( got$gaf_lines,       testGafLineCount        )
         expect_equal( got$uniqueGene,        TRUE )
         expect_equal( got$skipUnknownGene,   TRUE )
         expect_equal( got$gaf_models,        7    )
         expect_equal( got$gaf_models_unique, 4    )
         expect_equal( got$gaf_extract,      outputGeneModelsFile                 )
         expect_equal( got$gaf_extract_real, normalizePath( outputGeneModelsFile ))
         expect_equivalent( got$gaf_extract_md5,  testGafGeneModelsContentProxy   )

         # Clean up after ourselves
         if ( file.exists( outputGeneModelsFile )) {
            unlink( outputGeneModelsFile );
         }
      })

      describe("Errors with default parameters", {

         it( "stops with error if can't find specified input GAF?", {

            noSuchFile <- 'noSuchFile.gaf'
            noSuchOutputFile <- paste0( testGaf, ".geneModels" )

            wantErrorRE <- paste0( "Can't find the specified GAF: \"", noSuchFile, "\"" )
            expect_error( extractGeneModels( noSuchFile ), wantErrorRE )
            expect_false( file.exists( noSuchOutputFile ))
            expect_false( file.exists( outputGeneModelsFile ))
         })

         it( "stops with error if input GAF filename is unsafe/hacking attempt?", {
            # Only checks that some attempt to do this is implemented. Impossible
            # to test this for completeness. [TODO - implement in pure R to avoid
            # system call and need to test this here...]

            wantErrorRE <- 'Unsafe character in GAF filename!'
            expect_error( extractGeneModels('my.gaf;ls'), wantErrorRE )
            expect_false( file.exists( outputGeneModelsFile ))
         })

         it( "stops with error if output file already exists.", {
            extractGeneModels( testGaf )
            expect_true( file.exists( outputGeneModelsFile ))
            wantErrorRE <- paste0( "Output file already exists; ",
               "use force= TRUE to overwrite: \"",
               outputGeneModelsFile, "\""
            )
            expect_error( extractGeneModels( testGaf ), wantErrorRE )

            # Test the file wasn't corrupted during failed attempt to overwrite.
            gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
            expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )

            # Clean up after ourselves
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

      }) # END: Errors with default parameters

   }) # END: Behavior with default parameters

   describe( "The optional parameters", {

      describe( "The 'outfile=' parameter", {

         outFileName <- tempfile(
            pattern = "FusionExpressionPlotTestDeleteMe",
            tmpdir = tempdir(), fileext = ".geneModels"
         )

         it( 'determines output filename when specified', {
            extractGeneModels( testGaf, outFile= outFileName )
            expect_true( file.exists( outFileName ))
            gotGeneModelsContentProxy <- tools::md5sum( outFileName )
            expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )

            # Clean up after ourselves
            if ( file.exists( outFileName )) {
               unlink( outFileName );
            }
         })

         it("returns the correct stats", {
            got <- extractGeneModels( testGaf, outFile= outFileName )

            expect_is( got, "list")
            expect_equal( length(got), 11)

            expect_equal( got$gaf,             testGaf                 )
            expect_equal( got$gaf_real,        normalizePath( testGaf ))
            expect_equal( got$gaf_md5,         testGafContentProxy     )
            expect_equal( got$gaf_lines,       testGafLineCount        )
            expect_equal( got$uniqueGene,        TRUE )
            expect_equal( got$skipUnknownGene,   TRUE )
            expect_equal( got$gaf_models,        7    )
            expect_equal( got$gaf_models_unique, 4    )
            expect_equal( got$gaf_extract,          outFileName                    )
            expect_equal( got$gaf_extract_real,     normalizePath( outFileName    ))
            expect_equivalent( got$gaf_extract_md5, testGafGeneModelsContentProxy  )

            # Clean up after ourselves
            if ( file.exists( outFileName )) {
               unlink( outFileName );
            }
         })

         it( 'stops with error if specified output filename exists', {
            extractGeneModels( testGaf, outFile= outFileName )
            expect_true( file.exists( outFileName ))
            wantErrorRE <- paste0( "Output file already exists; ",
                                   "use force= TRUE to overwrite: \"",
                                   outFileName, "\""
            )
            expect_error( extractGeneModels( testGaf, outFile= outFileName ), wantErrorRE )

            # Test the file wasn't corrupted during failed attempt to overwrite.
            gotGeneModelsContentProxy <- tools::md5sum( outFileName )
            expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )

            # Clean up after ourselves
            if ( file.exists( outFileName )) {
               unlink( outFileName );
            }
         })

         it( "stops with error if output filename is unsafe/hacking attempt?", {
            # Only checks that some attempt to do this is implemented. Impossible
            # to test this for completeness. [TODO - implement in pure R to avoid
            # system call and need to test this here...]

            wantErrorRE <- 'Unsafe character in output geneModel filename!'
            expect_error( extractGeneModels(testGaf, outFile= 'my.gaf;ls'), wantErrorRE )
            expect_false( file.exists( 'my.gaf' ))
         })

      }) # END 'outfile='

      describe( "The 'force=' parameter", {

         it( "overwrites with warning if output file exists and force= TRUE.", {
            file.create( outputGeneModelsFile )
            expect_true( file.exists( outputGeneModelsFile ))
            expect_warning( extractGeneModels( testGaf, force= TRUE ),
               "Forcing overwrite of output file: \"data/mock.gaf.geneModels\""
            )
            gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
            expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )

            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

         it("returns the correct stats when force= TRUE", {
            got <- extractGeneModels( testGaf, force= TRUE )

            expect_is( got, "list")
            expect_equal( length(got), 11)

            expect_equal( got$gaf,             testGaf                 )
            expect_equal( got$gaf_real,        normalizePath( testGaf ))
            expect_equal( got$gaf_md5,         testGafContentProxy     )
            expect_equal( got$gaf_lines,       testGafLineCount        )
            expect_equal( got$uniqueGene,        TRUE )
            expect_equal( got$skipUnknownGene,   TRUE )
            expect_equal( got$gaf_models,        7    )
            expect_equal( got$gaf_models_unique, 4    )
            expect_equal( got$gaf_extract,      outputGeneModelsFile                 )
            expect_equal( got$gaf_extract_real, normalizePath( outputGeneModelsFile ))
            expect_equivalent( got$gaf_extract_md5,  testGafGeneModelsContentProxy   )

            # Clean up after ourselves
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

         it( "stops with error if output file exists and force= FALSE.", {
            extractGeneModels( testGaf )
            expect_true( file.exists( outputGeneModelsFile ))
            wantErrorRE <- paste0( "Output file already exists; ",
                                   "use force= TRUE to overwrite: \"",
                                   outputGeneModelsFile, "\""
            )
            expect_error( extractGeneModels( testGaf, force= FALSE ), wantErrorRE )

            # Clean up after ourselves
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }

         })

         it( "just output if file doesn't exists for force= TRUE or FALSE.", {

            expect_false( file.exists( outputGeneModelsFile ))
            extractGeneModels( testGaf, force= FALSE )
            expect_true( file.exists( testGaf ))
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
            expect_false( file.exists( outputGeneModelsFile ))
            extractGeneModels( testGaf, force= TRUE )
            expect_true( file.exists( outputGeneModelsFile ))
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

      }) # END: 'force='

      describe( "The 'uniqueGene=' and 'skipUnknownGene=' parameters", {

         it( "Ignores skipUnknownGene if uniqueGene=TRUE", {

            extractGeneModels( testGaf, uniqueGene= TRUE )
            expect_true( file.exists( outputGeneModelsFile ))
            gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
            expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }

            extractGeneModels( testGaf, uniqueGene= TRUE, skipUnknownGene= FALSE )
            expect_true( file.exists( outputGeneModelsFile ))
            gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
            expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }

            extractGeneModels( testGaf, uniqueGene= TRUE, skipUnknownGene= TRUE )
            expect_true( file.exists( outputGeneModelsFile ))
            gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
            expect_equivalent( gotGeneModelsContentProxy, testGafGeneModelsContentProxy )
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }

         })

         it( "generates warning if ignoring skipUnknownGene= TRUE", {
            wantWarnRE <- "uniqueGene=TRUE sets skipUnknownGene=TRUE"
            expect_warning(
               extractGeneModels( testGaf, uniqueGene= TRUE, skipUnknownGene= FALSE ),
               wantWarnRE
            )
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

         it( "Includes duplicate genes but not '?' genes when uniqueGene= FALSE and skipUnknownGene= TRUE", {
            expectGafGeneModels <- 'data/mock.gaf.extract.dupNoUnk.geneModels'
            expectGafGeneModelsContentProxy <- tools::md5sum( expectGafGeneModels )

            extractGeneModels( testGaf, uniqueGene= FALSE, skipUnknownGene= TRUE )
            gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
            expect_equivalent( gotGeneModelsContentProxy, expectGafGeneModelsContentProxy )
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

         it( "Includes duplicate genes and all '?' genes when uniqueGene= FALSE and skipUnknownGene= FALSE", {
            expectGafGeneModels <- 'data/mock.gaf.extract.dup.geneModels'
            expectGafGeneModelsContentProxy <- tools::md5sum( expectGafGeneModels )

            extractGeneModels( testGaf, uniqueGene= FALSE, skipUnknownGene= FALSE )
            gotGeneModelsContentProxy <- tools::md5sum( outputGeneModelsFile )
            expect_equivalent( gotGeneModelsContentProxy, expectGafGeneModelsContentProxy )
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

         it("returns the correct stats when uniqueGene= TRUE, skipUnknownGene= FALSE", {
            got <- extractGeneModels( testGaf, uniqueGene= TRUE, skipUnknownGene= FALSE )

            expect_is( got, "list")

            expect_equal( length(got), 11)

            expect_equal( got$gaf,             testGaf                 )
            expect_equal( got$gaf_real,        normalizePath( testGaf ))
            expect_equal( got$gaf_md5,         testGafContentProxy     )
            expect_equal( got$gaf_lines,       testGafLineCount        )
            expect_equal( got$uniqueGene,        TRUE )
            expect_equal( got$skipUnknownGene,   TRUE )
            expect_equal( got$gaf_models,        7    )
            expect_equal( got$gaf_models_unique, 4    )
            expect_equal( got$gaf_extract,      outputGeneModelsFile                 )
            expect_equal( got$gaf_extract_real, normalizePath( outputGeneModelsFile ))
            expect_equivalent( got$gaf_extract_md5,  testGafGeneModelsContentProxy   )

            # Clean up after ourselves
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

         it("returns the correct stats when uniqueGene= FALSE, skipUnknownGene= TRUE", {
            expectGafGeneModels <- 'data/mock.gaf.extract.dupNoUnk.geneModels'
            expectGafGeneModelsContentProxy <- tools::md5sum( expectGafGeneModels )
            got <- extractGeneModels( testGaf, uniqueGene= FALSE, skipUnknownGene= TRUE )

            expect_is( got, "list")
            expect_equal( length(got), 11)

            expect_equal( got$gaf,             testGaf                 )
            expect_equal( got$gaf_real,        normalizePath( testGaf ))
            expect_equal( got$gaf_md5,         testGafContentProxy     )
            expect_equal( got$gaf_lines,       testGafLineCount        )
            expect_equal( got$uniqueGene,        FALSE )
            expect_equal( got$skipUnknownGene,   TRUE )
            expect_equal( got$gaf_models,        7    )
            expect_equal( got$gaf_models_unique, NA    )
            expect_equal( got$gaf_extract,      outputGeneModelsFile                 )
            expect_equal( got$gaf_extract_real, normalizePath( outputGeneModelsFile ))
            expect_equivalent( got$gaf_extract_md5,  expectGafGeneModelsContentProxy   )

            # Clean up after ourselves
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

         it("returns the correct stats when uniqueGene= FALSE, skipUnknownGene= FALSE", {
            expectGafGeneModels <- 'data/mock.gaf.extract.dup.geneModels'
            expectGafGeneModelsContentProxy <- tools::md5sum( expectGafGeneModels )
            got <- extractGeneModels( testGaf, uniqueGene= FALSE, skipUnknownGene= FALSE )

            expect_is( got, "list")
            expect_equal( length(got), 11)

            expect_equal( got$gaf,             testGaf                 )
            expect_equal( got$gaf_real,        normalizePath( testGaf ))
            expect_equal( got$gaf_md5,         testGafContentProxy     )
            expect_equal( got$gaf_lines,       testGafLineCount        )
            expect_equal( got$uniqueGene,        FALSE )
            expect_equal( got$skipUnknownGene,   FALSE )
            expect_equal( got$gaf_models,        10    )
            expect_equal( got$gaf_models_unique, NA    )
            expect_equal( got$gaf_extract,      outputGeneModelsFile                 )
            expect_equal( got$gaf_extract_real, normalizePath( outputGeneModelsFile ))
            expect_equivalent( got$gaf_extract_md5,  expectGafGeneModelsContentProxy   )

            # Clean up after ourselves
            if ( file.exists( outputGeneModelsFile )) {
               unlink( outputGeneModelsFile );
            }
         })

      }) # END: 'uniqueGene=' and 'skipUnknownGene='

   }) # END: The optional parameters

}) # END: extractGeneModels()

describe( "extractTranscriptModels()", {

   testGaf <- 'data/mock.gaf'
   testGafTranscripts <- 'data/mock.gaf.extract.transcriptModels'
   testGafTranscriptsContentProxy <- tools::md5sum( testGafTranscripts )
   outputTranscriptsFile <- paste0( testGaf, ".transcriptModels" )

   describe( "Behavior with default parameters", {

      it( "outputs correctly named file in same directory as the GAF?", {
         extractTranscriptModels( testGaf )
         expect_true( file.exists( outputTranscriptsFile ))

         # Clean up after ourselves
         if ( file.exists( outputTranscriptsFile )) {
            unlink( outputTranscriptsFile );
         }
      })

      it( "selects all gene records and only genes records from the GAF?", {
         extractTranscriptModels( testGaf )
         gotTranscriptsContentProxy <- tools::md5sum( outputTranscriptsFile )
         # Need equivalent as md5 value named with filename.
         expect_equivalent( gotTranscriptsContentProxy, testGafTranscriptsContentProxy )

         # Clean up after ourselves
         if ( file.exists( outputTranscriptsFile )) {
            unlink( outputTranscriptsFile );
         }
      })

      it("returns the correct value", {
         got <- extractTranscriptModels( testGaf )
         expect_equal( got, 4)
         # Clean up after ourselves
         if ( file.exists( outputTranscriptsFile )) {
            unlink( outputTranscriptsFile );
         }
      })

      describe("Errors with default parameters", {

         it( "stops with error if can't find specified input GAF?", {

            noSuchFile <- 'noSuchFile.gaf'
            noSuchOutputFile <- paste0( testGaf, ".transcriptModels" )

            wantErrorRE <- paste0( "Can't find the specified GAF: \"", noSuchFile, "\"" )
            expect_error( extractTranscriptModels( noSuchFile ), wantErrorRE )
            expect_false( file.exists( noSuchOutputFile ))
            expect_false( file.exists( outputTranscriptsFile ))
         })

         it( "stops with error if input GAF filename is unsafe/hacking attempt?", {
            # Only checks that some attempt to do this is implemented. Impossible
            # to test this for completeness. [TODO - implement in pure R to avoid
            # system call and need to test this here...]

            wantErrorRE <- 'Unsafe character in GAF filename!'
            expect_error( extractTranscriptModels('my.gaf;ls'), wantErrorRE )
            expect_false( file.exists( outputTranscriptsFile ))
         })

         it( "stops with error if output file already exists.", {
            extractTranscriptModels( testGaf )
            expect_true( file.exists( outputTranscriptsFile ))
            wantErrorRE <- paste0( "Output file already exists; ",
                                   "use force= TRUE to overwrite: \"",
                                   outputTranscriptsFile, "\""
            )
            expect_error( extractTranscriptModels( testGaf ), wantErrorRE )

            # Test the file wasn't corrupted during failed attempt to overwrite.
            gotTranscriptsContentProxy <- tools::md5sum( outputTranscriptsFile )
            expect_equivalent( gotTranscriptsContentProxy, testGafTranscriptsContentProxy )

            # Clean up after ourselves
            if ( file.exists( outputTranscriptsFile )) {
               unlink( outputTranscriptsFile );
            }
         })

      }) # END: Errors with default parameters

   }) # END: Behavior with default parameters

   describe( "The optional parameters", {

      describe( "The 'outfile=' parameter", {

         outFileName <- tempfile(
            pattern = "FusionExpressionPlotTestDeleteMe",
            tmpdir = tempdir(), fileext = ".transcriptModels"
         )

         it( 'determines output filename when specified', {
            extractTranscriptModels( testGaf, outFile= outFileName )
            expect_true( file.exists( outFileName ))
            gotTranscriptsContentProxy <- tools::md5sum( outFileName )
            expect_equivalent( gotTranscriptsContentProxy, testGafTranscriptsContentProxy )

            # Clean up after ourselves
            if ( file.exists( outFileName )) {
               unlink( outFileName );
            }
         })

         it("returns the correct stats", {
            got <- extractTranscriptModels( testGaf, outFile= outFileName )

            expect_equal( got, 4)

            # Clean up after ourselves
            if ( file.exists( outFileName )) {
               unlink( outFileName );
            }
         })

         it( 'stops with error if specified output filename exists', {
            extractTranscriptModels( testGaf, outFile= outFileName )
            expect_true( file.exists( outFileName ))
            wantErrorRE <- paste0( "Output file already exists; ",
                                   "use force= TRUE to overwrite: \"",
                                   outFileName, "\""
            )
            expect_error( extractTranscriptModels( testGaf, outFile= outFileName ), wantErrorRE )

            # Test the file wasn't corrupted during failed attempt to overwrite.
            gotTranscriptsContentProxy <- tools::md5sum( outFileName )
            expect_equivalent( gotTranscriptsContentProxy, testGafTranscriptsContentProxy )

            # Clean up after ourselves
            if ( file.exists( outFileName )) {
               unlink( outFileName );
            }
         })

         it( "stops with error if output filename is unsafe/hacking attempt?", {
            # Only checks that some attempt to do this is implemented. Impossible
            # to test this for completeness. [TODO - implement in pure R to avoid
            # system call and need to test this here...]

            wantErrorRE <- 'Unsafe character in output transcriptModel filename!'
            expect_error( extractTranscriptModels(testGaf, outFile= 'my.gaf;ls'), wantErrorRE )
            expect_false( file.exists( 'my.gaf' ))
         })

      }) # END: The 'outfile=' parameter

      describe( "The 'force=' parameter", {

         it( "overwrites with warning if output file exists and force= TRUE.", {
            file.create( outputTranscriptsFile )
            expect_true( file.exists( outputTranscriptsFile ))
            expect_warning( extractTranscriptModels( testGaf, force= TRUE ),
                            "Forcing overwrite of output file: \"data/mock.gaf.transcriptModels\""
            )
            gotTranscriptsContentProxy <- tools::md5sum( outputTranscriptsFile )
            expect_equivalent( gotTranscriptsContentProxy, testGafTranscriptsContentProxy )

            if ( file.exists( outputTranscriptsFile )) {
               unlink( outputTranscriptsFile );
            }
         })

         it("returns the correct stats when force= TRUE", {
            got <- extractTranscriptModels( testGaf, force= TRUE )
            expect_equal(got, 4)

            # Clean up after ourselves
            if ( file.exists( outputTranscriptsFile )) {
               unlink( outputTranscriptsFile );
            }
         })

         it( "stops with error if output file exists and force= FALSE.", {
            extractTranscriptModels( testGaf )
            expect_true( file.exists( outputTranscriptsFile ))
            wantErrorRE <- paste0( "Output file already exists; ",
                                   "use force= TRUE to overwrite: \"",
                                   outputTranscriptsFile, "\""
            )
            expect_error( extractTranscriptModels( testGaf, force= FALSE ), wantErrorRE )

            # Clean up after ourselves
            if ( file.exists( outputTranscriptsFile )) {
               unlink( outputTranscriptsFile );
            }

         })

         it( "just output if file doesn't exists for force= TRUE or FALSE.", {

            expect_false( file.exists( outputTranscriptsFile ))
            extractTranscriptModels( testGaf, force= FALSE )
            expect_true( file.exists( testGaf ))
            if ( file.exists( outputTranscriptsFile )) {
               unlink( outputTranscriptsFile );
            }
            expect_false( file.exists( outputTranscriptsFile ))
            extractTranscriptModels( testGaf, force= TRUE )
            expect_true( file.exists( outputTranscriptsFile ))
            if ( file.exists( outputTranscriptsFile )) {
               unlink( outputTranscriptsFile );
            }
         })

      }) # END: The 'force=' parameter

   }) # END: The optional parameters

}) # END: extractTranscriptModels()

describe( "loadGafModels()", {

   describe( "parameter validation", {
      it("errors if file can't be found.", {
         modelFile <- "noSuchFile.I.hope"
         wantErrorRE <- paste0( "Can't find the specified GAF extract file: \"",
                                modelFile, "\"")
         expect_error( loadGafModels( modelFile ), wantErrorRE)
      })
   })
   describe( "Loading a unique gene model file", {

      modelFile <- 'data/mock.gaf.extract.geneModels'
      geneExonCounts <- c(14, 27, 11, 2)
      lastExonRows <- c(cumsum(geneExonCounts))
      firstExonRows <- c(1,lastExonRows[1:length(geneExonCounts) - 1] + 1)
      secondExonRows <- firstExonRows + 1

      df <- loadGafModels( modelFile )

      it( "is a data frame with correct columns and rows", {
         wantClass <- 'data.frame'
         expect_is(df, wantClass)

         wantRowCount <- sum(geneExonCounts)
         expect_equal( nrow(df), wantRowCount )

         wantColNames <- c("gene", "chr", "strand", "gstart", "gend",
                           "exon", "start", "end", "length")
         expect_equal( colnames(df), wantColNames)

      })

      it( 'has the correct "per model" data', {
         # gene, chr, strand, gstart, and gend are the same for all exons
         # in each model

         wantGene <- rep( c("TP53", "BRCA1", "SLC35E2", "SPANXB2"), geneExonCounts )
         expect_equal( df$gene, wantGene)

         wantChr <- rep( c("chr17", "chr17", "chr1", "chrX"), geneExonCounts )
         expect_equal( df$chr, wantChr)

         wantStrand <- rep( c("-", "-", "-", "+"), geneExonCounts )
         expect_equal( df$strand, wantStrand)

         wantGStart <- rep( c(7565097, 41196313, 1590990, 140084756), geneExonCounts )
         expect_equal( df$gstart, wantGStart)

         wantGEnd <- rep( c(7590863, 41322420, 1677431, 140085870), geneExonCounts )
         expect_equal( df$gend, wantGEnd)
      })

      it( 'gets all exon numbers correct', {
         # Exon numbers vary per
         wantExon <- c( 1:geneExonCounts[1], 1:geneExonCounts[2],
                        1:geneExonCounts[3], 1:geneExonCounts[4] )
         expect_equal( df$exon, wantExon)
      })

      it( 'gets first, second, and last values correct for per-exon values', {
         forRows <- c(firstExonRows[1], secondExonRows[1], lastExonRows[1])
         wantStarts  <- c(7565097, 7569421, 7590695)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(7565332, 7569562, 7590863)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[2], secondExonRows[2], lastExonRows[2])
         wantStarts  <- c(41196313, 41199660, 41322143)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(41197819, 41199720, 41322420)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[3], secondExonRows[3], lastExonRows[3])
         wantStarts  <- c(1590990, 1592940, 1677163)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(1591572, 1597458, 1677431)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[4], secondExonRows[4], lastExonRows[4])
         wantStarts  <- c(140084756, 140085596, 140085596)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(140084948, 140085870, 140085870)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )
      })
   }) # END: Unique gene model

   describe( "Loading a gene model file with dups (no unknowns)", {

      modelFile <- 'data/mock.gaf.extract.dupNoUnk.geneModels'

      geneExonCounts <- c(14, 27, 3, 10, 11, 2, 2)
      lastExonRows <- c(cumsum(geneExonCounts))
      firstExonRows <- c(1,lastExonRows[1:length(geneExonCounts)-1] + 1)
      secondExonRows <- firstExonRows + 1

      df <- loadGafModels( modelFile )

      it( "is a data frame with correct columns and rows", {
         wantClass <- 'data.frame'
         expect_is(df, wantClass)

         wantRowCount <- sum(geneExonCounts)
         expect_equal( nrow(df), wantRowCount )

         wantColNames <- c("gene", "chr", "strand", "gstart", "gend",
                           "exon", "start", "end", "length")
         expect_equal( colnames(df), wantColNames)

      })

      it( 'has the correct "per model" data', {
         # gene, chr, strand, gstart, and gend are the same for all exons
         # in each model

         wantGene <- rep( c("TP53", "BRCA1", "DUX4", "SLC35E2", "SLC35E2", "SPANXB2", "SPANXB2"), geneExonCounts )
         expect_equal( df$gene, wantGene)

         wantChr <- rep( c("chr17", "chr17", "GL000228.1", "chr1", "chr1", "chrX", "chrX"), geneExonCounts )
         expect_equal( df$chr, wantChr)

         wantStrand <- rep( c("-", "-", "+", "-", "-", "+", "+"), geneExonCounts )
         expect_equal( df$strand, wantStrand)

         wantGStart <- rep( c(7565097, 41196313, 95774, 1601214, 1590990, 140084756, 140096761), geneExonCounts )
         expect_equal( df$gstart, wantGStart)

         wantGEnd <- rep( c(7590863, 41322420, 97391, 1677431, 1677431, 140085870, 140097875), geneExonCounts )
         expect_equal( df$gend, wantGEnd)
      })

      it( 'gets all exon numbers correct', {
         # Exon numbers vary per
         wantExon <- c( 1:geneExonCounts[1], 1:geneExonCounts[2],
            1:geneExonCounts[3], 1:geneExonCounts[4], 1:geneExonCounts[5],
            1:geneExonCounts[6], 1:geneExonCounts[7])
         expect_equal( df$exon, wantExon)
      })

      it( 'gets first, second, and last values correct for per-exon values', {
         forRows <- c(firstExonRows[1], secondExonRows[1], lastExonRows[1])
         wantStarts  <- c(7565097, 7569421, 7590695)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(7565332, 7569562, 7590863)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[2], secondExonRows[2], lastExonRows[2])
         wantStarts  <- c(41196313, 41199660, 41322143)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(41197819, 41199720, 41322420)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[3], secondExonRows[3], lastExonRows[3])
         wantStarts  <- c(95774, 95787, 95937)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(95782, 95932, 97391)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[4], secondExonRows[4], lastExonRows[4])
         wantStarts  <- c(1601214, 1602948, 1677163)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(1601590, 1603068, 1677431)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[5], secondExonRows[5], lastExonRows[5])
         wantStarts  <- c(1590990, 1592940, 1677163)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(1591572, 1597458, 1677431)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[6], secondExonRows[6], lastExonRows[6])
         wantStarts  <- c(140084756, 140085596, 140085596)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(140084948, 140085870, 140085870)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[7], secondExonRows[7], lastExonRows[7])
         wantStarts  <- c(140096761, 140097601, 140097601)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(140096953, 140097875, 140097875)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )
      })

   }) # END: Dups gene model (no unk)

   describe( "Loading a gene model file with dups (no unknowns)", {
      modelFile <- 'data/mock.gaf.extract.dup.geneModels'
      geneExonCounts <- c(14, 27, 3, 10, 11, 6, 6, 6, 2, 2)
      lastExonRows <- c(cumsum(geneExonCounts))
      firstExonRows <- c(1,lastExonRows[1:length(geneExonCounts)-1] + 1)
      secondExonRows <- firstExonRows + 1

      df <- loadGafModels( modelFile )

      it( "is a data frame with correct columns and rows", {
         wantClass <- 'data.frame'
         expect_is(df, wantClass)

         wantRowCount <- sum(geneExonCounts)
         expect_equal( nrow(df), wantRowCount )

         wantColNames <- c("gene", "chr", "strand", "gstart", "gend",
                           "exon", "start", "end", "length")
         expect_equal( colnames(df), wantColNames)
      })

      it( 'has the correct "per model" data', {
         # gene, chr, strand, gstart, and gend are the same for all exons
         # in each model

         wantGene <- rep( c("TP53", "BRCA1", "DUX4", "SLC35E2", "SLC35E2",
            "?", "?", "?", "SPANXB2", "SPANXB2"), geneExonCounts )
         expect_equal( df$gene, wantGene)

         wantChr <- rep( c("chr17", "chr17", "GL000228.1", "chr1", "chr1",
            "chr9", "chr15", "chr15", "chrX", "chrX"), geneExonCounts )
         expect_equal( df$chr, wantChr)

         wantStrand <- rep( c("-", "-", "+", "-", "-", "+", "+", "+", "+", "+"),
            geneExonCounts )
         expect_equal( df$strand, wantStrand)

         wantGStart <- rep( c(7565097, 41196313, 95774, 1601214, 1590990,
               69425668,82647286,83023773,140084756, 140096761), geneExonCounts )
         expect_equal( df$gstart, wantGStart)

         wantGEnd <- rep( c(7590863, 41322420, 97391, 1677431, 1677431,
            69448861, 82708202, 83084727, 140085870, 140097875), geneExonCounts )
         expect_equal( df$gend, wantGEnd)
      })

      it( 'gets all exon numbers correct', {
         # Exon numbers vary per
         wantExon <- c( 1:geneExonCounts[1], 1:geneExonCounts[2],
            1:geneExonCounts[3], 1:geneExonCounts[4], 1:geneExonCounts[5],
            1:geneExonCounts[6], 1:geneExonCounts[7], 1:geneExonCounts[8],
            1:geneExonCounts[9], 1:geneExonCounts[10])
         expect_equal( df$exon, wantExon)
      })

      it( 'gets first, second, and last values correct for per-exon values', {
         forRows <- c(firstExonRows[1], secondExonRows[1], lastExonRows[1])
         wantStarts  <- c(7565097, 7569421, 7590695)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(7565332, 7569562, 7590863)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[2], secondExonRows[2], lastExonRows[2])
         wantStarts  <- c(41196313, 41199660, 41322143)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(41197819, 41199720, 41322420)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[3], secondExonRows[3], lastExonRows[3])
         wantStarts  <- c(95774, 95787, 95937)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(95782, 95932, 97391)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[4], secondExonRows[4], lastExonRows[4])
         wantStarts  <- c(1601214, 1602948, 1677163)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(1601590, 1603068, 1677431)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[5], secondExonRows[5], lastExonRows[5])
         wantStarts  <- c(1590990, 1592940, 1677163)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(1591572, 1597458, 1677431)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[6], secondExonRows[6], lastExonRows[6])
         wantStarts  <- c(69425668, 69432635, 69448743)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(69425817, 69432764, 69448861)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[7], secondExonRows[7], lastExonRows[7])
         wantStarts  <- c(82647286, 82692871, 82707733)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(82647542, 82692927, 82708202)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[8], secondExonRows[8], lastExonRows[8])
         wantStarts  <- c(83023773, 83069396, 83084258)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(83024029, 83069452, 83084727)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[9], secondExonRows[9], lastExonRows[9])
         wantStarts  <- c(140084756, 140085596, 140085596)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(140084948, 140085870, 140085870)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[10], secondExonRows[10], lastExonRows[10])
         wantStarts  <- c(140096761, 140097601, 140097601)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(140096953, 140097875, 140097875)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )
      })
   }) # END: Dups gene model (no unk)

   describe( "Loading a transcript model file", {
      modelFile <- 'data/mock.gaf.extract.transcriptModels'
      transcriptExonCounts <- c(2, 5, 7, 11)
      lastExonRows <- c(cumsum(transcriptExonCounts))
      firstExonRows <- c(1,lastExonRows[1:length(transcriptExonCounts)-1] + 1)
      secondExonRows <- firstExonRows + 1

      df <- loadGafModels( modelFile )

      it( "is a data frame with correct columns and rows", {
         wantClass <- 'data.frame'
         expect_is(df, wantClass)

         wantRowCount <- sum(transcriptExonCounts)
         expect_equal( nrow(df), wantRowCount )

         wantColNames <- c("gene", "chr", "strand", "gstart", "gend",
                           "exon", "start", "end", "length")
         expect_equal( colnames(df), wantColNames)
      })

      it( 'has the correct "per model" data', {
         # gene, chr, strand, gstart, and gend are the same for all exons
         # in each model

         wantTranscript <- rep( c("uc001acb.1", "uc009vjs.1", "uc002gig.1", "uc002gim.2"),
                          transcriptExonCounts )
         expect_equal( df$gene, wantTranscript)

         wantChr <- rep( c("chr1", "chr1", "chr17", "chr17"), transcriptExonCounts )
         expect_equal( df$chr, wantChr)

         wantStrand <- rep( c("+", "+", "-", "-"), transcriptExonCounts )
         expect_equal( df$strand, wantStrand)

         wantGStart <- rep( c(896829, 995083, 7565097, 7571720), transcriptExonCounts )
         expect_equal( df$gstart, wantGStart)

         wantGEnd <- rep( c(897858, 997436, 7579937, 7590863), transcriptExonCounts )
         expect_equal( df$gend, wantGEnd)
      })

      it( 'gets all exon numbers correct', {
         # Exon numbers vary per
         wantExon <- c( 1:transcriptExonCounts[1], 1:transcriptExonCounts[2],
                        1:transcriptExonCounts[3], 1:transcriptExonCounts[4] )
         expect_equal( df$exon, wantExon)
      })

      it( 'gets first, second, and last values correct for per-exon values', {
         forRows <- c(firstExonRows[1], secondExonRows[1], lastExonRows[1])
         wantStarts  <- c(896829, 897206, 897206)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(897130, 897858, 897858)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[2], secondExonRows[2], lastExonRows[2])
         wantStarts  <- c(995083, 995657, 997229)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(995162, 995773, 997436)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[3], secondExonRows[3], lastExonRows[3])
         wantStarts  <- c(7565097, 7577499, 7579839)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(7565332, 7577608, 7579937)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )

         forRows <- c(firstExonRows[4], secondExonRows[4], lastExonRows[4])
         wantStarts  <- c(7571720, 7573927, 7590695)
         expect_equal( df$start[forRows], wantStarts )
         wantEnds    <- c(7573008, 7574033, 7590863)
         expect_equal( df$end[forRows], wantEnds )
         wantLengths <- abs(wantStarts-wantEnds) + 1
         expect_equal( df$length[forRows], wantLengths )
      })

   }) # END: transcript model

}) # END: loadGafModels()