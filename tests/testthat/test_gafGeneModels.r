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