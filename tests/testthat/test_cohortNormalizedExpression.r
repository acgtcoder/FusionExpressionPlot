context("cohortNormalizedExpression.R")

describe( "loadCohortDefinition()", {

   cohortArray.default <- c(
      "# This is a fake cohort.",
      "# It is used only for testing",
      paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
      paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
      paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
      paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
      paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
   )
   cohortFile.default <- tempfile( "cohortFile.default" )
   writeLines( cohortArray.default, cohortFile.default )

   describe( "default behavior", {
      it( 'Generates the expected data frame for a simple cohort file', {
         df <- loadCohortDefinition( cohortFile.default )

         expect_is(df, 'data.frame')
         wantNames <- c("sample", "exonExpressionFile", "extra" )
         got <- names(df)
         expect_equal(got, wantNames)

         rowDat <- c("1", "2", "4")
         wantSamples <- paste0( "bob-", rowDat )
         expect_equal( df$sample, wantSamples )

         wantExonExpressionFiles <- paste0( "bob-", rowDat, "/expression.txt" )
         expect_equal( df$exonExpressionFile, wantExonExpressionFiles )

         wantExtra <- as.numeric(paste0( rowDat, rowDat ))
         expect_equal( df$extra, wantExtra )
      })
   })

   describe( "Options", {
      describe( "file=", {
         it ("Errors if file not found", {
            missingFile <- "NoSuchFile"
            expect_false( file.exists( missingFile ))
            wantErrorRE <- "Can't find the cohort definition file: \".*NoSuchFile.*\""
            expect_error( loadCohortDefinition( missingFile ), wantErrorRE )
         })
      })

      describe( "samples=", {
         describe( "Selection of rows specified", {
            it( "Select all rows with explicit default settings", {
               df <- loadCohortDefinition( file= cohortFile.default, samples= NULL )
               wantSamples <- c( "bob-1", "bob-2", "bob-4" )
               expect_equal( df$sample, wantSamples )
            })
            it( "Selecting correct rows if only listing some samples", {
               wantSamples <- "bob-1"
               df <- loadCohortDefinition( file= cohortFile.default, samples= wantSamples )
               expect_equal( df$sample, wantSamples )

               wantSamples <- c( "bob-2", "bob-4" )
               df <- loadCohortDefinition( file= cohortFile.default, samples= wantSamples )
               expect_equal( df$sample, wantSamples )
            })
            it( "Select all rows if listing all samples", {
               wantSamples <- c( "bob-1", "bob-2", "bob-4" )
               df <- loadCohortDefinition( file= cohortFile.default, samples= wantSamples )
               expect_equal( df$sample, wantSamples )
            })
         })
         describe( "Handling of badly specified sample list", {
            it('Selects specified rows when some samples case sensitively different', {
               badSampleList <- c( "bob-1", "BOB-2", "bob-4" )
               df <- loadCohortDefinition( file= cohortFile.default, samples= badSampleList )
               wantSamples <- c( "bob-1", "bob-4" )
               expect_equal( df$sample, wantSamples )

            })
            it('Selects all rows when some additional samples are present', {
               extraSamplesList <- c( "bob-1", "bob-2", "bob-3", "bob-4" )
               df <- loadCohortDefinition( file= cohortFile.default, samples= extraSamplesList )
               wantSamples <- c( "bob-1", "bob-2", "bob-4" )
               expect_equal( df$sample, wantSamples )
            })
            it('Warns if sample list contains unknown samples', {
               wantErrorRE <- "Ignoring missing samples specified in the \"samples\" parameter list:.*bob-3"
               extraSamplesList <- c( "bob-1", "bob-2", "bob-3", "bob-4" )
               expect_warning( loadCohortDefinition( cohortFile.default, samples= extraSamplesList ),
                             wantErrorRE )
            })
            it('Warning with unknown samples lists all unknown entries', {
               wantErrorRE <- "NotASample.*bob-3"
               extraSamplesList <- c( "bob-1", "NotASample", "bob-2", "bob-3"  )
               expect_warning( df <- loadCohortDefinition( cohortFile.default, samples= extraSamplesList ),
                               wantErrorRE )
               wantSamples <- c( "bob-1", "bob-2" )
               expect_equal( df$sample, wantSamples )
            })
            it('Warns if duplicates are present in the sample selection list', {
               wantErrorRE <- "Duplicates in \"samples\" parameter filter list ignored"
               extraSamplesList <- c( "bob-1", "bob-1", "bob-2"  )
               expect_warning( df <- loadCohortDefinition( cohortFile.default, samples= extraSamplesList ),
                               wantErrorRE )
               wantSamples <- c( "bob-1", "bob-2" )
               expect_equal( df$sample, wantSamples )
            })
            it('All warnings needed are issued', {
               wantErrorRE <- "Duplicates in \"samples\" parameter filter list ignored"
               extraSamplesList <- c( "bob-1", "NotASample", "bob-1", "bob-2"  )
               expect_warning( df <- loadCohortDefinition( cohortFile.default, samples= extraSamplesList ),
                               wantErrorRE )
               wantSamples <- c( "bob-1", "bob-2" )
               expect_equal( df$sample, wantSamples )
               wantErrorRE <- "Ignoring missing samples specified in the \"samples\" parameter list:.*NotASample"
               expect_warning( df <- loadCohortDefinition( cohortFile.default, samples= extraSamplesList ),
                               wantErrorRE )
               wantSamples <- c( "bob-1", "bob-2" )
               expect_equal( df$sample, wantSamples )
            })
         })
      })

      describe( "dataRoot=", {
         it("Does nothing with explicit defaults", {
            df <- loadCohortDefinition( cohortFile.default, dataRoot= NULL )
            rowDat <- c("1", "2", "4")
            wantExonExpressionFiles <- paste0( "bob-", rowDat, "/expression.txt" )
            expect_equal( df$exonExpressionFile, wantExonExpressionFiles )

            df <- loadCohortDefinition( cohortFile.default, dataRoot= NA )
            rowDat <- c("1", "2", "4")
            wantExonExpressionFiles <- paste0( "bob-", rowDat, "/expression.txt" )
            expect_equal( df$exonExpressionFile, wantExonExpressionFiles )

            df <- loadCohortDefinition( cohortFile.default, dataRoot= "" )
            rowDat <- c("1", "2", "4")
            wantExonExpressionFiles <- paste0( "bob-", rowDat, "/expression.txt" )
            expect_equal( df$exonExpressionFile, wantExonExpressionFiles )
         })
         it("Prefixis dataRoot to expressionFilePath column", {
            df <- loadCohortDefinition( cohortFile.default, dataRoot= "/some/root" )
            rowDat <- c("1", "2", "4")
            wantExonExpressionFiles <- paste0( "/some/root/bob-", rowDat, "/expression.txt" )
            expect_equal( df$exonExpressionFile, wantExonExpressionFiles )
         })
      })

      describe( "comment.char=", {
         it ("Drops comments if explicitly use comment.char= '#'", {
            df <- loadCohortDefinition( cohortFile.default, comment.char= '#' )
            rowDat <- c("1", "2", "4")
            wantSamples <- paste0( "bob-", rowDat )
            expect_equal( df$sample, wantSamples )
         })
         it ("Can change comment.char= ';'", {
            cohortArray.semi <- c(
               "; This is a fake cohort.",
               "; It is used only for testing",
               paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
               paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
               paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
               paste( ";bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
               paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
            )
            cohortFile.semi <- tempfile( "cohortArray.percent" )
            writeLines( cohortArray.semi, cohortFile.semi )
            df <- loadCohortDefinition( cohortFile.semi, comment.char= ';' )
            rowDat <- c("1", "2", "4")
            wantSamples <- paste0( "bob-", rowDat )
            expect_equal( df$sample, wantSamples )

         })
         it ("Reads comment lines if set comment.char= ''", {
               cohortArray.noComment <- c(
                  paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
                  paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
                  paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
                  paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
                  paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
               )
               cohortFile.noComment <- tempfile( "cohortArray.noComment" )
               writeLines( cohortArray.noComment, cohortFile.noComment )
               df <- loadCohortDefinition( cohortFile.noComment, comment.char= '' )
               wantSamples <- c("bob-1", "bob-2","#bob-3", "bob-4")
               expect_equal( df$sample, wantSamples )
         })
         it ("Reads comment lines if comment.char= '#' has leading whitespace", {
            cohortArray.oops <- c(
               paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
               paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
               paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
               paste( " #bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
               paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
            )
            cohortFile.oops <- tempfile( "cohortArray.oops" )
            writeLines( cohortArray.oops, cohortFile.oops )
            df <- loadCohortDefinition( cohortFile.oops, comment.char= '' )
            wantSamples <- c("bob-1", "bob-2", " #bob-3", "bob-4")
            expect_equal( df$sample, wantSamples )
         })
      })
   })

   describe( "Validating input data file values", {
      it ("Errors without sample column", {
         cohortArray.noSample <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "SAMPLE", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
         )
         cohortFile.noSample <- tempfile( "cohortFile.noSample" )
         writeLines( cohortArray.noSample, cohortFile.noSample )
         wantErrorRE <- "The cohort file has no sample column."
         expect_error( df <- loadCohortDefinition( cohortFile.noSample ), wantErrorRE )
      })
      it ("Errors without an exonExpressionFile column", {
         cohortArray.noExpression <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFiles", "extra", sep= "\t"),
            paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
         )
         cohortFile.noExpression <- tempfile( "cohortFile.noExpression" )
         writeLines( cohortArray.noExpression, cohortFile.noExpression )
         wantErrorRE <- "The cohort file has no exonExpressionFile column."
         expect_error( df <- loadCohortDefinition( cohortFile.noExpression ), wantErrorRE )
      })
      it ("Errors without strings in sample column", {
         cohortArray.sampleInt <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  "2", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "4", "bob-4/expression.txt",  "44", sep= "\t")
         )
         cohortFile.sampleInt <- tempfile( "cohortFile.sampleInt" )
         writeLines( cohortArray.sampleInt, cohortFile.sampleInt )
         wantErrorRE <- "The sample and exon expression file columns must contain text"
         expect_error( df <- loadCohortDefinition( cohortFile.sampleInt ), wantErrorRE )
      })
      it ("Errors without strings in exonExpressionFile column", {
         cohortArray.expressionInt <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "1",  "11", sep= "\t"),
            paste(  "bob-2", "2",  "22", sep= "\t"),
            paste( "bob-3", "3",  "33", sep= "\t"),
            paste(  "#bob-4", "4",  "44", sep= "\t")
         )
         cohortFile.expressionInt <- tempfile( "cohortFile.expressionInt" )
         writeLines( cohortArray.expressionInt, cohortFile.expressionInt )
         wantErrorRE <- "The sample and exon expression file columns must contain text"
         expect_error( df <- loadCohortDefinition( cohortFile.expressionInt ), wantErrorRE )
      })
      it ("Errors with duplicate sample values", {
         cohortArray.dupSample <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-1", "bob-4/expression.txt",  "44", sep= "\t")
         )
         cohortFile.dupSample <- tempfile( "cohortFile.dupSample" )
         writeLines( cohortArray.dupSample, cohortFile.dupSample )
         wantErrorRE <- "The cohort file can not contains duplicate samples"
         expect_error( df <- loadCohortDefinition( cohortFile.dupSample ), wantErrorRE )
      })
      it ("Errors with duplicate exon expression file values", {
         cohortArray.dupExpression <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-4", "bob-2/expression.txt",  "44", sep= "\t")
         )
         cohortFile.dupExpression <- tempfile( "cohortFile.dupExpression" )
         writeLines( cohortArray.dupExpression, cohortFile.dupExpression )
         wantErrorRE <- "The cohort file can not contain duplicate cohort expression files"
         expect_error( df <- loadCohortDefinition( cohortFile.dupExpression ), wantErrorRE )
      })
      it ("Errors with empty string in sample column", {
         cohortArray.emptySample <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  "", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
         )
         cohortFile.emptySample <- tempfile( "cohortFile.emptySample" )
         writeLines( cohortArray.emptySample, cohortFile.emptySample )
         wantErrorRE <- "All sample and exon expression file entries must be non-empty text"
         expect_error( df <- loadCohortDefinition( cohortFile.emptySample ), wantErrorRE )
      })
      it ("Errors with NA string in sample column", {
         cohortArray.naSample <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  NA, "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
         )
         cohortFile.naSample <- tempfile( "cohortFile.naSample" )
         writeLines( cohortArray.naSample, cohortFile.naSample )
         wantErrorRE <- "All sample and exon expression file entries must be non-empty text"
         expect_error( df <- loadCohortDefinition( cohortFile.naSample ), wantErrorRE )
      })
      it ("Errors with empty string in exonExpressionFile column", {
         cohortArray.emptyExpression <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "",  "11", sep= "\t"),
            paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-4", "bob-4/expression.txt",  "44", sep= "\t")
         )
         cohortFile.emptyExpression <- tempfile( "cohortFile.emptyExpression" )
         writeLines( cohortArray.emptyExpression, cohortFile.emptyExpression )
         wantErrorRE <- "All sample and exon expression file entries must be non-empty text"
         expect_error( df <- loadCohortDefinition( cohortFile.emptyExpression ), wantErrorRE )
      })
      it ("Errors with na string in exonExpressionFile column", {
         cohortArray.naExpression <- c(
            "# This is a fake cohort.",
            "# It is used only for testing",
            paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
            paste(  "bob-1", "bob-1/expression.txt",  "11", sep= "\t"),
            paste(  "bob-2", "bob-2/expression.txt",  "22", sep= "\t"),
            paste( "#bob-3", "bob-3/expression.txt",  "33", sep= "\t"),
            paste(  "bob-4", NA,  "44", sep= "\t")
         )
         cohortFile.naExpression <- tempfile( "cohortFile.naExpression" )
         writeLines( cohortArray.naExpression, cohortFile.naExpression )
         wantErrorRE <- "All sample and exon expression file entries must be non-empty text"
         expect_error( df <- loadCohortDefinition( cohortFile.naExpression ), wantErrorRE )
      })
   })

})

describe( "loadExonExpressionFile()", {

   describe( "A default run", {
      it( "returns a data.frame", {

      })
      it( "returns the correct columns", {

      })
      it( "returns the correct column contents for column '...'", {

      })
   })

   describe( "Options", {
      describe( "file=", {

      })
      describe( "type=", {

      })
   })
})