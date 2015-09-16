context("cohortNormalizedExpression.R")

describe( "loadCohortDefinition()", {

   cohortArray.default <- c(
      "# This is a fake cohort.",
      "# It is used only for testing",
      paste( "sample", "exonExpressionFile", "extra", sep= "\t"),
      paste( "bob-1", "bob-1/expression.txt", "11", sep= "\t"),
      paste( "bob-2", "bob-2/expression.txt", "22", sep= "\t"),
      paste( "#bob-3", "bob-3/expression.txt", "33", sep= "\t"),
      paste( "bob-4", "bob-4/expression.txt", "44", sep= "\t")
   )
   cohortFile.default <- tempfile( "cohortFile.default" )
   writeLines( cohortArray.default, cohortFile.default )

   describe( "default behavior", {
      it( 'Generates the expected data frame for a simple cohort file', {
         df <- loadCohortDefinition( cohortFile.default )

         expect_is(df, 'data.frame')
         wantNames <- c("sample", "exonExpressionFile", "extra")
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
            wantErrorRE <- "Can't find the cohort definition file: \"NoSuchFile\""
         })
      })

      describe( "samples=", {
         describe( "In context of default sample column", {
            describe( "Selection of rows specified", {
               it( "Select all rows with explicit default settings", {
               })
               it( "Selecting correct rows if only listing some samples", {
               })
               it( "Select all rows if listing all samples", {
               })
            })
            describe( "Handling of badly specified sample list", {
               it('Selects specified rows when some samples case sensitively different', {
               })
               it('Selects all rows when some additional samples are present', {
               })
               it('Warns if sample list contains unknown samples', {
                  wantErrorRE <- "Ignoring missing samples specified in the \"samples\" parameter list:"
               })
               it('Warning with unknown samples lists all unknown entries', {
                  wantErrorRE <- "Bad sample name 1"
                  wantErrorRE <- "Bad sample name 2"
               })
               it('Warns if duplicates are present in the sample selection list', {
                  wantErrorRE <- "Duplicates in \"samples\" parameter filter list ignored"
               })
            })
         })
         describe( "In context of specified sample column", {
            describe( "Selection of rows specified", {
               it( "Select all rows with explicit default settings", {
               })
               it( "Selecting correct rows if only listing some samples", {
               })
               it( "Select all rows if listing all samples", {
               })
            })
            describe( "Handling of badly specified sample list", {
               it('Selects specified rows when some samples case sensitively different', {
               })
               it('Selects all rows when some additional samples are present', {
               })
               it('Warns if sample list contains unknown samples', {
                  wantErrorRE <- "Ignoring missing samples specified in the \"samples\" parameter list:"
               })
               it('Warning with unknown samples lists all unknown entries', {
                  wantErrorRE <- "Bad sample name 1"
                  wantErrorRE <- "Bad sample name 2"
               })
               it('Warns if duplicates are present in the sample selection list', {
                  wantErrorRE <- "Duplicates in \"samples\" parameter filter list ignored"
               })
            })
         })
      })

      describe( "columns=", {
         it( "Works with explicit default", {

         })
         it( "Can choose different columns", {

         })
         describe( "Handling invalid columns parameter", {
            wantErrorRE <- "Columns must use strings for names."
            wantErrorRE <- "Must specify both the sample and exonExpressionColumn headers"
            wantErrorRE <- "Columns must specify non-empty strings for column names"
            wantErrorRE <- "May not specify the same column for sample and exonExpressionFile data"
         })
      })

      describe( "dataRoot=", {
         describe( "In context of default expressionFilePath column", {
            it("Does nothing with explicit defaults", {
            })
            it("Prefixis dataRoot to expressionFilePath column", {
            })
         })
         describe( "In context of specified expressionFilePath column", {
            it("Does nothing with explicit defaults", {
            })
            it("Prefixis dataRoot to expressionFilePath column", {
            })
         })
      })

      describe( "comment.char=", {
         it ("Drops comments if explicitly use comment.char= '#'", {
         })
         it ("Can change comment.char= ';'", {
         })
         it ("Reads comment lines if set comment.char= ''", {
         })
         it ("Reads comment lines if comment.char= '#' has leading whitespace", {
         })
      })
   })

   describe( "Validating data frame values", {
      wantErrorRE <- "The cohort file has no sample column \"sample\""
      wantErrorRE <- "The cohort file has no exon expression file column \"exonExpressionColumn\""
      wantErrorRE <- "The sample and exon expression file columns must contain strings"
      wantErrorRE <- "The cohort file can not contains duplicate samples"
      wantErrorRE <- "The cohort file can not contain duplicate cohort expression files"
      wantErrorRE <- "All sample and exon expression file entries must be non-empty strings"
   })

})