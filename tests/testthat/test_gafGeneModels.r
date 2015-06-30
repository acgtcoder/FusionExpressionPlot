context("gafGeneModels.R file.")


test_that( "extractGeneModels(); Generates correct default data file?", {
   makesFile <- 'data/mock.gaf.geneModels'
   extractGeneModels( 'data/mock.gaf' )
   expect_true( file.exists( makesFile ))
   expect_equivalent( tools::md5sum( 'data/expect.geneModels' ),
                      tools::md5sum( makesFile )
   )
   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }
})

test_that( "extractGeneModels(); Correct error handling with default parameters?", {
   makesFile <- 'data/mock.gaf.geneModels'
   expect_error( extractGeneModels('noSuchFile.gaf'),
                 "Can't find the specified gaf file: 'noSuchFile.gaf'"
   )
   expect_false( file.exists( makesFile ))

   expect_error( extractGeneModels('my.gaf;ls'),
                 'Unsafe character in gaf file name!'
   )
   expect_false( file.exists( makesFile ))

   extractGeneModels( 'data/mock.gaf' )
   expect_true( file.exists( makesFile ))
   expect_error( extractGeneModels( 'data/mock.gaf' ),
                 "Output file already exists; use force=TRUE to overwrite: 'data/mock.gaf.geneModels'"
   )

   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }

})

test_that( 'extractGeneModels() can specify output file?', {
   makesFile <- 'data/deleteMe.geneModels'
   extractGeneModels( 'data/mock.gaf', outFile= makesFile )
   expect_true( file.exists( makesFile ))
   expect_equivalent( tools::md5sum( 'data/expect.geneModels' ),
                      tools::md5sum( makesFile )
   )
   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }
})

test_that( 'extractGeneModels() error handling with output file?', {
   makesFile <- 'data/deleteMe.geneModels'
   expect_error( extractGeneModels('noSuchFile.gaf', outFile= makesFile ),
                 "Can't find the specified gaf file: 'noSuchFile.gaf'"
   )
   expect_false( file.exists( makesFile ))

   expect_error( extractGeneModels('my.gaf;ls', outFile= makesFile ),
                 'Unsafe character in gaf file name!'
   )
   expect_false( file.exists( makesFile ))

   expect_error( extractGeneModels('data/mock.gaf', outFile= 'data/deleteMe;ls'),
                 'Unsafe character in output geneModel file name!'
   )
   expect_false( file.exists( 'data/deleteMe;ls' ))

   extractGeneModels( 'data/mock.gaf', outFile= makesFile )
   expect_true( file.exists( makesFile ))
   expect_error( extractGeneModels( 'data/mock.gaf', outFile= makesFile ),
                 "Output file already exists; use force=TRUE to overwrite: 'data/deleteMe.geneModels'"
   )

   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }

})

test_that( 'extractGeneModels() can overwrite output file?', {
   makesFile <- 'data/mock.gaf.geneModels'
   extractGeneModels( 'data/mock.gaf' )
   expect_true( file.exists( makesFile ))
   expect_warning( extractGeneModels( 'data/mock.gaf', force= TRUE ),
                   "Forcing overwrite of output file 'data/mock.gaf.geneModels'"
   )
   expect_equivalent( tools::md5sum( 'data/expect.geneModels' ),
                      tools::md5sum( makesFile )
   )
   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }
})

test_that( "extractTranscriptModels(); Generates correct default data file?", {
   makesFile <- 'data/mock.gaf.transcriptModels'
   extractTranscriptModels( 'data/mock.gaf' )
   expect_true( file.exists( makesFile ))
   expect_equivalent( tools::md5sum( 'data/expect.transcriptModels' ),
                      tools::md5sum( makesFile )
   )
   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }

})

test_that( "extractTranscriptModels(); Correct error handling with default parameters?", {

   makesFile <- 'data/mock.gaf.transcriptModels'
   expect_error( extractTranscriptModels('noSuchFile.gaf'),
                 "Can't find the specified gaf file: 'noSuchFile.gaf'"
   )
   expect_false( file.exists( makesFile ))

   expect_error( extractTranscriptModels('my.gaf;ls'),
                 'Unsafe character in gaf file name!'
   )
   expect_false( file.exists( makesFile ))

   extractTranscriptModels( 'data/mock.gaf' )
   expect_true( file.exists( makesFile ))
   expect_error( extractTranscriptModels( 'data/mock.gaf' ),
                 "Output file already exists; use force=TRUE to overwrite: 'data/mock.gaf.transcriptModels'"
   )

   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }

})

test_that( 'extractTrascriptModels() can specify output file?', {
   makesFile <- 'data/deleteMe.transcriptModels'
   extractTranscriptModels( 'data/mock.gaf', outFile= makesFile )
   expect_true( file.exists( makesFile ))
   expect_equivalent( tools::md5sum( 'data/expect.transcriptModels' ),
                      tools::md5sum( makesFile )
   )
   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }

})

test_that( 'extractTranscriptModels() error handling with output file?', {
   makesFile <- 'data/deleteMe.transcriptModels'
   expect_error( extractTranscriptModels('noSuchFile.gaf', outFile= makesFile ),
                 "Can't find the specified gaf file: 'noSuchFile.gaf'"
   )
   expect_false( file.exists( makesFile ))

   expect_error( extractTranscriptModels('my.gaf;ls', outFile= makesFile ),
                 'Unsafe character in gaf file name!'
   )
   expect_false( file.exists( makesFile ))

   expect_error( extractTranscriptModels('data/mock.gaf', outFile= 'data/deleteMe;ls'),
                 'Unsafe character in output transcriptModel file name!'
   )
   expect_false( file.exists( 'data/deleteMe;ls' ))

   extractTranscriptModels( 'data/mock.gaf', outFile= makesFile )
   expect_true( file.exists( makesFile ))
   expect_error( extractTranscriptModels( 'data/mock.gaf', outFile= makesFile ),
                 "Output file already exists; use force=TRUE to overwrite: 'data/deleteMe.transcriptModels'"
   )

   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }

})

test_that( 'extractTranscriptModels() can overwrite output file?', {
   makesFile <- 'data/mock.gaf.transcriptModels'
   extractTranscriptModels( 'data/mock.gaf' )
   expect_true( file.exists( makesFile ))
   expect_warning( extractTranscriptModels( 'data/mock.gaf', force= TRUE ),
                   "Forcing overwrite of output file 'data/mock.gaf.transcriptModels'"
   )
   expect_equivalent( tools::md5sum( 'data/expect.transcriptModels' ),
                      tools::md5sum( makesFile )
   )
   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }
})

test_that( 'extractTranscriptModels() returns correct value', {
   makesFile <- 'data/mock.gaf.transcriptModels'
   expect_equal( extractTranscriptModels( 'data/mock.gaf' ), 1)
   # Clean up after ourselves
   if ( file.exists( makesFile )) {
      unlink( makesFile );
   }
})