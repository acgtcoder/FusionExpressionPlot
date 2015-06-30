context("Testing the test infrastructure")

test_that( "Testing works?", {
   succeed( "Testing with testthat" )
})

test_that( "Test data is available?", {
   expect_true(  dir.exists(  'data' ))
   expect_true(  file.exists( 'data/mock.gaf' ))
   expect_true(  file.exists( 'data/expect.geneModels' ))
   expect_true(  file.exists( 'data/expect.transcriptModels' ))
})