context("Unit tests for functions in the dataRangeReductionUtils.R file.")

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
