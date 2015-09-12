context("Unit tests for functions in the dataRangeReductionUtils.R file.")

describe( "mapLabels()", {

   it( 'calculates correct values for a data vector.', {
      ends = c(-Inf, -1, -0.01, 0, 0.01, 1, Inf)
      labels = c('---', '--', '-', '+', '++', '+++')

      data = c(-Inf, -3e5, -1.001, -1, -0.5, -0.01, -0.001, -0, 0)
      got = mapLabels(data, ends, labels)
      want = c('---', '---', '---', '---', '--', '--', '-', '-', '-')
      expect_equal(got, want)

      data = c(Inf, 3e5, 1.001, 1, 0.5, 0.01, 0.001, +0, 0)
      got = mapLabels(data, ends, labels)
      want = c('+++', '+++', '+++', '++', '++', '+', '+', '-', '-')
      expect_equal(got, want)
   })

   it( 'handles single value data.', {
      ends = c(-Inf, -1, -0.01, 0, 0.01, 1, Inf)
      labels = c('---', '--', '-', '+', '++', '+++')

      got = mapLabels(7, ends, labels)
      want = '+++'
      expect_equal(got, want)
   })

   it( 'handles missing/out of range values.', {
      ends = c(-10, -1, 0, 1, 10)
      labels = c('--', '-', '+', '++')
      data = c(NA_integer_, 2, 23)

      got = mapLabels(data, ends, labels)
      want = c(NA_character_, '++', NA_character_)
      expect_equal(got, want)
   })

   it( 'sorts input bin values.', {
      ends = c(Inf,-Inf,1,-1)
      labels = c('-', '0', '+' )
      data = c(-Inf,-5,-1,-0.5,0,0.5,1,5,Inf)

      got = mapLabels(data, ends, labels)
      want = c('-','-','-','0','0','0','0','+','+')
      expect_equal(got, want)
   })

   it( 'errors if label count and bin number are not the same', {
      wantErrorRE <- 'Too many labels specified'
   })
})
#
# mapColors()
#

test_that( 'mapColors; adds or returns correct column for default pallette', {
   x <- c(-100, -1, 0, 1, 1000)
   binEnds = c(-Inf,-10,10,Inf)
   colSet <- brewer.pal(3, 'RdBu')

   # Forward
   got <- mapColors( x, binEnds=binEnds )
   want <- c(colSet[1], colSet[2], colSet[2], colSet[2], colSet[3])
   expect_equal(got, want)

   # Reverse
   got <- mapColors( x, binEnds=binEnds, reverse=TRUE )
   want <- c(colSet[3], colSet[2], colSet[2], colSet[2], colSet[1])
   expect_equal(got, want)

   colSet <- c('red', 'violet', 'blue')
   got <- mapColors( x, binEnds=binEnds, colors= colSet )
   want <- c(colSet[1], colSet[2], colSet[2], colSet[2], colSet[3])
   expect_equal(got, want)

})

