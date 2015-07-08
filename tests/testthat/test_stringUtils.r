context("Unit testing functions in the stringUtils.R file.")

#
# regexprNamedMatches
#

test_that( "regexprNamedMatches returns correct results", {

   # Normal matching
   regExp <- "(?<key>.+?)\\s*=\\s*(?<value>.+)"
   data <- c('name = Stuart R. Jefferys', 'email=srj@unc.edu')

   matchResults <- regexpr(regExp, data, perl= TRUE)
   got<-regexprNamedMatches(matchResults, data)
   want<-do.call(cbind, list(key=c('name', 'email'), value=c('Stuart R. Jefferys', 'srj@unc.edu')))
   expect_equal(got, want)

   # If only one string
   regExp <- "(?<key>.+?)\\s*=\\s*(?<value>.+)"
   data <- c('name = Stuart R. Jefferys')

   matchResults <- regexpr(regExp, data, perl= TRUE)
   got<-regexprNamedMatches(matchResults, data)
   want<-do.call(cbind, list(key=c('name'), value=c('Stuart R. Jefferys')))
   expect_equal(got, want)

   # If only one capture group
   regExp <- ".*?(?<num>\\d+).*"
   data <- c('Embedded 123 number', '0234')

   matchResults <- regexpr(regExp, data, perl= TRUE)
   got<-regexprNamedMatches(matchResults, data)
   want<-do.call(cbind, list(num=c('123', '0234')))
   expect_equal(got,want)

   # If doesn't match?
   regExp <- "(?<key>.+?)\\s*=\\s*(?<value>\\d+)"
   data <- c('id = 0123', 'email=srj@unc.edu')

   # Want ''
   matchResults <- regexpr(regExp, data, perl= TRUE)
   got<-regexprNamedMatches(matchResults, data)
   want<-do.call(cbind, list(key=c('id', ''), value=c('0123', '')))
   expect_equal(got, want)

   # Want NAs
   matchResults <- regexpr(regExp, data, perl= TRUE)
   got<-regexprNamedMatches(matchResults, data, use.na=TRUE)
   want<-do.call(cbind, list(key=c('id', NA), value=c('0123', NA)))
   expect_equal(got, want)

   # If matches with substring empty?
   regExp <- "(?<key>.+?)\\s*=\\s*(?<value>\\d*)"
   data <- c('id = 0123', 'email=srj@unc.edu')

   # Want ''
   matchResults <- regexpr(regExp, data, perl= TRUE)
   got<-regexprNamedMatches(matchResults, data)
   want<-do.call(cbind, list(key=c('id', 'email'), value=c('0123', '')))
   expect_equal(got, want)

   # Still want '', not NAs
   matchResults <- regexpr(regExp, data, perl= TRUE)
   got<-regexprNamedMatches(matchResults, data, use.na=TRUE)
   want<-do.call(cbind, list(key=c('id', 'email'), value=c('0123', '')))
   expect_equal(got, want)

})