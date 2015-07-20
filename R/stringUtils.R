# This file provides utility functions for working with strings.

#' Extract matched named substrings
#'
#' Extracts the substrings matched by named capture groups from the provided
#' match result output from \code{\link[base]{regexpr}} (with \code{perl=
#' TRUE}). Will return a matrix of strings with one column for each named
#' capture group and one row for each string in the vector matched against. By
#' default will return empty strings if match fails, but can be set to return
#' NAs if desired.
#'
#' This is intended for use with \code{\link[base]{regexpr}} to parse a string
#' and extract substrings via named capture groups, similar to how
#' \code{\link[base]{regmatches}} is used. If only one string is matched
#' against, then returned matrix will have one row only.
#'
#' Note that regExp with multiple capture groups will need to use greedy and
#' non-greedy matching carefully if the capture groups are to work correctly and
#' not interfering with each other or not-capture components of the regExp.
#'
#' @param matchResults The results of a match performed using
#'   \code{\link[base]{regexpr}(regExp, matchText, perl= TRUE)} where
#'   \code{regExp} has named capture groups like \code{(?<theName>...)}.
#'
#' @param matchText The text originally matched against, a vector of strings.
#'
#' @param use.na By default returns empty strings for all capture groups if the
#'   regExp fails to match. That can not be distinguished from a match with all
#'   capture groups matching nothing, e.g. \code{(?<num>\\d*)}. Setting this
#'   \code{TRUE} causes a failing match to return all NA values instead.
#'
#' @return A matrix with one column for each named capture group (with matching
#'   column name) and one row for each string in the text vector matched to. The
#'   value of each cell is the text matched by the named capture group. If any
#'   capture group does not match, all returned strings are empty for that text
#'   vector element (row), or \code{NA} if \code{use.na= TRUE}
#'
#' @examples
#' regExp <- "(?<key>.+?)\\s*=\\s*(?<value>.+)"
#' data <- c('name = Stuart R. Jefferys', 'email=srj@@unc.edu')
#' matchResults <- regexpr(regExp, data, perl= TRUE)
#' regexprNamedMatches(matchResults, data)
#' #=>      key     value
#' #=> [1,] "name"  "Stuart R. Jefferys"
#' #=> [2,] "email" "srj@@unc.edu"
#'
#' @export
regexprNamedMatches <- function( matchResults, matchText, use.na=FALSE ) {

   captureNames <- attr(matchResults,'capture.names')
   nrows <- length(matchText)
   ncols <- length(captureNames)
   retMat <- matrix(character(nrows*ncols), nrow = nrows, ncol = ncols, dimnames=list(rep(NULL,nrows),captureNames))
   captureStarts <- attr(matchResults,'capture.start')
   captureLengths <- captureStarts + attr(matchResults,'capture.length') - 1
   for (captureName in captureNames) {
      retMat[,captureName] = substr(matchText,captureStarts[,captureName], captureLengths[,captureName])
   }

   # Simple but possibly inefficient to just reset values afterwards.
   if (use.na) {
      for (row in 1:nrows) {
         if (matchResults[row] == -1) {
            retMat[row,] <- rep(NA, ncols)
         }
      }
   }
   return(retMat)
}

#' Evaluate and fill string templates
#'
#' Given a vector of strings containing {{variables}} and returns a copy
#' replacing the templated fields with the value of the specified variables.
#' Variables must be defined in the calling environment (or the one passed in),
#' or an error will occur. If as.R is set TRUE, then any {{R code}} can be used
#' and it will be evaluated to obtain a return value. That is, of course,
#' dangerous if you don't trust the source of your template. This code is
#' executed as if run in new function called in the current environment (or the
#' one passed in).
#'
#' @param x Vector of strings with fields containing variables or code to be
#'   interpolated.
#'
#' @param delim Vector of two string, the first used to signal the start of a
#' template section and the second used to signal the end. These may not be
#' the same, nor have one embedded in the other/ By default the delimiers are
#' c( "{{", "}}" )
#'
#' @param as.R Set TRUE to allow full R code evaluation. By default is FALSE and
#'   only allows variable substitution. Setting this true is a security risk
#'   when you don't trust the provider of the template text as much as you trust
#'   the person who provided your R code, so it generates a warning.
#'
#' @param envir The execution environment to be used. Can be used to pass in the
#'   an environment in which variables are defined for use in interpolation. If
#'   not specified, then by default this will be a new environment whose parent
#'   is the callers environment, as returned by parent.frame(). Variables
#'   visible in the calling function (or set there) will be avialble for use in
#'   the template. Note that although R code will normally only set or change
#'   variables in this frame when evaluated, it can set or change variables at
#'   any level, hence malicous or careless as.R evaluated templates can leak or
#'   interfere with other R variables in your code (or indeed in any other
#'   package or even system code). With great power comes great responsibility.
#'
#' @return A copy of the original vector of strings, but with variable names
#'   replaced with their values, or with the result of evaluating the
#'   interpolated string as R code. Note that everything is returned as a
#'   string, so '{1+1}' is returned as '2'.
#'
#' @examples
#' # Template is single text element (could be multi-line)
#' templateText <- "Dear {{name}}: Please call me at {{phone}}."
#' name<-"John Doe"
#' phone<-"555-555-5555"
#' templateFill(templateText)
#' #=> [1] "Dear John Doe: Please call me at 555-555-5555."
#'
#' # Delimiters can be changed
#' templateText <- "Dear -<[name]>-: Please contact me at -<[email]>-."
#' name<-"John"
#' email<-"the.bobs@@layoffs.com"
#' templateFill(templateText, delim=c('-<[', ']>-'))
#' #=> [1] "Dear John: Please contact me at the.bobs@@layoffs.com."
#'
#' # Multiple text elements (each could be multi line)
#' templateText <- c("ID: {{id}}", "Item: {{name}}", "Description: {{desc}}")
#' id <- "0001-12"
#' name <- "widget"
#' desc <- "Widget to foo the bar."
#' templateFill(templateText)
#' #=> [1] "ID: 0001-12"
#' #=> [2] "Item: widget"
#' #=> [3] "Description: Widget to foo the bar."
#'
#' # Evaluating R code
#' x <- 21
#' y <- 'Helloooo'
#' templateText <- c(
#'     "Simple: {{1 + 1}}",
#'     "Variables are accessible: {{x *2}}",
#'     "Complex: {{ echo <- function(x) { paste(x,x,sep='...') }; echo(y) }}",
#'     "Code environment is shared: {{ echo( 'Goodbyyyy' ) }}"
#' )
#' templateFill(templateText, as.R=TRUE)
#' #=> [1] "Simple: 2"
#' #=> [2] "Variables are accessible: 42"
#' #=> [3] "Complex: Helloooo...Helloooo"
#' #=> [4] "Code environment is shared: Goodbyyyy...Goodbyyyy"
#' #=> Warning message:
#' #=> In templateFill(templateText, as.R = TRUE) :
#' #=>    Potential security risk: templateFill() is evaluating user-provided
#' #=>    R code If you trust where the template is coming from, you can
#' #=>    suppress this message with suppressWarnings().
#' @export
templateFill <- function( x,
                          delim = c( '{{', '}}' ),
                          as.R = FALSE, envir = new.env( parent= parent.frame() )
) {
   if (length(delim) != 2) {
      stop("delim= must have exactly two elements.")
   }
   if (delim[1] == delim[2]) { stop("delim= must have different open and close elements") }
   if (    grepl(delim[1], delim[2], fixed= TRUE)
        || grepl(delim[2], delim[1], fixed= TRUE)
   ) {
      stop("Can't have one of the delimiters embeded in the other.")
   }
   if (as.R) {
      warning( "Potential security risk:",
               " templateFill() is evaluating user-provided R code",
               " If you trust where the template is coming from,",
               " you can suppress this message with suppressWarnings().")
   }

   # Find delimiter positions
   starts <- gregexpr(delim[1], x, fixed= TRUE)
   ends <- gregexpr(delim[2], x, fixed= TRUE)

   # Pre-allocate the returned vector of strings
   retVal <- character(length(x))

   # Process each string in the input vector (possibly 0)
   for (stringNum in 1:length(x)) {
      # Any string without BOTH delimiters is just returned as is
      if (starts[[stringNum]] == -1 || ends[[stringNum]] == -1) {
         retVal[stringNum] <- x[stringNum]
         next
      }
      # If any string has both delimiters, but has a mismatched number of open
      # and closed delimiters, fail for the whole thing.
      if (length(starts[[stringNum]]) > length(ends[[stringNum]])) {
         stop("Too many ", delim[1], " found in template text element ", stringNum, ".",
              " Probably missing one or more ", delim[2], ".")
      }
      else if (length(starts[[stringNum]]) < length(ends[[stringNum]])) {
         stop("Too many ", delim[2], " found in template text element ", stringNum, ".",
              " Probably missing one or more ", delim[1], ".")
      }
      # Have equal number of paired delimiters, so ready to start. Haven't
      # verified delimiter come in correct order, that will be done as we
      # process them in pairs.

      # Split this string into pieces at each delimiter (begin AND and). Some
      # string pieces may be 0 if the string begins and/or ends with a delimiter
      pieces <- character(2 * length(starts[[stringNum]]) + 1)

      # First piece is string up to first open deliminter
      pieces[1] <- substr(x[stringNum], 1, starts[[stringNum]][1]-1)

      # Remianng pieces come in pairs: between open and close delimiter (the
      # text to process as a template), and, except for the last close delimiter,
      # the part between the close delimiter and the next close delimiter.
      for (fieldNum in 1:length(starts[[stringNum]])) {
         # This pair of delimiters comes in the correct order, or die.
         if (starts[[stringNum]][fieldNum] > ends[[stringNum]][fieldNum]) {
            stop(delim[2], " before ", delim[1], " in string ", stringNum, ".")
         }
         # This is not the last pair of delimiters, so check for nexted delimiters
         # (the *next* start can't come before this open delimiter's paired close)
         if (length(starts[[stringNum]]) > fieldNum) {
            if (starts[[stringNum]][fieldNum + 1] < ends[[stringNum]][fieldNum]) {
               stop("Nested delimiters not allowed: ", delim[1], " occurs again before ", delim[2], " in string ", stringNum, ".")
            }
         }
         # Yay, we finally have a guaranteed good delimiter pair. Get the contents
         # of the string between these delimiters (the template text) as "field"
         fieldStart <- starts[[stringNum]][fieldNum] + attr(starts[[stringNum]], 'match.length')[fieldNum]
         fieldEnd <- ends[[stringNum]][fieldNum] - 1
         field <- substr(x[stringNum], fieldStart, fieldEnd)

         if (as.R) {
            # Evalute template text as R code
            pieces[2*fieldNum] <- eval(parse(text=field), envir= envir, enclos=envir)
         }
         else {
            # Evaluate template text as a variable name
            pieces[2*fieldNum] <- get(field, envir= envir, inherits=TRUE)
         }
         nonfieldStart <- ends[[stringNum]][fieldNum] + attr(ends[[stringNum]], 'match.length')[fieldNum]
         if (length(starts[[stringNum]]) > fieldNum) {
            nonfieldEnd <- starts[[stringNum]][fieldNum+1] - 1
         }
         else {
            nonfieldEnd <- nchar(x[stringNum])
         }
         pieces[2*fieldNum + 1] <- substr(x[stringNum], nonfieldStart, nonfieldEnd)
      }
      retVal[stringNum] <- paste0( pieces, collapse="")
   }
   return(retVal)
}