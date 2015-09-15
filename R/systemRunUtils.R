
#' Verify that a fileName is safe for use in shell commands.
#'
#' Validates filenames for use in system commands. Allows for either blacklist
#' and whitelist filtering, and independently controls several characters in
#' common use in filenames that could be dangerous in some contexts. Optionally
#' may also require a filename to exist or not to exist.
#
#' @section Valid names:
#'
#'    Filtering at a single character level will catch most bad names:
#'
#'    \tabular{ll}{
#'       Always bad \tab \command{` $ ( ) | ; & > <}\cr
#'       Always ok \tab \command{ A:Z a:z 0:9 _ . - / }\cr
#'       Can set bad \tab \command{" "} (space)\cr
#'       Can set bad \tab \command{~} (tilde will be allowed anywhere)\cr
#'       Can set ok \tab \command{\\ , :}\cr
#'       May not start with \tab \command{+ -}\cr
#'    }

#'    In addition, no filename may contain either \code{" -"} or \code{" +"} ( the
#'    \code{-} or \code{+} character preceded by a space.) This prevents injecting
#'    new options.
#'
#'    By default, whitelist validation is used. Files may only use charcters
#'    explicitly allowed. Can switch to blacklist validation by setting \code{weak=
#'    TRUE}, which allows any character not explicitly banned.
#'
#' @param names A vector of filenames to check.
#'
#' @param exists NULL Can set \code{TRUE} to ensure file exists, or \code{FALSE}
#' to ensure file does NOT exist. (Note, don't rely on this as only error check for
#' reading/writting file as that introduces a race condition).
#'
#' @param weak default this is \code{FALSE} and
#' validation is done using white-listed characters. Set \code{TRUE} to do
#' validation using blacklisted characters. See Valid names below for
#' additional info
#
#' @param okSpace Set \code{TRUE} to ban a space. Allows " " by default.
#'
#' @param okTilde Set \code{TRUE} to ban a tilde. Allows "~" by default
#'   (anywhere!).
#'
#' @param okBackslash Set \code{TRUE} to allow a backslash. No \code{\\} are
#'   allowed by default.
#'
#' @param okComma Set \code{TRUE} to allow a comma. No \code{,} are allowed by
#'   default.
#'
#' @param okColon  Set \code{TRUE} to allow a colon. No \code{:} are allowed by
#'   default.
#'
#' @return Returns a vector of booleans the same length as the input names,
#'   TRUE if the file is safe, FALSE if it fails any test.
#'
#' @export
is.safeFileName <- function( names,
                             exists= NULL, weak=FALSE,
                             okSpace=TRUE, okTilde=TRUE,
                             okBackslash=FALSE, okComma=FALSE, okColon=FALSE
) {

   badChar = c( "`", "$", "(", ")", "|", ";", "&", ">", "<" );

   # Add selectable bad char
   maybeBadChar = c(" ", "~", "\\", ",", ":");
   selector = c(! okSpace, ! okTilde, ! okBackslash, ! okComma, !okColon);
   badChar = c(badChar, maybeBadChar[selector]);

   # Create regexp character class matching bad characters
   badCharClass = paste( c("[", badChar, "]"), collapse= "");

   # Create regexp that identifes things that might look like options.
   noInjectOptions = " [-+]|^[-+]";

   # Create logical list, anything matching this regexp is bad
   matchBadRE =  paste( c(badCharClass, "|", noInjectOptions), collapse= "");
   nameSelector = logical();
   nameSelector <- ! grepl( matchBadRE, names );

   # If using white list, also set bad those that don't match this
   # regexp. All maybeBad characters are allowed here, if not to be
   # allowed they are already marked bad above.
   # NOTE: must start class with "-" this to work correctly!
   if (! weak) {
      okChar = c("A-Z", "a-z", "0-9");
      okChar = c( "-", okChar, "_", ".", "/", maybeBadChar );
      okCharClass = paste( c("[", okChar, "]"), collapse= "");
      matchGoodRE = paste( c("^", okCharClass, "+$"), collapse= "");
      nameSelector = nameSelector & grepl( matchGoodRE, names );
   }

   # If exists= is set (true/false), file allowed only if exists/does not exist.
   if (! is.null(exists)) {
      if (exists) {
         nameSelector = nameSelector & file.exists(names);
      }
      else {
         nameSelector = nameSelector & ( ! file.exists(names) );
      }
   }

   return( nameSelector );
}
