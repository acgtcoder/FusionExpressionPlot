
#' Read in the cohort definition file
#'
#' Creates a data frame from a tab-delimited file with a header describing a
#' sample cohort. The cohort file must specify at minimum the sample and
#' associated exon expression file. Allows selecting which columns are which and
#' also only selecting a subset of the samples to keep. Will skip comment lines.
#'
#' @param file The name of the cohort file to read in (tab-delimited with
#'   header)
#'
#' @param samples The list of samples to read in, as a string vector. The rows
#'   matching these are the only ones read in. By default, this is set to
#'   \code{NULL}, meaning it reads all the samples. If any sample in the
#'   provided list is not found, a warning will be generated. Must match exactly
#'   (case sensitive).
#'
#' @param columns The names of the columns to read the sample name and exon
#'   expression filename from. By default this is \code{c("sample",
#'   "exonExpressionFile")} If provided, must specify both, and must match
#'   exactly (case sensitive).
#'
#' @param dataRoot A directory to prepended to the \code{exonExpressionFile}
#'   filenames. By default this is \code{NULL} and nothing is prepended. It is
#'   ok if \code{dataRoot} does not exist at this time as the file locations are
#'   not verified at this point, but \code{NA} and empty strings \code{""} are
#'   ignored.
#'
#' @param comment.char The character starting comment lines, by default
#'   \code{'#'}. Must be a single character. Any line begining with this
#'   character is considered to be a comment line and is ignored. Must be the
#'   first character in a line (no leading white space is allowed). Set to the
#'   empty string, \code{""}, to skip comment line filtering.
#'
#' @return A data frame with two columns:
#'
#' \tabular{ll}{
#'    \code{sample}
#'       \tab The sample names read from the file, by default from the
#'             \code{sample} column\cr
#'   \code{exonExpressionFile}
#'       \tab The exon expression file names as read from the file, by default
#'             from the \code{exonExpressionFile} column\cr
#' }
#'
#' @section Errors:
#'
#' \describe{
#'    \item{
#'       \command{Can't find the cohort definition file: \var{file}}
#'    }{
#'       The file you are trying to read in can't be found - case, permission,
#'       relative directory, and misspellings are all possible reasons.
#'    }
#'    \item{
#'       \command{Columns must use strings for names.}
#'    }{
#'       You specified something other than a string when setting the
#'       \code{columns} parameter. Perhaps you need something like \code{"1"}
#'       instead of \code{1}.
#'    }
#'    \item{
#'       \command{Must specify both the sample and exonExpressionColumn headers}
#'    }{
#'       You can't specify only one of these and expect the other to be the
#'       default. You have to specify both.
#'    }
#'    \item{
#'       \command{Columns must specify non-empty strings for column names}
#'    }{
#'       Can't specify an empty string as a column heading. What are you trying
#'       to do? This must be a real column heading from your file that can be
#'       used to identify the column after being read in.
#'    }
#'    \item{
#'       \command{May not specify the same column for sample and exonExpressionFile data}
#'    }{
#'       Why would you use the same heading for both needed columns? Don't do that.
#'    }
#'    \item{
#'       \command{The cohort file has no sample column \var{sample}}
#'    }{
#'       You have to identify which column in the cohort file can be used as
#'       sample names
#'    }
#'    \item{
#'       \command{The cohort file has no exon expression file column
#'          \var{exonExpressionColumn}}
#'    }{
#'       You have to identify which column in the cohort file has the file names
#'       for the exon expression data.
#'    }
#'    \item{
#'       \command{The cohort file can not contains duplicate samples}
#'    }{
#'       The cohort file contains multiple lines with the same sample name, even
#'       after filtering out any samples you didn't want considered. That's not
#'       ok.
#'    }
#'    \item{
#'       \command{The cohort file can not contain duplicate cohort expression files}
#'    }{
#'       The cohort file contains multiple lines with the same cohort expression
#'       file name, even after filtering out any samples you didn't want
#'       considered. That's not ok. If two different samples really have the
#'       same expression file name, you'll have to put them in different
#'       directories.
#'    }
#' }
#'
#' @section Warnings:
#'
#' \describe{
#'    \item{
#'       \command{Duplicates in "samples" parameter filter list ignored}
#'    }{
#'       The list of samples you provided to filter the cohort data frame by
#'       contained the same sample name more than once. You can only include
#'       a sample once, so the duplicates are just ignored. Warning you allows
#'       this to be changed for future use (if you care).
#'    }
#'    \item{
#'       \command{Ignoring missing samples specified in the "samples" parameter
#'                list: \var{sample, sample,...}}
#'    }{
#'       The list of samples you provided to filter the cohort data frame
#'       contained sample names that were not actually in the cohort file.
#'       Perhaps you should double check your cohort file?
#'    }
#'}
#'
#' @export
loadCohortDefinition <- function ( file, samples=NULL, comment.char= '#',
   columns= c("sample", "exonExpressionFile"), dataRoot= NULL ) {

   if (! file.exists(file)) {
      stop("Can't find the cohort definition file: ", file)
   }

   if (! is.null(samples)) {
      if (anyDuplicated(samples)) {
         samples <- unique(samples)
         warning( "Duplicates in \"samples\" parameter filter list ignored" )
      }
   }

   if (! is.null(dataRoot) ) {
      if (is.na(dataRoot) || dataRoot == "") {
         dataRoot <- NULL
      }
   }

   if (! class(columns) == 'string') {
      stop ("Columns must use strings for names.")
   } else if (length(columns) != 2) {
      stop ("Must specify both the sample and exonExpressionColumn headers")
   } else if (nchar(columns[1]) < 1 || nchar(columns[2]) < 1) {
      stop( "Columns must specify non-empty strings for column names")
   } else if (columns[1] != columns[2]) {
      stop("May not specify the same column for sample and exonExpressionFile data")
   }

   df <- read.delim(
      file, comment.char= comment.char, stringsAsFactors= FALSE
   )

   cohortHeadings <- names(df)
   if (! columns[1] %in% cohortHeadings) {
      stop("The cohort file has no sample column", columns[1])
   }
   if (! columns[2] %in% cohortHeadings) {
      stop("The cohort file has no exonExpressionFile column", columns[2])
   }

   if (! is.null(samples)) {
      samplesFound <- samples %in% df[,columns[1]]
      if (! all(samplesFound)) {
         missing <- samples[ ! samplesFound]
         missing <- paste(missing, collapase= ", ")
         warning( "Ignoring missing samples specified in the \"samples\" parameter list:\n",
                  missing )
      }
      df <- df[columns[1] %in% samples, ]
   }

   if (anyDuplicated(df[,columns[1]])) {
      stop( "The cohort file can not contains duplicate samples")
   }
   if (anyDuplicated(df[,columns[2]])) {
      stop( "The cohort file can not contain duplicate cohort expression files")
   }
   return( df );
}

# Loads an exon expression data file column as a data frame
loadExonExpressionFile <- function (path, type= "rpkm" ) {
   ###
   # Loads a single exon expression data file.
   ###
   #     PARAM path: path to the data file to load
   #     PARAM type: column from data file to load
   ###
   if (! file.exists( path )) {
      stop( "No such file: ", path );
   }

   okType = c( "rpkm", "count", "coverage" );
   if (! type %in% okType) {
      stop( "Only type values allowed are: ", paste(okType, collapse = ", " ));
   }

   columnNamesRead <- c( "exon", "count", "coverage", "rpkm" );
   columnTypesRead <- c("character", "integer", "numeric", "numeric");
   df <- read.delim( path, header=FALSE, col.names=columnNamesRead, colClasses = columnTypesRead, stringsAsFactors=FALSE );
   exonExpression = data.frame(
      chr      = sub( ":.+",     "", df$exon ),
      start    = sub( "^.+:",    "", sub( "-.+[-+]$", "", df$exon )),
      end      = sub( ":.+$",    "", sub( "^.+[^:]-", "", df$exon )),
      strand   = sub( "^.+:.+:", "", df$exon ),
      count    = df$count,
      coverage = df$coverage,
      rpkm     = df$rpkm,
      stringsAsFactors= FALSE);
   exonExpression <- exonExpression[ , c( "chr", "start", "end", "strand", type )];
   return( exonExpression );
}

# Loads an exon expression data file and appends it to an existing exon
# expression data frame
addExonExpression <- function( exonModels, id, path, type="rpkm" ) {
   ###
   # Adds the expression data from a single exon data file to a list of gene
   # models. Can be used to grow a file of all exon data on sample at a time.
   ###
   #     PARAM: exomModels - the exon models data frame to append to
   #     PARAM: id - the name of the column when appended, must be using
   # in the data frame
   #     PARAM: path - the path to the data file to add a column from
   #     PARAM: type - the name of the data column in the file to add
   ###

   if (id %in% colnames( exonModels )) {
      stop( "column ", id, " already in the models data set" );
   }
   # path and type validated in subroutine loadExonExpressionFile()
   exonExpression <- loadExonExpressionFile( path, type );

   colnames( exonExpression )[ which( colnames( exonExpression ) == type )] = id;

   exonExpression = merge(
      exonModels, exonExpression, by=c( "chr", "start", "end", "strand" ), all.x=TRUE
   );
   if (any(is.na(exonExpression[,id]))) {
      warning( "Expression file for ", id, "is missing some exons" );
   }

   return( exonExpression );
}

# Loads data
getCohortExonExpressionData <- function ( geneModels, cohortFiles, progress=TRUE, type="rpkm" ) {
   # Get all the expression files for a cohort of samples and build a single
   # table. This takes a while, so progress by sample is reported by default.
   ###
   #   PARAM: geneModels - the gene models data fram from load gaf models
   #   PARAM: cohortFiles - the data frame giving samples and data file names
   #   PARRAM: progress - show a progress bar
   #   PARAM: type - the name of the column from the files in the cohort to load
   ###
   #   Returns a data frame with the gene models + one column per sample in the
   # cohort Files
   ###

   exonExpression <- geneModels;
   for ( rowNum in 1:nrow(cohortFiles) ) {
      if (progress) {
         progressMessage = paste(
            c( rowNum, " of ", nrow(cohortFiles),
               ". Adding exon expression data for ",  cohortFiles[rowNum, "samples"]
            ), collapse = ""
         );
         print( progressMessage );
      }
      exonExpression <- addExonExpression(
         exonModels= exonExpression,
         id=   cohortFiles[rowNum, "samples"],
         path= cohortFiles[rowNum, "exon_expression_files"],
         type= type
      );
   }
   sortSelect <- order(
      exonExpression$gene,
      exonExpression$chr,
      exonExpression$strand,
      exonExpression$start,
      exonExpression$end
   );
   exonExpression <- exonExpression[sortSelect,];

   return( exonExpression )

   # TODO: Need to fix headers at some point back from "." to "-".
}


# Normalizes data
normExpressionData <- function( exonExpressionData ) {

   # Data columns are 10+
   sampleColumnNums = 10:ncol( exonExpressionData );
   infoColumnNums = 1:9;

   # scale data columns per gene
   geneScaledExpression <- by( exonExpressionData,
                               exonExpressionData$gene,
                               function(df) { scale( df[ , sampleColumnNums ], center=T, scale=T ); }
   );

   # convert back into data frame
   geneScaled <- data.frame( do.call( "rbind", geneScaledExpression ));

   # Scaling introduced NaN if 0 variance, want these to be all 0.
   # Non-uniformity strikes again is.nan different from is.na, no
   # data frame method!

   is.nan.data.frame <- function(x) {
      do.call(cbind, lapply(x, is.nan))
   }

   geneScaled[is.nan(geneScaled)] <- 0;

   # Bind rest of data back to get complete df for possible use later.
   geneScaled <- data.frame( exonExpressionData[ , infoColumnNums ], geneScaled );

   # Scale the exons
   # Since scaling all exons (not applying by factor), this is easier,
   # although scale only works by columns, so need to transpose before and after.
   cohortScaled <- data.frame( geneScaled[ , infoColumnNums],
                               t( scale( t( geneScaled[ , sampleColumnNums] ), center=T, scale=T)));

   # This is a hack; Need to figure out why converts names and be REALLY sure it doesn't change order.
   names(cohortScaled) <- names(exonExpressionData)

   return(cohortScaled)
}
