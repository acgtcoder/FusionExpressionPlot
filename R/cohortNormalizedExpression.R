
#' Read in the cohort definition file
#'
#' Creates a data frame from the cohort description file. The file read must be
#' tab-delimited, one line per sample, with header. Lines beginning with a '#'
#' are ignored. It must include the columns 'sample' and 'exonExpressionFile'
#' (case sensitive). Allows loading only a subset of the samples. All columns
#' from the file are loaded.
#'
#' @param file The name of the cohort file to read in (tab-delimited with
#'   header)
#'
#' @param samples The list of samples to read in, as a character vector. The rows
#'   matching these are the only ones read in. By default, this is set to
#'   \code{NULL}, meaning it reads all the samples. If any sample in the
#'   provided list is not found, a warning will be generated. Must match exactly
#'   (case sensitive).
#'
#' @param dataRoot A directory to prepended to the \code{exonExpressionFile}
#'   filenames. By default this is \code{NULL} and nothing is prepended.
#'   \code{NA} and empty text \code{""} are liekwise ignored. It is ok if
#'   \code{dataRoot} does not exist at this time as the file locations are not
#'   verified at this point.
#'
#' @param comment.char The character starting comment lines, by default
#'   \code{'#'}. Must be a single character. Any line beginning with this
#'   character is considered to be a comment line and is ignored. Must be the
#'   first character in a line (no leading white space is allowed). Set to the
#'   empty text, \code{""}, to skip comment line filtering.
#'
#' @return A data frame representing the sample cohort with at least two
#'   columns:
#'
#' \tabular{ll}{
#'    \code{sample}
#'       \tab The sample names read from the input file's \code{sample} column\cr
#'   \code{exonExpressionFile}
#'       \tab The exon expression file names read from the input file's
#'             \code{exonExpressionFile} column\cr
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
#'       \command{The cohort file has no sample column \var{sample}}
#'    }{
#'       The cohort file must have a column \code{sample} with the sample names
#'    }
#'    \item{
#'       \command{ The cohort file has no expression file column \var{exonExpressionFile}}
#'    }{
#'       The cohort file must have a column \code{exonExpressionFile} giving the file names
#'       of the exon expression data for the sample in the same row.
#'    }
#'    \item{
#'       \command{The sample and exon expression file columns must contain text}
#'    }{
#'       Either or both of these columns was read in as something other than
#'       a character vector.
#'    }
#'    \item{
#'       \command{The cohort file can not contains duplicate samples}
#'    }{
#'       The cohort file contains multiple lines with the same sample name, even
#'       after filtering out any samples you didn't want considered. That's not
#'       ok.
#'    }
#'    \item{
#'       \command{The cohort file can not contain duplicate cohort expression
#'          files}
#'    }{
#'       The cohort file contains multiple lines with the same cohort expression
#'       file name, even after filtering out any samples you didn't want
#'       considered. That's not ok. If two different samples really have the
#'       same expression file name, you'll have to put them in different
#'       directories.
#'    }
#'    \item{
#'       \command{All sample and exon expression file entries must be non-empty
#'          text}
#'    }{
#'       Can not have missing or empty text for samples as this is the
#'       primary key for later work. Can not have missing or empty text for
#'       exon expression files (after filtering) as such samples can not then
#'       be part of the cohort.
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
#'          list: \var{sample, sample,...}}
#'    }{
#'       The list of samples you provided to filter the cohort data frame
#'       contained sample names that were not actually in the cohort file.
#'       Perhaps you should double check your cohort file?
#'    }
#'}
#'
#' @section Todo:
#' \itemize{
#'    \item{Consider general validator tool with key column and column
#'    validation selectable or specifiable. }
#' }
#' @export
loadCohortDefinition <- function ( file, samples=NULL, comment.char= '#',
                                   dataRoot= NULL ) {

   if (! file.exists(file)) {
      stop("Can't find the cohort definition file: \"", file, "\"")
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

   df <- read.delim(
      file, comment.char= comment.char, stringsAsFactors= FALSE
   )

   cohortHeadings <- names(df)
   if (! "sample" %in% cohortHeadings) {
      stop("The cohort file has no sample column.")
   }
   if (! "exonExpressionFile" %in% cohortHeadings) {
      stop("The cohort file has no exonExpressionFile column.")
   }
   if (class(df$sample) != 'character' || class(df$exonExpressionFile) != 'character') {
      stop("The sample and exon expression file columns must contain text")
   }

   if (! is.null(samples)) {
      samplesFound <- samples %in% df$sample
      if (! all(samplesFound)) {
         missing <- samples[ ! samplesFound]
         missing <- paste(missing, collapase= ", ")
         warning( "Ignoring missing samples specified in the \"samples\" parameter list:\n",
                  missing )
      }
      df <- df[df$sample %in% samples, ]
   }

   if (anyDuplicated(df$sample)) {
      stop( "The cohort file can not contains duplicate samples")
   }
   if (anyDuplicated(df$exonExpressionFile)) {
      stop( "The cohort file can not contain duplicate cohort expression files")
   }
   if (any(is.na(df$sample)) || any(nchar(df$sample) < 1)
       || any(is.na(df$exonExpressionFile)) || any(nchar(df$exonExpressionFile) < 1)) {
      stop("All sample and exon expression file entries must be non-empty text")
   }

   if (! is.null(dataRoot)) {
      df$exonExpressionFile <- file.path(dataRoot, df$exonExpressionFile)
   }
   return( df );

}

#' Load a TCGA exon expression data file as a data frame
#'
#' Loads an RNA exon expression data file as generated for the TCGA, returning
#' its contents as a data frame. The input file has four tab-delimited columns,
#' no header. The first column in the file describes the exon as
#' <chr>:<start>-<end><strand>, the next three give count, coverage, and rpkm
#' values.
#'
#' @param path The path to the exon expression data file to load
#'
#' @return A data frame with the following columns from the exon expression
#'   file:
#'
#' \tabular{ll}{
#'    \code{chr}      \tab The exon's chromosome\cr
#'    \code{start}    \tab The genomic start coordinate for the exon\cr
#'    \code{end}      \tab The genomic end coordinate for the exon\cr
#'    \code{strand}   \tab The strand the chromosome is on, one of \code{+ | - | *}\cr
#'    \code{count}   \tab One of the exon expression level columns\cr
#'    \code{coverage} \tab One of the exon expression level columns\cr
#'    \code{rpkm}     \tab One of the exon expression level columns\cr
#' }
#'
#' @section Errors:
#'
#' The following fatal errors can occur:
#'
#' \describe{
#'    \item{
#'       \command{No such file: \var{path}}
#'    }{
#'       The exon expression file named can not be found. Case, spelling
#'       permission and wrong working dir for relative paths are all common
#'       errors.
#'    }
#' }
#'
#' @export
loadExonExpressionFile <- function (path) {
   if (! file.exists( path )) {
      stop( "No such exonExpressionFile: ", path );
   }

   columnNamesRead <- c( "exon", "count", "coverage", "rpkm" );
   columnTypesRead <- c("character", "integer", "numeric", "numeric");
   df <- read.delim( path, header=FALSE, col.names=columnNamesRead, colClasses = columnTypesRead, stringsAsFactors=FALSE );
   exonExpression = data.frame(
      chr      = sub( ":.+",     "", df$exon ),
      start    = as.integer(sub( "^.+:",    "", sub( "-.+[-+]$", "", df$exon ))),
      end      = as.integer(sub( ":.+$",    "", sub( "^.+[^:]-", "", df$exon ))),
      strand   = sub( "^.+:.+:", "", df$exon ),
      count    = df$count,
      coverage = df$coverage,
      rpkm     = df$rpkm,
      stringsAsFactors= FALSE);
   return( exonExpression );
}

#' Add exon expression file data to an existing data frame
#'
#' Adds a column from a single exon data file to gene models data frame or an
#' exon expression data frame. Can be used to grow a cohort-wide file of
#' expression data one sample at a time by recursively using the returned data
#' frame as the input \code{exonModels} file.
#'
#' @param exonModels The existing exon models data frame to append to, with or
#'   without additional columns. The data.frame returned by this function is a
#'   valid input for this parameter. See also \code{\link{loadCohortDefinition}}
#'
#' @param id The name of the column in the data frame after it is appended
#'
#' @param path The path to the exon expression data file being added.
#'
#' @param type The name of the column in the exon expression data file that is
#'   extracted and appended. By default this is \code{"rpkm"}
#'
#' @return A data frame with one row per exon with the initial columns
#' describing the exons (as loaded by \code{\link{loadCohortDefinition}}), any
#' following columns of data, and then the specified column from the given
#' exon expression data file appended.
#'
#' @export
addExonExpression <- function( exonModels, id, path, type="rpkm" ) {

   str(exonModels)

   message("Incoming Header: ", paste(header, collapse="," ))
   if (id %in% header) {
      stop( "column ", id, " already in the models data set" );
   }

   # Get a data frame with exon descriptions and only the columns wanted.
   keepColumns <- c( "chr", "start", "end", "strand", type)
   # path validated in subroutine loadExonExpressionFile()
   exonExpression <- loadExonExpressionFile( path )[keepColumns]
   colnames(exonExpression[type]) <- id
   message("exonExpression Loaded: ", exonExpression[1:3,])
   str(exonExpression)

   exonExpression <- merge(
      exonModels, exonExpression,
      by=c( "chr", "start", "end", "strand" ), all.x=TRUE
   )
   message("exonExpression after merge: ", exonExpression[1:3,])
   str(exonExpression)

   if (any(is.na(exonExpression[,id]))) {
      warning( "Expression file for ", id, "is missing some exons" );
   }
   return( exonExpression );
}

#' Generate a data frame of cohort exon expression data
#'
#' Load all the exon expression files for a cohort of samples and build a single
#' data frame with exon information from a geneModels data frame and one column
#' of expression information per sample in the cohort. This takes a while, so
#' progress by sample is reported by default.
#'
#' @param geneModels A data frame describing the exon structure of the genes
#'   (gene models). See \code{\link{loadGafModels}}. This must correspond to the
#'   exon expression file exon names.
#'
#' @param cohortFiles A data frame describing the samples and their exon
#'   expression files. See \code{\link{loadCohortDefinition}}. All files must
#'   exist, a simple check is \code{all( file.exists(
#'   \var{cohortFiles}$exonExpressionFile))}.
#'
#' @param progress By default shows a progress bar. Set this to \code{FALSE} if
#'   don't want this (e.g. if you are logging output when running this in a
#'   batch job).
#'
#' @param type The type of expression being used. This is the name of the column
#'   from the exon expression files to use. By default uses the column named
#'   "rpkm". See \code{\link{loadExonExpressionFile}}
#'
#' @return Returns a data frame that concatenates the columns from the gene
#' models data frame with one column from each sample in the cohort data frame.
#'
#' @export
getCohortExonExpressionData <- function ( geneModels, cohortFiles, progress=TRUE, type="rpkm" ) {

   exonExpression <- geneModels;
   for ( rowNum in 1:nrow(cohortFiles) ) {
      if (progress) {
         progressMessage = paste(
            c( rowNum, " of ", nrow(cohortFiles),
               ". Adding exon expression data for ",  cohortFiles[rowNum, "sample"]
            ), collapse = ""
         );
         print( progressMessage );
         str(exonExpression)
      }
      exonExpression <- addExonExpression(
         exonModels= exonExpression,
         id=   cohortFiles[rowNum, "sample"],
         path= cohortFiles[rowNum, "exonExpressionFile"],
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


#' Normalize a cohort of exon expression data
#'
#' @param exonExpressionData A data frame of exon expression data, with one
#'   column per sample in the cohort. See
#'   \code{\link{getCohortExonExpressionData}}.
#'
#' @return A data frame of exon expression data, but with the data normalized
#'   (scaled to a standard normal) first within one sample by gene, then across
#'   all samples by exon.
#'
#' @export
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
