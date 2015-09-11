

# Given a file of sample names, returns those names as a vector
getSampleCohort <- function( file ) {
   ###
   #     PARAM: file= (string) A file of samples names, one per line. Blank and
   # lines starting with "#" are ignored.
   ###
   #     RETURNS: Vector of sample names
   ###

   samples <- scan( file, what= character(0), comment.char= "#" );
   if (any(duplicated(samples))) {
      stop("Sample file may not contain duplicates: ", file, "\nDuplicated: ",
           paste( samples[duplicated( samples )], collapse=", " ));
   }
   return( samples );
}

# Read file of sample, data file name and return a data frame containing those
# samples given in "samples".
getCohortFiles <- function ( file, samples) {
   ###
   #     PARAM: file= (string) A file with one row per sample. at least one
   # column and a header line. The first coloumn must contain sample names
   # (mathcing sample). The remaining columns should contain the full path file
   # name of the data files. Columns should be tab separated.
   #     PARAM: samples= (character vec) A vector of strings giving the sample
   # names to extract rows from the list for
   ###
   #     RETURNS: A data frame with one column for eacn column in the input file
   # and one row for each row in the data file matching a sample name in the
   # samples list.
   ###
   allDataFile <- read.delim( file, stringsAsFactors=FALSE );
   select <- allDataFile[,1] %in% samples
   cohortDataFiles <- allDataFile[select,];
   return( cohortDataFiles );
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
