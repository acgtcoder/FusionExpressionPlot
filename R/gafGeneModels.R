
#' Extract gene models from the GAF
#'
#' This extracts a tab delimited file of all canonical exon gene models from the
#' TCGA gene annotation file or \acronym{GAF}. The canonical model is intended to be at
#' least a common starting point for a simplified definition of "the gene called
#' X", based on the union of all transcript. By default each model in the output
#' file (row) is uniquely identified by gene name.
#'
#' This function is implemented using a unix system command and requires the
#' "grep" program, so this only works on linux/mac systems. [TODO - reimplement
#' in pure R.] The GAF version this works with is the version used for the
#' TCGA RNAseq expression data files. It is available for download at the NCI
#' uncompressed:
#' \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf}{TCGA.hg19.June2011.gaf}
#' or gzipped:
#' \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf.gz}{TCGA.hg19.June2011.gaf.gz}
#'
#' @param gaf The full-path name of the GAF [REQ]
#'
#' @param outFile The output filename. By default will be the same as the
#'   GAF, with ".geneModels" appended (hence it will be created in the same
#'   directory by default). This will not overwrite an existing file unless
#'   \code{force = TRUE}.
#'
#' @param force If the output file exists, setting this \code{TRUE} will allow
#' overwriting it. Doing so generates a warning.
#'
#' @param uniqueGene By default this is \code{TRUE} and only one copy of every
#'   gene will be kept. This makes the gene name a unique key. The GAF contains
#'   additional versions of some genes which are skipped. See the GAF gene names
#'   section below for more information.
#'
#' @param skipUnknownGene By default this is \code{TRUE} and the unknown genes
#'   (those whose names in the GAF are "?") are dropped. See the GAF gene names
#'   section below for more information. Note that this setting is over-ridden
#'   to \code{TRUE} with a warning if \code{uniqueGene= TRUE} as all these genes
#'   have the same gene name ("?").
#'
#' @return The main output is the GAF gene mode extract file, which is by
#'   default just the genes with unique names. See
#'   \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF_bundle_Feb2011/docs/GAF_v2_file_description.docx}{GAF_v2_file_description.docx}
#'    However, this function will returns a list with some summary info about
#'   the GAF and the generated geneModels file:
#'
#'   \tabular{ll}{
#'      \code{$gaf}        \tab The GAF filename used as a parameter (relative filename).\cr
#'      \code{$gaf_real}   \tab The absolute full path filename to the input GAF.\cr
#'      \code{$gaf_md5}    \tab The md5 checksum of the GAF, as a string.\cr
#'      \code{$uniqueGene} \tab The \code{uniqueGene} parameter setting used.\cr
#'      \code{$skipUnknownGene} \tab The \code{skipUnknownGene} parameter setting used.\cr
#'      \code{$gaf_lines}       \tab The number of lines in the GAF.\cr
#'      \code{$gaf_models}      \tab The number of models in the GAF. [Currently this
#'         does not include unknown genes if \code{skipUnknownGene= TRUE}].\cr
#'      \code{$gaf_models_unique} \tab The number of unique models in the GAF. Will
#'         be \code{NA} if \code{uniqueGene= FALSE}.\cr
#'      \code{$gaf_extract}      \tab The output filename, based on the input GAF by default.\cr
#'      \code{$gaf_extract_real} \tab The absolute full path filename of the output file.\cr
#'      \code{$gaf_extract_md5}  \tab The md5 checksum of the output file.\cr
#'   }
#'
#' @section GAF gene names:
#'   Some genes in the GAF have multiple variants with the same name. These are
#'   annotated like "1ofN", "2ofN", etc. By default only the "...1of" variant
#'   is kept. Setting \code{uniqueGene= FALSE} will keep all variants. Gene name will
#'   not then be a unique key.
#'
#'   Some genes in the GAF have no known name and are annotated as "?". Some of
#'   these have multiple versions also. There are 32 such genes in the GAF. By
#'   default non of these are kept. If \code{uniqueGene= FALSE} is set, these
#'   will still be skipped unless \code{skipUnknownGene= FALSE} is also set.
#'   There is no way to keep these genes while skipping the variants of the
#'   named genes.
#'
#'   The gene "SLC35E2" is present twice, but without a numbered annotation. By
#'   default, only the larger (encompassing) version of "SLC35E2" is kept.
#'   If \code{uniqueGene= FALSE} is set, then both versions of this gene will be
#'   kept, regardless of the \code{skipUnknownGene} setting.
#
#' @section Errors:
#'    These errors are fatal and will terminate processing.
#'
#' \describe{
#'    \item{
#'       \command{Unsafe character in GAF filename!}
#'    }{
#'       An invalid characters was passed as part of the GAF filename. This is
#'       important as the filename is used in a system command as a parameter
#'       and could be used for command injection.
#'    }
#'    \item{
#'       \command{Can't find the specified GAF: "\var{file}"}
#'    }{
#'       The specified GAF doesn't seem to exist on the file system.
#'       Probably have the name wrong or are using a relative name from the
#'       wrong directory, but could also be that permissions are hiding it.
#'    }
#'    \item{
#'       \command{Unsafe character in output geneModel filename!}
#'    }{
#'       An invalid characters was passed as part of the output filename. This
#'       is important as the filename is used in a system command as a parameter
#'       and could be used for command injection.
#'    }
#'    \item{
#'       \command{Output file already exists; use force= TRUE to overwrite: "\var{file}"}
#'    }{
#'       The specified GAF gene extract output file already exists. You probably
#'       don't want to overwrite it. However, you can set \code{force= TRUE}
#'       to allow this. It will still generate a warning.
#'    }
#' }
#'
#' @section Warnings:
#'
#' \describe{
#'    \item{
#'       \command{Forcing overwrite of output file: "\var{file}"}
#'    }{
#'       Just letting you know an existing file is actually being overwritten.
#'       This won't happen unless explicitly allowed by setting \code{force=
#'       TRUE}). Having a warning allows distinguishing between the cases where
#'       an overwrite occurred vs those where one was allowed but did not occur.
#'    }
#'    \item{
#'       \command{uniqueGene=TRUE sets skipUnknownGene=TRUE}
#'    }{
#'       You can't have unique gene names if you keep the genes without names. I
#'       could just blow up, but I'm just going to assume that since you asked
#'       for unique gene names, that's what you really want. That means I have
#'       to ignore your request to keep the unknown genes. Not what you wanted?
#'       That's why I'm warning you.
#'    }
#'    \item{
#'       Various warnings from failed system commands
#'    }{
#'       System commands are used for several things in this function. If they
#'       fail, error messages are returned as warnings.
#'    }
#'
#' }
#'
#' @section Todo:
#'    \itemize{
#'       \item{Add test for corner case - one exon gene}
#'    }
#'
#' @examples
#' \dontrun{
#'
#' # Extract gene models to the default output file
#' stats <- extractGeneModels( 'path/to/GAF' )
#'
#' # Same, with all defaults made explicit
#' stats <- extractGeneModels(
#'    gaf= 'path/to/GAF', outFile= 'path/to/GAF.geneModels', force= FALSE,
#'    uniqueGene= TRUE, skipUnknownGene= TRUE
#' )
#'
#' # Extract gene models to gaf.genes in run directory
#' stats <- extractGeneModels( 'path/to/GAF', outFile= 'gaf.genes' )
#'
#' # Overwrite outFile if it exists (here using the default name)
#' stats <- extractGeneModels( 'path/to/GAF', force= TRUE )
#'
#' # Extract all gene models, including duplicates and unknowns.
#' stats <- extractGeneModels(
#'    gaf= 'path/to/GAF', uniqueGene= FALSE, skipUnknownGene= FALSE
#' )
#'
#' # Extract all gene models except unknown (includes duplicates)
#' stats <- extractGeneModels( gaf= 'path/to/GAF', uniqueGene= FALSE)
#' )
#' }
#'
#' @export
extractGeneModels <- function( gaf,
                               outFile= paste0(gaf, ".geneModels"), force= FALSE,
                               uniqueGene= TRUE, skipUnknownGene= TRUE
) {


   # Using filename as a parameter in a system command so need to be careful.
   if (! is.safeFileName( gaf )) {
      stop( "Unsafe character in GAF filename!" );
   }
   if (! file.exists(gaf)) {
      stop( "Can\'t find the specified GAF: \"", gaf, "\"")
   }
   if (! is.safeFileName( outFile )) {
      stop( "Unsafe character in output geneModel filename!" );
   }
   if (file.exists(outFile)) {
      if (! force) {
         stop( "Output file already exists; use force= TRUE to overwrite: \"", outFile, "\"")
      } else {
         warning( "Forcing overwrite of output file: \"", outFile, "\"");
         unlink(outFile);
      }
   }

   if (uniqueGene && ! skipUnknownGene) {
      warning(" uniqueGene=TRUE sets skipUnknownGene=TRUE")
      skipUnknownGene=TRUE;
   }

   info = list();
   info$gaf = gaf;
   info$gaf_real = normalizePath(gaf);
   info$gaf_md5 = tools::md5sum(gaf);
   info$uniqueGene = uniqueGene;
   info$skipUnknownGene = skipUnknownGene;

   # Get number of lines in original GAF
   cmd = paste( c( "wc -l ", gaf ), collapse="" );
   got = system( cmd, intern=TRUE );
   info$gaf_lines = as.integer( sub( " .+", "", sub( "^[ ]+", "", got )));

   geneModelsFile = outFile;
   if (uniqueGene) {
      geneModelsFile = tempfile( "TEMP_gafExtractGenome" );
   }

   # Load header line into output
   cmd = paste(c("head -n1 ", gaf, " > ", geneModelsFile ), collapse="");
   system( cmd, intern=TRUE );

   dropUnknownGene = " | grep -v -E \"\t[?][|][0-9]\"";
   if (! skipUnknownGene) {
      dropUnknownGene = ""
   }

   # Extract gene models.
   cmd = paste( c( "grep genome ", gaf,
                   " | grep calculated",
                   " | grep chr",
                   dropUnknownGene,
                   " >> ", geneModelsFile ),
                collapse= "" );
   system( cmd, intern=TRUE );
   cmd = paste(c("wc -l ", geneModelsFile), collapse="");
   got = system( cmd, intern=TRUE );
   # Subtract one as has header in file
   info$gaf_models = as.integer( sub( " .+", "", sub( "^[ ]+", "", got ))) - 1;

   # If want truely unique models, need to ignore the 2of3, etc lines and the
   # duplicate SLC35E2 gene. Will keep the [1 of x] and the larger (encompasing)
   # version of SLC35E2.
   info$gaf_models_unique = NA
   if (uniqueGene) {
      # As grep for what doesn not match, will keep header line automatically
      cmd = paste(
         c( "grep -v -E \"[02-9]of|1[0-9]of|SLC35E2.9906\" ", geneModelsFile, " > ", outFile ),
         collapse="" );
      system( cmd, intern=TRUE );

      # Done with tempCalculate intermediate; calculate line count in out file
      unlink(geneModelsFile)
      cmd = paste(c("wc -l ", outFile), collapse="");
      got = system( cmd, intern=TRUE );
      # Subtract one as has header in file
      info$gaf_models_unique = as.integer( sub( " .+", "", sub( "^[ ]+", "", got ))) - 1;
   }

   info$gaf_extract= outFile;
   info$gaf_extract_real= normalizePath(outFile);
   info$gaf_extract_md5 = tools::md5sum(outFile);
   return( info )
}

#' Extract transcript models from the GAF
#'
#' This extracts a tab delimited file of all canonical transcript models from the
#' GAF. Each model in the output file (row) can be refereed to
#' uniquely by transcript name. This is much simpler than extracting genes
#'
#' This function is implemented using a unix system command and requires the
#' "grep" program, so this only works on linux/mac systems. [TODO - reimplement
#' in pure R.] The GAF version this works with is the version used for the
#' TCGA RNAseq expression data files. It is available for download at the NCI
#' uncompressed:
#' \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf}{TCGA.hg19.June2011.gaf}
#' or gzipped:
#' \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf.gz}{TCGA.hg19.June2011.gaf.gz}
#'
#' @param gaf The full-path name of the GAF [REQ]
#'
#' @param outFile The output filename. By default will be the same as the
#'   GAF, with ".transcriptModels" appended (hence it will be created in the same
#'   directory by default). This will not overwrite an existing file unless
#'   \code{force = TRUE}.
#'
#' @param force If the output file exists, setting this \code{TRUE} will allow
#' overwriting it. Doing so generates a warning.
#'
#' @return The number of transcripts in the output file (should be 73,707)
#'
#' @section GAF transcript names:
#'
#' Transcript names correspond to the UCSCgene.Dec2009 release. Transcript
#' names are unique.
#'
#' @section Errors:
#'    These errors are fatal and will terminate processing.
#'
#' \describe{
#'    \item{
#'       \command{Unsafe character in GAF filename!}
#'    }{
#'       An invalid characters was passed as part of the GAF filename. This is
#'       important as the filename is used in a system command as a parameter
#'       and could be used for command injection.
#'    }
#'    \item{
#'       \command{Can't find the specified GAF: "\var{file}"}
#'    }{
#'       The specified GAF doesn't seem to exist on the file system.
#'       Probably have the name wrong or are using a relative name from the
#'       wrong directory, but could also be that permissions are hiding it.
#'    }
#'    \item{
#'       \command{Unsafe character in output transcriptModel filename!}
#'    }{
#'       An invalid characters was passed as part of the output filename. This
#'       is important as the filename is used in a system command as a parameter
#'       and could be used for command injection.
#'    }
#'    \item{
#'       \command{Output file already exists; use force= TRUE to overwrite: "\var{file}"}
#'    }{
#'       The specified GAF transcript extract output file already exists. You
#'       probably don't want to overwrite it. However, you can set \code{force=
#'       TRUE} to allow this. It will still generate a warning.
#'    }
#' }
#'
#' @section Warnings:
#'
#' \describe{
#'    \item{
#'       \command{Forcing overwrite of output file: "\var{file}"}
#'    }{
#'       Just letting you know an existing file is actually being overwritten.
#'       This won't happen unless explicitly allowed by setting \code{force=
#'       TRUE}). Having a warning allows distinguishing between the cases where
#'       an overwrite occurred vs those where one was allowed but did not occur.
#'    }
#'    \item{
#'       Various warnings from failed system commands
#'    }{
#'       System commands are used for several things in this function. If they
#'       fail, error messages are returned as warnings.
#'    }
#'
#' }
#'
#' @section Todo:
#'    \itemize{
#'       \item{Add test for corner case - one exon transcript}
#'    }
#'
#' @examples
#' \dontrun{
#'
#' # Extract transcripts to default output file
#' count <- extractTranscriptModels( 'path/to/GAF' )
#'
#' # Same, with all defaults made explicit
#' count <- extractTranscriptModels(
#'    gaf= 'path/to/GAF', outFile= 'path/to/GAF.transcriptModels', force= FALSE
#' )
#'
#' # Extract transcripts to gaf.transcripts in run directory
#' count <- extractTranscriptModels( 'path/to/GAF', 'gaf.transcripts' )
#'
#' # Overwrite outFile if it exists (here using the default name)
#' count <- extractTranscriptModels( 'path/to/GAF', force= TRUE )
#' }
#'
#' @export
extractTranscriptModels <- function(
   gaf,
   outFile= paste0(gaf, ".transcriptModels"),
   force= FALSE
) {

   # Using filename as a parameter in a system command so need to be careful.
   if (! is.safeFileName( gaf )) {
      stop( "Unsafe character in GAF filename!" );
   }
   if (! file.exists(gaf)) {
      stop( "Can't find the specified GAF: \"", gaf, "\"")
   }
   if (! is.safeFileName( outFile )) {
      stop( "Unsafe character in output transcriptModel filename!" );
   }
   if (file.exists(outFile)) {
      if (! force) {
         stop( "Output file already exists; use force= TRUE to overwrite: \"", outFile, "\"")
      } else {
         warning( "Forcing overwrite of output file: \"", outFile, "\"");
         unlink(outFile);
      }
   }

   # Load header line into output
   cmd = paste(c("head -n1 ", gaf, " > ", outFile ), collapse="");
   system( cmd, intern=TRUE );

   # Extract gene models.
   cmd = paste( c( "grep transcript ", gaf,
                   " | grep UCSCgene.Dec2009.fa",
                   " | grep -v calculated",
                   " >> ", outFile ),
                collapse= ""
   );
   system( cmd, intern=TRUE );

   # Count returned lines
   cmd = paste(c("wc -l ", outFile), collapse="");
   cmdOut = system( cmd, intern=TRUE );
   lineCount = as.integer( sub( " .+", "", sub( "^[ ]+", "", cmdOut )));

   return(lineCount - 1);
}

#' Load the gene models file into a data frame.
#'
#' @param file The name of the gene models file extracted from the GAF
#'
#' @return Data frame with the data from the GAF gene models:
#'
#'   \tabular{ll}{
#'       \code{$gene} \tab Feature name\cr
#'       \code{$chr}  \tab chromosome as chr#, chr##, chrX, chrY, chrM, or chrM_rCRS,
#'                         or possibly a contig name, like 'GL000228.1'.\cr
#'       \code{$strand} \tab strand as *, +, or -\cr
#'       \code{$gstart} \tab first base of gene\cr
#'       \code{$gend} \tab last base of gene\cr
#'       \code{$exon} \tab Exon number\cr
#'       \code{$start} \tab first base of exon\cr
#'       \code{$end} \tab Last base of exon\cr
#'       \code{$length} \tab length of the exon in bases\cr
#'    }
#'
#' @export
loadGeneModels <- function ( file ) {
   .Deprecated(loadGafModels)
   return(loadGafModels(file))
#    geneModelData <- read.delim( file, stringsAsFactors = FALSE );
#
#    allGenes <- sub("[|].+$", "", geneModelData$Gene);
#    chr <- sub(":.+:[-+]", "", geneModelData$CompositeCoordinates);
#    exons <- sub( "^.+:", "", sub( ":[-+]$", "", geneModelData$CompositeCoordinates ));
#    gstart <- sub("-.+","", exons);
#    gend <- sub(".+[-]","", exons);
#    exons <- strsplit(exons,",")
#    strand <- sub("^.+:.+:","", geneModelData$CompositeCoordinates);
#    getStartsF <- function(x) {sub("-.+$", "", x )}
#    getEndsF <- function(x) {sub("^[^-]+-", "", x )}
#    getNumF <- function(x) {1:length(x)}
#    modelsMatrix <- do.call( rbind, mapply( cbind,
#                                            allGenes, chr, strand, gstart, gend, lapply( exons, getNumF ), lapply( exons, getStartsF ), lapply( exons, getEndsF )
#    ));
#    geneModels <- data.frame(
#       gene= modelsMatrix[,1],
#       chr= modelsMatrix[,2],
#       strand= modelsMatrix[,3],
#       gstart=  as.integer(modelsMatrix[,4]),
#       gend=  as.integer(modelsMatrix[,5]),
#       exon= as.integer(modelsMatrix[,6]),
#       start= as.integer(modelsMatrix[,7]),
#       end= as.integer(modelsMatrix[,8]),
#       stringsAsFactors= FALSE
#    );
#    geneModels$length <- as.integer(geneModels$end - geneModels$start + 1);
#
#    return(geneModels);
}

#' Load a GAF models file into a per exon data frame.
#'
#' @param file The file of models extracted from the GAF.
#'
#' @return Returns a per-exon data frame of the models extracted from the GAF.
#'   Most models (e.g. genes or transcripts) are composed of multiple exons with
#'   each features-exon in separate row. For convenience, the table replicates a
#'   lot of information that is the same for every exon in a model as a string,
#'   i.e columns \code{gene}, \code{chr}, \code{strand}, \code{gstart}, and
#'   \code{gend}).
#'
#'   The data frame contains the following columns
#'   \tabular{ll}{
#'      \code{gene}   \tab Feature name - gene name, transcript id, etc.\cr
#'      \code{chr}    \tab Chromosome, as chr#, chr##, chrX, chrY, chrM, or chrM_rCRS
#'                         or possibly a contig name, like 'GL000228.1'.\cr
#'      \code{strand} \tab Strand, as *, +, or -\cr
#'      \code{gstart} \tab Genomic coordinate of first base of feature\cr
#'      \code{gend}   \tab Genomic coordinate of last base of feature\cr
#'      \code{exon}   \tab Exon number within feature\cr
#'      \code{start}  \tab Genomic coordinate of first base of exon\cr
#'      \code{end}    \tab Genomic coordinate of last base of exon\cr
#'      \code{length} \tab length of the exon in bases\cr
#'   }
#'
#' @section Todo:
#'   \itemize{
#'      \item{change column named "gene" to "name"}
#'      \item{add column "gafId"}
#'      \item{consider an alias column with | separated names}
#'      \item{change gstart and gend to mstart and mend}
#'   }
#'
#' @section Errors:
#'    These errors are fatal and will terminate processing.
#'
#'    \describe{
#'       \item{
#'          \command{Can't find the specified GAF extract file: "\var{file}"}
#'       }{
#'          The specified GAF doesn't seem to exist on the file system.
#'          Probably have the name wrong or are using a relative name from the
#'          wrong directory, but could also be that permissions are hiding it.
#'       }
#'    }
#'
#' @examples
#'
#' \dontrun{
#'
#' extractGeneModels( 'TCGA.hg19.June2011.gaf' )
#' gafGeneModelsFile <- 'TCGA.hg19.June2011.gaf.geneModels'
#' geneExonsDF <- loadGafModels( gafTranscriptModelsFile )
#'
#' extractTranscriptModels( 'TCGA.hg19.June2011.gaf' )
#' gafTranscriptModelsFile <- 'TCGA.hg19.June2011.gaf.transcriptModels'
#' transcriptExonsDF <- loadGafModels( gafTranscriptModelsFile )
#' }
#'
#' @export
loadGafModels <- function ( file ) {

   if (! file.exists(file)) {
      stop( "Can't find the specified GAF extract file: \"", file, "\"")
   }

   df <- read.delim( file, stringsAsFactors = FALSE );

   id <- sub("[|].+$", "", df$FeatureID);
   chr <- sub(":.+:[-+]", "", df$CompositeCoordinates);
   exons <- sub( "^.+:", "", sub( ":[-+]$", "", df$CompositeCoordinates ));
   gstart <- sub("-.+","", exons);
   gend <- sub(".+[-]","", exons);
   exons <- strsplit(exons,",")
   strand <- sub("^.+:.+:","", df$CompositeCoordinates);
   getStartsF <- function(x) {sub("-.+$", "", x )}
   getEndsF <- function(x) {sub("^[^-]+-", "", x )}
   getNumF <- function(x) {1:length(x)}

   modelsMatrix <- do.call(
      rbind, mapply( cbind,
         id, chr, strand, gstart, gend,
         lapply( exons, getNumF ),
         lapply( exons, getStartsF ),
         lapply( exons, getEndsF )
   ));
   geneModels <- data.frame(
      gene= modelsMatrix[,1],
      chr= modelsMatrix[,2],
      strand= modelsMatrix[,3],
      gstart=  as.integer(modelsMatrix[,4]),
      gend=  as.integer(modelsMatrix[,5]),
      exon= as.integer(modelsMatrix[,6]),
      start= as.integer(modelsMatrix[,7]),
      end= as.integer(modelsMatrix[,8]),
      stringsAsFactors= FALSE
   );
   geneModels$length <- as.integer(geneModels$end - geneModels$start + 1);

   return(geneModels);
}

# Get the list of all gene names in the model data frame
getAllGeneNames <-function( geneModelsDF ) {
   ###
   #     PARAM geneModelsDF the data frame returned by laodGeneModels
   ###
   #     Returns a vector of gene names (guaranteed to be unique, even if
   # repeated in the data frame.
   ###
   return(unique(geneModelsDF$gene));
}
