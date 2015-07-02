

# Verifies that a fileName is safe for use in shell commands.
# Needed by extractGeneModels
is.safeFileName <- function( names,
                             exists= NULL, weak=FALSE,
                             okSpace=TRUE, okTilde=TRUE,
                             okBackslash=FALSE, okComma=FALSE, okColon=FALSE
) {
   ###
   # Validates file names, including for use in system commands. Uses both
   # blacklist and whitelist filtering, allowing several characters to be
   # optionally allowed or banned. Optionally may also require a filename
   # to exist or not to exist. Returns a vector of TRUE for validated file
   # names or FALSE if fail.
   #
   # Some characters are never allowed:
   #     "`", "$", "(", ")", "|", ";", "&", ">", "<"
   # Some characters are always allowed:
   #     "A" - "Z", "a" - "z", "0" - "9"
   #     "_", ".", "-", "/"
   # Some characters selectably bad are allowed by default:
   #     " ", "~" (~ allowed anywhere).
   # Some characters selectably bad are dis-allowed by default:
   #     "\", ",", ":"
   # Additionally, no file name may start with "+" or "-", nor may
   # it contain " -" of " +". This prevents injecting new options.
   ###
   #     PARAM: names - vector of file names to check.
   #     PARAM: weak=FALSE - If set TRUE, only check for known problems
   # (i.e. only do blacklist validation). If FALSE (the default)
   # then only allows a limited range of filenames (i.e. also
   # do whitelist validation.)
   #     PARAM: okSpace= TRUE - Allows " " by default.
   #     PARAM: okTilde= TRUE - Allows "~" by default (anywhere!).
   #     PARAM: okBackslash= FALSE - No "\" allowed by default.
   #     PARAM: okComma=     FALSE - No "," allowed by default.
   #     PARAM: okColon=     FALSE - No ":" allowed by default.
   #     PARAM: exists= NULL - If set TRUE, file must exist. If set FALSE
   # file may not exist. (Note, don't rely on this as only error check for
   # reading/writting file as that introduces a race condition).
   ###
   #     RETURNS: A vector of booleans the same length as the input
   # names, TRUE if the file is safe, FALSE if it fails any test.
   ###

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

test.is.safeFileName <- function() {
   goodFileNames=c("x", "a.txt", "3", "/me/and/you", "this_or-that");
   got <- is.safeFileName( goodFileNames );
   if (! all( got )) {
      print(got);
      stop("is.safeFileName() Failed good file names test!");
   }
   got <- is.safeFileName( goodFileNames, weak=TRUE );
   if (! all( got )) {
      print(got);
      stop("is.safeFileName() Failed good file names weak test!");
   }

   badFileNames=c("`ls`", "$(USER)", "$", "(", ")", "f.txt;ls", "> x.txt", "<bad.txt", "|ls", "&ls", "+opt", "-opt", "arg +opt", "arg -opt");
   got <- is.safeFileName( badFileNames );
   if ( any( got )) {
      print(got);
      stop("is.safeFileName() failed bad file names test!");
   }
   got <- is.safeFileName( badFileNames, weak=TRUE );
   if ( any( got )) {
      print(got);
      stop("is.safeFileName() failed bad file names weak test");
   }

   defaultsCheckNames=c("blank ok", "any~ok", "no,", "no:", "no\\", "no\\\\");
   expect = c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE);
   got <- is.safeFileName( defaultsCheckNames );
   if ( ! identical( got, expect )) {
      print(got);
      print(expect);
      stop("is.safeFileName() failed defaults file names test");
   }
   got <- is.safeFileName( defaultsCheckNames, weak=TRUE );
   if ( ! identical( got, expect )) {
      print(got);
      print(expect);
      stop("is.safeFileName() failed defaults file names weak test");
   }

   expect = !expect;
   got <- is.safeFileName( defaultsCheckNames, okSpace=FALSE, okTilde=FALSE, okComma=TRUE, okColon=TRUE, okBackslash=TRUE );
   if ( ! identical( got, expect )) {
      print(got);
      print(expect);
      stop("is.safeFileName() failed anti-defaults file names test");
   }
   got <- is.safeFileName( defaultsCheckNames, okSpace=FALSE, okTilde=FALSE, okComma=TRUE, okColon=TRUE, okBackslash=TRUE, weak=TRUE );
   if ( ! identical( got, expect )) {
      print(got);
      print(expect);
      stop("is.safeFileName() failed anti-defaults weak file names test");
   }

   name = "Not banned for weak%*@'=[]{}?!";
   got <- is.safeFileName( name, weak=TRUE );
   if ( ! got) {
      stop("is.safeFileName() failed weak file names test");
   }
   got <- is.safeFileName( name, weak=TRUE, okSpace=FALSE );
   if ( got ) {
      stop("is.safeFileName() failed weak file name negative test");
   }

   return( TRUE );
}

#' Extract gene models from the gaf
#'
#' This extracts a tab delimited file of all canonical
#' exon gene models from the GAF. The cannonical model is intended to
#' be at least a common starting point for a simplified definition of "the gene
#' called X", based on the union of all transcript. By default each model in the
#' output file (row) can be refereed to uniquely by gene name.

#' This requires the "grep" program, so only works on linux/mac systems. The Gaf
#' file version this works with is the version used for the TCGA RNAseq
#' expression data files. It is available for download at the NCI uncompressed:
#' \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf}{TCGA.hg19.June2011.gaf}
#' or gzipped:
#' \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf.gz}{TCGA.hg19.June2011.gaf.gz}
#'
#' @export
#' @param gafFile The full-path name of the GAF file [REQ]
#' @param outFile The output file name, by default the same as the
#'    gafFile, with ".geneModels" appended. Will not overwrite unless force= is
#'    set true
#' @param force If the output file exists, overwrite (generates
#'   a warning).
#' @param uniqueGene The gaf file contains several versions of
#' some genes, these are annotated as  "1ofN", "2ofN", etc. By default, keeps
#' only the "...1of" variant. Also, SLC35E2 is present twice, but without a
#' numbered annotation. Keeps only the the larger (encompasing) version of
#' SLC35E2. Set false to keep all versions of all genes (gene name will not
#' then be a unique key.). Also implies skipUnknownGene
#' @param skipUnknownGene Drops unknown genes, those whose names
#' in the gaf are "?". There are 32 in the version of the GAF this targets.
#' Set false to keep the unknown genes. Note, this setting is over-ridden to
#' TRUE (skipping these genes) when uniqueGene= TRUE as all these genes have
#' the same gene name, "?".
#'
#' @return The main output is the gaf extract file, which is just the
#' gene models and by default just the genes with unique names. See
#' \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/other/GAF/GAF_bundle_Feb2011/docs/GAF_v2_file_description.docx}{GAF_v2_file_description.docx}
#' However, this function also returns a list with some summary info about the
#' gaf and the generated geneModels file:
#'
#'     info$gaf - The gaf file name used as a parameter (relative file name)
#'     info$gaf_real - The absolute full path filename to the input gaf file
#'     info$gaf_md5 - The md5 checksum of the gaf file, as a string
#'     info$uniqueGene - The uniqueGene parameter setting used.
#'     info$skipUnknownGene - The skipUnknownGene parameter setting used.
#'     info$gaf_lines - The number of lines in the gaf file
#'     info$gaf_models - The number of models in the gaf file
#'     info$gaf_models_unique - The number of unique models in the gaf file
#'     info$gaf_extract - The geneModels file name, based on the input gaf
#'     info$gaf_extract_real - The absolute full path filename of the geneModels file
#'     info$gaf_extract_md5 - The md5 checksum of the geneModels file, as a # string
extractGeneModels <- function( gafFile,
                               outFile= paste0(gafFile, ".geneModels"), force=FALSE,
                               uniqueGene= TRUE, skipUnknownGene= TRUE
) {


   # Output line counts include 1 header line in all files.

   # Using filename as a parameter in a system command so need to be careful.
   if (! is.safeFileName( gafFile )) {
      stop( "Unsafe character in gaf file name!" );
   }
   if (! file.exists(gafFile)) {
      stop( 'Can\'t find the specified gaf file: \'', gafFile, '\'')
   }
   if (! is.safeFileName( outFile )) {
      stop( "Unsafe character in output geneModel file name!" );
   }
   if (file.exists(outFile)) {
      if (! force) {
         stop( "Output file already exists; use force=TRUE to overwrite: '", outFile, "'")
      } else {
         warning( "Forcing overwrite of output file '", outFile, "'");
         unlink(outFile);
      }
   }

   if (uniqueGene && ! skipUnknownGene) {
      warning(" uniqueGene=TRUE sets ignoreUnknownGene=TRUE")
      skipUnknownGene=TRUE;
   }

   info = list();
   info$gaf = gafFile;
   info$gaf_real = normalizePath(gafFile);
   info$gaf_md5 = tools::md5sum(gafFile);
   info$uniqueGene = uniqueGene;
   info$skipUnknownGene = skipUnknownGene;

   # Get number of lines in original gaf
   cmd = paste( c( "wc -l ", gafFile ), collapse="" );
   got = system( cmd, intern=TRUE );
   info$gaf_lines = as.integer( sub( " .+", "", sub( "^[ ]+", "", got )));

   geneModelsFile = outFile;
   if (uniqueGene) {
      geneModelsFile = tempfile( "TEMP_gafFileExtractGenome" );
      if (file.exists(geneModelsFile)) {
         if (! force) {
            stop( "Temporary file exists, delete or use force: ", geneModelsFile );
         } else {
            warning( "Forcing overrite of pre-existing temporary file: ", geneModelsFile );
         }
         unlink( geneModelsFile );
      }
   }

   # Load header line into output
   cmd = paste(c("head -n1 ", gafFile, " > ", geneModelsFile ), collapse="");
   system( cmd );

   dropUnknownGene = " | grep -v -E \"\t[?][|][0-9]\"";
   if (! skipUnknownGene) {
      dropUnknownGene = ""
   }

   # Extract gene models.
   cmd = paste( c( "grep genome ", gafFile,
                   " | grep calculated",
                   " | grep chr",
                   dropUnknownGene,
                   " >> ", geneModelsFile ),
                collapse= "" );
   system( cmd );
   cmd = paste(c("wc -l ", geneModelsFile), collapse="");
   got = system( cmd, intern=TRUE );
   info$gaf_models = as.integer( sub( " .+", "", sub( "^[ ]+", "", got )));

   # If want truely unique models, need to ignore the 2of3, etc lines and the
   # duplicate SLC35E2 gene. Will keep the [1 of x] and the larger (encompasing)
   # version of SLC35E2.
   if (uniqueGene) {
      # As grep for what doesn not match, will keep header line automatically
      cmd = paste(
         c( "grep -v -E \"[02-9]of|1[0-9]of|SLC35E2.9906\" ", geneModelsFile, " > ", outFile ),
         collapse="" );
      system( cmd );

      # Done with tempCalculate intermediate; calculate line count in out file
      unlink(geneModelsFile)
      cmd = paste(c("wc -l ", outFile), collapse="");
      got = system( cmd, intern=TRUE );
      info$gaf_models_unique = as.integer( sub( " .+", "", sub( "^[ ]+", "", got )));
   }

   info$gaf_extract= outFile;
   info$gaf_extract_real= normalizePath(outFile);
   info$gaf_extract_md5 = tools::md5sum(outFile);
   return( info )
}

demo.extractGeneModels <- function( gaf ) {
   info <- extractGeneModels( gaf, force=TRUE );
   return( info );
}

#' Extract transcript models from the gaf
#'
#' This extracts a tab delimited file of all canonical transcript models from the
#' GAF. Each model in the output file (row) can be refereed to
#' uniquely by transcript name.
#'
#' @param gafFile The full-path name of the GAF file [REQ]
#' @param outFile The output file name, by default the same as the gafFile, with
#'   ".geneModels" appended. Will not overwrite unless force= is set true
#' @param force If the output file exists, overwrite (generates a warning).
#'
#' @return The number of transcripts in the output file (should be 73,707)
#'
#' @export
extractTranscriptModels <- function(
   gafFile,
   outFile= paste0(gafFile, ".transcriptModels"),
   force=FALSE
) {

   # Using filename as a parameter in a system command so need to be careful.
   if (! is.safeFileName( gafFile )) {
      stop( "Unsafe character in gaf file name!" );
   }
   if (! file.exists(gafFile)) {
      stop( 'Can\'t find the specified gaf file: \'', gafFile, '\'')
   }
   if (! is.safeFileName( outFile )) {
      stop( "Unsafe character in output transcriptModel file name!" );
   }
   if (file.exists(outFile)) {
      if (! force) {
         stop( "Output file already exists; use force=TRUE to overwrite: '", outFile, "'")
      } else {
         warning( "Forcing overwrite of output file '", outFile, "'");
         unlink(outFile);
      }
   }

   # Load header line into output
   cmd = paste(c("head -n1 ", gafFile, " > ", outFile ), collapse="");
   system( cmd );

   # Extract gene models.
   cmd = paste( c( "grep transcript ", gafFile,
                   " | grep UCSCgene.Dec2009.fa",
                   " | grep -v calculated",
                   " >> ", outFile ),
                collapse= ""
   );
   system( cmd );

   # Count returned lines
   cmd = paste(c("wc -l ", outFile), collapse="");
   cmdOut = system( cmd, intern=TRUE );
   lineCount = as.integer( sub( " .+", "", sub( "^[ ]+", "", cmdOut )));

   return(lineCount - 1);
}

#' Load the gene models file into a data frame.
#'
#' @param file The name of the gene models file extracted from the gaf file
#'
#' @return Data frame with the data from the gaf gene models:
#'
#'       $gene = gene name
#'       $chr = chromosome as chr#, chr##, chrX, chrY, chrM, or chrM_rCRS
#'       $strand= strand as *, +, or -
#'       $gstart=  first base of gene
#'       $gend=  last base of gene
#'       $exon= Exon number
#'       $start= first base of exon
#'       $end= Last base of exon
#'       $length= length of the exon in bases
#'
#' @export
loadGeneModels <- function ( file ) {

   geneModelData <- read.delim( file, stringsAsFactors = FALSE );

   allGenes <- sub("[|].+$", "", geneModelData$Gene);
   chr <- sub(":.+:[-+]", "", geneModelData$CompositeCoordinates);
   exons <- sub( "^.+:", "", sub( ":[-+]$", "", geneModelData$CompositeCoordinates ));
   gstart <- sub("-.+","", exons);
   gend <- sub(".+[-]","", exons);
   exons <- strsplit(exons,",")
   strand <- sub("^.+:.+:","", geneModelData$CompositeCoordinates);
   getStartsF <- function(x) {sub("-.+$", "", x )}
   getEndsF <- function(x) {sub("^[^-]+-", "", x )}
   getNumF <- function(x) {1:length(x)}
   modelsMatrix <- do.call( rbind, mapply( cbind,
                                           allGenes, chr, strand, gstart, gend, lapply( exons, getNumF ), lapply( exons, getStartsF ), lapply( exons, getEndsF )
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

#' Load a gaf feature file into an exon data frame.
#'
#' @param file The name of the transcrpt models file extracted from the gaf file
#'
#' @return Data frame with the data from the gaf transcrpt models:
#' /describe{
#'    {$gene}{gene name}
#'    {$chr}{chromosome as chr#, chr##, chrX, chrY, chrM, or chrM_rCRS}
#'    {$strand}{strand as *, +, or -}
#'    {$gstart}{first base of gene}
#'    {$gend}{last base of gene}
#'    {$exon}{Exon number}
#'    {$start}{first base of exon}
#'    {$end}{Last base of exon}
#'    {$length}{length of the exon in bases}
#' }
#' @export
loadGafModels <- function ( file ) {

   df <- read.delim( file, stringsAsFactors = FALSE );

   id <- sub("[|].+$", "", df$FeatureID);
   chr <- sub(":.+:[-+]", "", df$CompositeCoordinates);
   exons <- sub( "^.+:", "", sub( ":[-+]$", "", df$CompositeCoordinates ));
   start <- sub("-.+","", exons);
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
getAllGeneNames <-function( geneModesDF ) {
   ###
   #     PARAM geneModelsDF the data frame returned by laodGeneModels
   ###
   #     Returns a vector of gene names (guaranteed to be unique, even if
   # repeated in the data frame.
   ###
   return(unique(geneModels$gene));
}
