isDebug = TRUE

getCandidateFusions <- function (
   fusionDF,
   sample, gene1, gene2
) {
   ###
   #      Select out the candidate fusions for the given samples as a new data.
   # For each samples, include any fusion involving either gene. May not care
   # about a specific sample-gene1-gene2 combo, but that is for somwthing else
   # to filter.
   ###
   matchSample = fusionDF$sample %in% sample;
   endInG1 = fusionDF$gene1 %in% c(gene1,gene2);
   endInG2 = fusionDF$gene2 %in% c(gene1,gene2);
   selectVec = matchSample & (endInG1 | endInG2);

   df <- fusionDF[ selectVec,  ];
   return(df);
}

getMapSpliceCohortFusionData <- function( file, all.columns=FALSE ) {
   ###
   # Loads a mapplice cohort fusion file.
   ###
   #     PARAM file: path to the data file to load
   #     PARAM all.columns: Set TRUE to load all columns from data file. By
   # default (all.columns=FALSE) only loads:
   #   "D_gene", renamed to "gene1";
   #   "A_gene", renamed to "gene2";
   #   "D_end", renamed to "gene1pos";
   #   "A_start", renamed to "gene2pos";
   #   "tcga_id", renamed to "sample";
   # These are renamed regardless of all.columns= setting.
   ###
   #      RETURNS: A data frame with the renamed columns, plus rest if
   # all.columns=TRUE. An "id" column is added at left (row number from
   # original list)

   ###
   if (! file.exists( file )) {
      stop( "No such file: ", file );
   }

   # Get headers separately as there are more columns than headers.
   headers= read.delim( file, nrow= 1, header=FALSE, stringsAsFactors= FALSE );
   df = read.delim( file, skip=1, header= FALSE, stringsAsFactors= FALSE);
   names(df) = headers;
   id = as.integer(row.names(df))
   df <- cbind(id, df)

   # There has to be a better way to do this...
   names(df)[ names(df) == "D_gene"  ] <- "gene1";
   df$gene1 <- sub(",$", "", df$gene1)
   names(df)[ names(df) == "A_gene"  ] <- "gene2";
   df$gene2 <- sub(",$", "", df$gene2)
   names(df)[ names(df) == "D_end"   ] <- "gene1pos";
   names(df)[ names(df) == "A_start" ] <- "gene2pos";
   names(df)[ names(df) == "tcga_id" ] <- "sample";

   if (! all.columns) {
      df <- df[ ,c("id","sample","gene1","gene2","gene1pos","gene2pos" )];
   }

   return(df)
}

filterFusions <- function( fusions, sample=c(), gene1=c(), gene2=c(), genePairing= "or") {
   ###
   # Sub-select fusion lines from a list of fusions based on sample name or gene.
   # If you need more complex filtering,
   # just do it yourself. This is a convienience function.
   ###
   #    PARAM: sample = list of sample names. If given, only fusions for the
   # samples with exactly mathcing sample names will be kept. Gene filtering
   # will select only from these samples.
   #    PARAM: gene1 = list of gene names to filter fusions based on the
   # upstream gene name.
   #    PARAM: gene2 = list of gene names to filter fusions based on the
   # downstream gene name.
   #    PARAM: genePairing= "or". Logic to use if both gene1 and gene2 given.
   # "or": If fusion gene1 in gene1 param or fusion gene2 in gene2 param, keep the fusion.
   # "and": If fusion gene1 in gene1 param and fusion gene2 in gene2 param, keep the fusion.
   # "pair" Keep fusion if param gene1[i] = fusion gene1 and param gene2[i] =
   # fusion gene2, for any gene param i (gene1 and gene2 must be same length)

   df <- fusions;  # all fusions
   if (length(sample) > 0) {
      sampleSelect = df$sample %in% sample;
      df <- df[sampleSelect, ];
   }

   # df contains only fusions for kept samples
   if (length(gene1) > 0 && length(gene2) > 0) {
      # Need to determine how to report for both columns
      if (genePairing == "or") {
         gene1Select = df$gene1 %in% gene1;
         gene2Select = df$gene2 %in% gene2;
         df <- df[gene1Select | gene2Select, ];
      } else if (genePairing == "and" ) {
         gene1Select = df$gene1 %in% gene1;
         gene2Select = df$gene2 %in% gene2;
         df <- df[gene1Select & gene2Select, ];
      } else if (genePairing == "pair" ) {
         if ( length(gene1) != length(gene2) ) {
            stop( "Can't use genePairing = 'pair' if gene1 and gene2 have different lengths.");
         }
         select = rep(FALSE, nrow(fusions));
         for (i in 1:length(gene1)) {
            select = (df$gene1 == gene1[i] & df$gene2 == gene2[i]) | select
         }
         df <- df[select, ];
      } else {
         stop( paste0( "unknown gene pairing", genePairing))
      }
   } else if (length(gene1) > 0) {
      # Only reporting gene1, as can't get here if also have gene2 specified
      gene1Select = df$gene1 %in% gene1;
      df <- df[gene1Select, ];
   } else if (length(gene2) > 0) {
      # Only reporting gene2, as can't get here if also have gene1 specified
      gene2Select = df$gene2 %in% gene2;
      df <- df[gene2Select, ];
   }
   return(df)
}


#    selectedFusions <- fusionData[
#       fusionData$tcga_id == "TCGA-2A-A8VL-01A-21R-A37L-07" & (
#          (fusionData$D_gene == gene1 & fusionData$A_gene == gene2)
#          | (fusionData$D_gene == gene2 & fusionData$A_gene == gene1)
#       ),
#       c("tcga_id", "D_gene", "A_gene", "D_end", "A_start")
#       ];
#    fusions = list();
#    fusions[[gene1]] <- c( selectedFusions[selectedFusions$D_gene == gene1, "D_end"],
#                           selectedFusions[selectedFusions$A_gene == gene1, "A_start"] );
#    fusions[[gene2]] <- c( selectedFusions[selectedFusions$D_gene == gene2, "D_end"],
#                           selectedFusions[selectedFusions$A_gene == gene2, "A_start"] );
#
#    return(fusions);
# }

# formatMapspliceFusions <- function ( fusions ) {
#    names(fusions)[ names(fusions) == "D_gene"  ] <- "gene1";
#    names(fusions)[ names(fusions) == "A_gene"  ] <- "gene2";
#    names(fusions)[ names(fusions) == "D_end"   ] <- "gene1pos";
#    names(fusions)[ names(fusions) == "A_start" ] <- "gene2pos";
#    names(fusions)[ names(fusions) == "tcga_id" ] <- "sample";
#    names(fusions)[ names(fusions) == "X"       ] <- "id";
#
#    return (fusions);
# }

do.allPlots <- function( fusions, normalizedCohortExpression ) {
   plotCount = 0;
   fusionLineCount = 0;
   # Iterate through each sample in the fusions data frame
   samples = unique(fusions$sample)
   for (sample in samples) {
      if (isDebug) { print("sample"); print(sample) }
      sampleSelectedFusions = fusions[ fusions$sample == sample, ]

      # For the fusions for this sample, iterate through each different gene 1
      # Note: expect intergenic fusions to have a gene name of "-"
      genes1 = unique(sampleSelectedFusions$gene1)
      for (gene1 in genes1) {
         if (isDebug) { print("gene1"); print(gene1) }
         sampleGene1SelectedFusions = sampleSelectedFusions[ sampleSelectedFusions$gene1 == gene1, ]

         # For the fusions for this sample and this gene1, itterate through each different gene 2
         genes2 = unique(sampleGene1SelectedFusions$gene2)
         for (gene2 in genes2) {
            if (isDebug) { print("gene2"); print(gene2) }
            selectedFusions = sampleGene1SelectedFusions[ sampleGene1SelectedFusions$gene2 == gene2, ]
            # selectedFusions now contains only the fusion pair to plot, possible
            # involving an intergenic gene.

            # Select the exon data columns and the expression data for sample from
            # the normalizedCohortExpression data frame for genes 1 and 2
            expressionRowSelect = normalizedCohortExpression$gene %in% c(gene1, gene2)
            expressionColumnNames = c(names(normalizedCohortExpression)[1:9], sample)
            selectedExpression = normalizedCohortExpression[expressionRowSelect, expressionColumnNames];

            # Plot one gene if either gene1 or gene 2 is "-", otherwise plot pair.
            # Fu
            if (gene1 == '-') {
               if (isDebug) { print("Plotting one gene (gene2)"); }
               plot.fusionExpressionOne( gene2, selectedFusions$gene2pos, sample, selectedExpression )
            } else if (gene2 == '-') {
               if (isDebug) { print("Plotting one gene (gene1)"); }
               plot.fusionExpressionOne( gene1, selectedFusions$gene1pos, sample, selectedExpression )
            } else {
               if (isDebug) { print("Plotting two genes"); }
               plot.fusionExpressionPair( gene1, gene2,
                                          selectedFusions$gene1pos, selectedFusions$gene2pos,
                                          sample, selectedExpression )
            } # end if/else (select how to plot data)
            plotCount = plotCount + 1;
            fusionLineCount = fusionLineCount + nrow(selectedFusions)
         } # end loop over unique fusions$gene2
      } # end loop over unique fusions$gene1
   } # end loop over unique fusions$sample
   return( list( plotCount= plotCount, fusionLineCount= fusionLineCount ) );
}

plotGeneExpression <- function( samples, genes, cohortExonExpression) {
   ###
   #     DESC: For each pair of sample and gene, plot the expression for that
   # gene in that sample.
   #     PARAM: samples = vector of sample names giving the expression column
   # from the cohortExonExpression to plot. Must be the same length as the
   # gene parameter
   #     PARAM: genes = vector of gene names giving the exon rows from the
   # cohortExonExpression to plot. Must be the same length as the sample
   # parameter
   #     PARAM: cohortExonExpression = cohort expression data frame, with the
   # first 9 columns giving the gene model with one exon per row, and the
   # remaining columns giving exon expression for different samples with their
   # headings the sample names.
   #     RETURNS: nothing
   ###
   if (length(genes) != length(samples)) { stop("Lengths of gene and sample vectors must be the same."); }

   for (plotNum in 1:length(samples)) {
      selectedCols <- c( names(cohortExonExpression)[1:9], samples[plotNum]);
      selectedRows <- c( cohortExonExpression$gene == genes[plotNum]);
      selectedExpression <- cohortExonExpression[selectedRows, selectedCols];
      plot.fusionExpressionOne( genes[plotNum], NULL, samples[plotNum], selectedExpression )
   }

   return( NULL );
}

do.plots <- function( fusionsDf, samples, geneOrderList, normalizedCohortExpressionDf ) {
   plotCount = 0;
   fusionCount = 0;
   singleFusionCount = 0;
   if (length(intersect(c("-"), geneOrderList)) > 0 ) { stop("Gene list may not contain the intergenic symbol \'-\'"); }
   for (sample in samples) {
      print(paste0("Sample ", sample));
      sampleFusionDf <- fusionsDf[fusionsDf$sample == sample,];
      normalizedSampleExpressionDf <- normalizedCohortExpressionDf[,c(names(normalizedCohortExpressionDf)[1:9], samples)]
      for (geneOfInterest in geneOrderList) {
         print("Main gene:");
         print(geneOfInterest);
         sampleFusionGeneDf <- sampleFusionDf[sampleFusionDf$gene1 == geneOfInterest | sampleFusionDf$gene2 == geneOfInterest,]
         print("sampleFusionGeneDf:");
         print(sampleFusionGeneDf);
         altGenes <- unique( c(sampleFusionGeneDf[sampleFusionGeneDf$gene1 == geneOfInterest,"gene2"],
                               sampleFusionGeneDf[sampleFusionGeneDf$gene2 == geneOfInterest,"gene1"] ));
         print("altGenes:");
         print(altGenes);
         for (altGene in altGenes) {
            print("altGene:");
            print(altGene);
            sampleFusionGeneAltDf <- sampleFusionGeneDf[sampleFusionGeneDf$gene1 == altGene | sampleFusionGeneDf$gene2 == altGene,]
            print("sampleFusionGeneAltDf:");
            print(sampleFusionGeneAltDf);
            fusionsOfInterest <- c(sampleFusionGeneAltDf[sampleFusionGeneAltDf$gene1 == geneOfInterest,"gene1Pos"],
                                   sampleFusionGeneAltDf[sampleFusionGeneAltDf$gene2 == geneOfInterest,"gene2Pos"])
            print("fusionsOfInterest:");
            print(fusionsOfInterest);
            if (altGene == '-') {
               print("Plot single");
               normalizedSampleGeneExpressionDf <- normalizedSampleExpressionDf[normalizedSampleExpressionDf$gene == geneOfInterest, ];
               plot.fusionExpressionOne( geneOfInterest, fusionsOfInterest, sample,
                                         normalizedSampleGeneExpressionDf )
               plotCount = plotCount + 1;
               fusionCount = fusionCount + length(fusionsOfInterest);
               singleFusionCount = singleFusionCount + length(fusionsOfInterest)
            }
            else {
               print("Plot pair");
               normalizedSampleGeneExpressionDf <- normalizedSampleExpressionDf[  normalizedSampleExpressionDf$gene == geneOfInterest
                                                                                  | normalizedSampleExpressionDf$gene == altGene , ];
               altFusions <- c(sampleFusionGeneAltDf[sampleFusionGeneAltDf$gene1 == altGene,"gene1Pos"],
                               sampleFusionGeneAltDf[sampleFusionGeneAltDf$gene2 == altGene,"gene2Pos"]);
               print("altFusions:");
               print(altFusions);
               plot.fusionExpressionPair( geneOfInterest, altGene,
                                          fusionsOfInterest, altFusions, sample,
                                          normalizedSampleGeneExpressionDf )
               plotCount = plotCount + 1;
               fusionCount = fusionCount + length(fusionsOfInterest);
            }
         }
      }
   }
   back = c(plotCount,fusionCount,singleFusionCount);
   names(back) <- c("plots", "fusionPairs", "fusionSingles");
   return(back);
}

