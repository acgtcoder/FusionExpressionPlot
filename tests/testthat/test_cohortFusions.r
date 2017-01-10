context("cohortFusions.R")

describe( "do.allPlots", {

#    oneGeneFusionAList <- list( geneName1= TP53, fusions1= 7590775,
#                                     sample= the_sample-1, cohortExpressionDF= fakeExpression
#    )
#    mock_plotFusionExpressionPair <- function(
#       geneName1=NA, geneName2=NA, fusions1=NA, fusions2=NA,
#       sample=NA, cohortExpressionDF=NA ) {
#       return(
#          list(
#             geneName1= geneName1, geneName2= geneName2,
#             fusions1= fusions1,   fusions2= fusions2,
#             sample= sample, cohortExpressionDF= cohortExpressionDF
#          )
#       )
#    }
#    mock_plotFusionExpressionOne <- function(
#       geneName1=NA, fusions1=NA,
#       sample=NA, cohortExpressionDF=NA ) {
#       return(
#          list(
#             geneName1= geneName1, fusions1= fusions1,
#             sample= sample, cohortExpressionDF= cohortExpressionDF
#          )
#       )
#    }
#
#    describe( "Plotting one sample", {
#       describe( "Plotting one fusion", {
#          describe("Plots acceptor side only gene fusion", {
#             wantCalledWith <- list(
#                geneName1= oneGeneFusionA$gene1[1], fusions1= oneGeneFusionA$gene1pos[1],
#                sample= oneGeneFusionA$sample[1], cohortExpressionDF= fakeExpression
#             )
#             with_mock(
#                `FusionExpressionPlot::plotFusionExpressionOne` = mock_plotFusionExpressionOne,
#                `FusionExpressionPlot::plotFusionExpressionPair` = mock_plotFusionExpressionPair,
#                expect_equal( do.allPlots( oneGeneFusionA, fakeExpression ))
#             )
#          })
#          describe("Plots donor side only gene fusion", {
#
#          })
#          describe( "Plots 2 gene fusion", {
#
#          })
#
#       })
#       describe( "Plotting two fusions", {
#          describe("Different plots (4 genes)", {
#
#          })
#          describe("Different plots, 1 shared gene (3 genes)", {
#
#          })
#          describe( "Same plot (2 shared genes), 4 different enpoints", {
#
#          })
#          describe( "Same plot (2 shared genes), 3 different enpoints", {
#
#          })
#       })
#    })
#    describe( "Multiple samples, multiple fusions, and multiple genes", {
#
#    })
})