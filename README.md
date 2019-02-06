# DGETools
Differential Gene Expression Analysis Tools

Currenty this includes only the Volcano function for rendering Volcano plots, and writing out dge lists from limma based
model fits, including contrast based models.

Here is a more complete example of the code leading up to a volcano plot.

    d <- DGEList(counts=filteredCounts, genes=filteredGenes)
    d <- calcNormFactors(d)
    keepRows <- rowSums(round(cpm(d$counts)) >= 1) >= cut.filter*ncol(filteredCounts)
    table(keepRows)
    subsetMono <- grepl("DP|SP", filteredDesign$cellTypeShort)
    testDGE <- d[keepRows,subsetMono]
    testDGE <- calcNormFactors(testDGE)
    testDesign <- filteredDesign[subsetMono,]

    combo <- factor(paste(testDesign$cellTypeShort, testDesign$timepoint, sep="."))

    dm <- model.matrix(~0+combo)
    colnames(dm) <- levels(combo)

    cm <- makeContrasts(
      SPTwovsOne   = SP.day7 - SP.day3,
      SPThreevsOne = SP.day14 - SP.day3,
      SPThreevsTwo  = SP.day14 - SP.day7,

      DPTwovsOne   = DP.day7 - DP.day3,
      DPThreevsOne = DP.day14 - DP.day3,
      DPThreevsTwo  = DP.day14 - DP.day7,

      DPandSPTwovsOne   = (DP.day7 + SP.day7) - (DP.day3 + SP.day3),
      DPandSPThreevsOne = (DP.day14 + SP.day14) - (DP.day3 + SP.day3),
      DPandSPThreevsTwo = (DP.day14 + SP.day14) - (DP.day7 + SP.day7),
      levels=dm
    )

    vv <- voomWithQualityWeights(testDGE, design=dm)
    corfit <- duplicateCorrelation(vv, dm, block=testDesign$donorId)
    #corfit$consensus.correlation # [1] 0.1501686
    vv <- voomWithQualityWeights(testDGE, design=dm, block=testDesign$donorId, correlation=corfit$consensus.correlation)
    fit <- lmFit(vv, dm, block=testDesign$donorId, correlation=corfit$consensus.correlation)
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
  
    Volcano(fit=fit2, target="DPandSPTwovsOne",
            title="Human Airway Cells in BPD:\nDP and SP from First to Second Timepoints",
            labelList=c("IL1A", "IL1B", "IL1RN"),
            pointSize=2, labelSize=6)

