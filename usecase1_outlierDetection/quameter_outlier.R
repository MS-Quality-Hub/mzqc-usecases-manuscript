# Step 1: Convert QuaMeter output to mzQC
#system("python quameter_to_mzqc.py --input_file Mtb-120-outlier-metrics.tsv --location_prefix ftp://massive.ucsd.edu/v01/MSV000079012/")

# Step 2: Load mzQC and extract relevant metrics
if (!requireNamespace("rmzqc", quietly = TRUE) ||
    packageVersion("rmzqc") < "0.5.6") {
  stop("rmzqc (>= 0.5.6) is required. Please install or update it.")
}
library(rmzqc)
Mtb <- readMZQC("Mtb-120-outlier-metrics.mzQC")


# Simplified helper to get the value of a named metric
extract_metric_by_name <- function(run, metric_name) {
  run$getMetric(name = metric_name)[[1]]$value
}


# Extract metrics and convert to data.frame
metrics <- lapply(Mtb$runQualities, function(run) {
  list(
    Filename = run$metadata$inputFiles[[1]]$name,
    XIC50 = extract_metric_by_name(run, "XIC50 fraction"),
    XICFWHMQ2 = extract_metric_by_name(run, "XIC-FWHM quantiles")[[2]],
    XICHeightQ2 = extract_metric_by_name(run, "XIC-Height quartile ratios")[[2]],
    TICQuartRTMiddle = sum(unlist(extract_metric_by_name(run, "TIC quarters RT fraction")[2:3])),
    MS1QuartRTMiddle = sum(unlist(extract_metric_by_name(run, "MS1 quarter RT fraction")[2:3])),
    MS2QuartRTMiddle = sum(unlist(extract_metric_by_name(run, "MS2 quarter RT fraction")[2:3])),
    MS1TICDeltaQuartRatioQ2 = extract_metric_by_name(run, "MS1 TIC-change quartile ratios")[[2]],
    MS1TICQuartRatioQ2 = extract_metric_by_name(run, "MS1 TIC quartile ratios")[[2]],
    MS1Count = extract_metric_by_name(run, "number of MS1 spectra"),
    MS2Count = extract_metric_by_name(run, "number of MS2 spectra"),
    MS1DensityQ2 = extract_metric_by_name(run, "MS1 density quantiles")[[2]],
    MS2DensityQ2 = extract_metric_by_name(run, "MS2 density quantiles")[[2]],
    MS2PreZ2 = extract_metric_by_name(run, "MS2 known precursor charges fractions")[["UO:0000191"]][[2]],
    MS2PreZ3 = extract_metric_by_name(run, "MS2 known precursor charges fractions")[["UO:0000191"]][[3]],
    MS1FastRate = extract_metric_by_name(run, "fastest frequency for MS level 1 collection"),
    MS2FastRate = extract_metric_by_name(run, "fastest frequency for MS level 2 collection")
  )
})



# Combine the list of named lists into a data frame
MtbDF <- do.call(rbind.data.frame, metrics)

rownames(MtbDF) <- MtbDF$Filename
MtbDF$Filename <- NULL


# Step 3: Robust PCA
library(MASS)
library(ggplot2)
library(ggfortify)

set.seed(123)
#set.seed(456)
#set.seed(789)

robust.cov <- cov.rob(MtbDF)  ## uses random sampling, not reproducible without fixed seed!
robust.cor <- cov2cor(robust.cov$cov)

robust.cor.1 <- robust.cov
robust.cor.1$cov <- robust.cor

metrics.pca.1 <- sweep(MtbDF, 2, sqrt(diag(robust.cov$cov)), FUN = "/")
robust.cor.1$center <- robust.cov$center / sqrt(diag(robust.cov$cov))

pc.cr <- princomp(metrics.pca.1, covmat = robust.cor.1, scores = TRUE)
summary(pc.cr)

# Step 4: Compute distances and medians
distMatrix <- as.matrix(dist(pc.cr$scores[,1:5], method = "euclidean"))
Medians <- apply(distMatrix, 1, median, na.rm = TRUE)
names(Medians) <- rownames(distMatrix)



# Step 5: Write output mzQC
library(rmzqc)
inputFile <- rmzqc::MzQCinputFile(name = "Mtb-120-outlier-metrics.mzQC", 
                                   location = "https://github.com/HUPO-PSI/mzQC/blob/main/specification_documents/examples/Mtb-120-outlier-metrics.mzQC",
                                   fileFormat = getCVTemplate("MS:1003160"))

Structure <- getQualityMetricTemplate("MS:4000005", getCVSingleton()) # "table"
Structure$value <- list(Filename = names(Medians), Distance = Medians)

set_qc <- MzQCsetQuality$new(
  metadata = MzQCmetadata$new(
    label = "Outlier Analysis by robust PCA",
    inputFiles = list(inputFile),
    analysisSoftware = list(toAnalysisSoftware(id = "MS:1000799", version = "20250124")) # custom unreleased software tool)
  ),
  qualityMetrics = list(Structure)
)

mzQC_document <- MzQCmzQC$new(
  version = "1.0.0",
  creationDate = MzQCDateTime$new(),
  contactName = "David Tabb",
  contactAddress = "d.l.tabb@rug.nl",
  description = "An example of setQuality used for outliers",
  runQualities = list(),
  setQualities = list(set_qc),
  controlledVocabularies = list(getCVInfo())
)

writeMZQC("Mtb-outlier-PCA.mzQC", mzQC_document)
