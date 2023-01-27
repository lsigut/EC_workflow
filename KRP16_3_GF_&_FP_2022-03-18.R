### Description ================================================================

# Eddy covariance workflow part 3/4 (https://github.com/lsigut/EC_workflow) This
# code primarily aims at u*-filtering, gap-filling (GF), flux partitioning (FP)
# and preparation of input data for summaries. Original and final data are
# plotted and further statistics are computed. Computation of bootstrap u*
# thresholds and estimation of standard deviation based on look-up tables
# provides further information about measurement uncertainty. This workflow part
# relies heavily on REddyProc package. For documentation of output variable
# names please visit MPI Online Tool website:
# https://bgc.iwww.mpg.de/5622399/REddyProc
#
# Code developed by Ladislav Sigut (sigut.l@czechglobe.cz).

### Set working directory to the folder where this document is saved ===========

# This expects you are working in RStudio and this document is saved in the root
# of already existing folder structure
# - see structure_eddy()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Install and load required packages and functions ===========================

# Load eddy covariance workflow utility functions
source("utilities.R")

# Attach packages from GitHub
# - you might need to have RTools for Windows machine to install openeddy:
#   https://cran.r-project.org/bin/windows/Rtools
# - uses attach_pkg() function saved in utilities.R
attach_pkg("openeddy", github = "lsigut/openeddy")

# Attach packages from CRAN
# - see https://github.com/bgctw/REddyProc if mlegp not available
packages <- c("REddyProc", "bigleaf", "mlegp")
invisible(lapply(packages, attach_pkg))

# Workflow is currently aligned only with specific package version
# - package version should be neither higher or lower
if (packageVersion("openeddy") < "0.0.0.9006")
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9006' and above")

if (packageVersion("REddyProc") < "1.3.0")
  warning("this version of workflow works reliably with REddyProc ",
          "version '1.3.0' or higher")

# REddyProc documentation:
# https://github.com/bgctw/REddyProc
# https://bgc.iwww.mpg.de/5624551/REddyProc-Rpackage

# Data formatting description:
# https://bgc.iwww.mpg.de/5624884/Data-Formats

# Article:
# https://www.biogeosciences.net/15/5015/2018/

### Provide metadata and set file paths and arguments ==========================

# Contact information
name <- "Ladislav Sigut" # person that performed processing
mail <- "sigut.l@czechglobe.cz" # mail of the person that performed processing

# Edit the siteyear
# - included in file names
siteyear <- "KRP16"

# Edit the year
# - used when plotting to console (single year allowed)
year <- 2016

# Specify site metadata
lat <- 49.5732575 # edit site latitude
long <- 15.0787731 # edit site longtitude
tzh <- 1 # timezone hour

# Load the list of folder structure paths
# - automated, no input required if proposed folder structure is followed
paths <- structure_eddy()

# Specify meteo variables that will be plotted, gap-filled and exported
meteo <- c('Rg', 'Tair', 'Tsoil', 'VPD')

# Temperature used for flux partitioning ('Tsoil' or 'Tair') 
# - default: FP_temp <- 'Tsoil'
# - note that MDS gap-filling is based on 'Tair'
FP_temp <- 'Tsoil'

# Show precheck plots in the console? 
plot_to_console <- FALSE

# Save figures as "png" (default) or "pdf"
# - NEEvsUStar plots are fixed to "pdf"
plot_as <- "png" 

# Specify the time shift (in seconds) to be applied to the date-time information
# in order to represent the center of averaging period
shift.by <- -900

# The path where input for gap-filling is located (automated)
input <- grep(paste0(siteyear, ".*txt"), 
              list.files(paths$Input_for_GF, full.names = TRUE),
              value = TRUE)

# The path where file with essential variables is located (automated)
ess_in <- grep(paste0(siteyear, ".*essentials.*csv"), 
               list.files(paths$Input_for_GF, full.names = TRUE),
               value = TRUE)

# Timestamp of the computation
# - automated, will be included in file names
Tstamp <- format(Sys.time(), "%Y-%m-%d") 

# Choose ustar threshold estimation method resolution
# - either seasonal ustar threshold (seasonal_ustar <- TRUE; default)
# - or annual thresholds (seasonal_ustar <- FALSE)
# - seasonal UT is recommended as it keeps more data (see Wutzler et al., 2018)
seasonal_ustar <- TRUE

# Should change-point detection (CPD) method (Barr et al. 2013) be used 
# (use_CPD <- TRUE) instead of classical approach (Reichstein et al. 2005, 
# Papale et al. 2006) using binned classes of ustar and temperature? 
# (use_CPD <- FALSE; default)
# The changepoint is estimated based on the entire subset within one season 
# and one temperature class; currently using argument 'isUsingCPTSeveralT' 
# of function usControlUstarEst()
use_CPD <- FALSE

### Prepare REddyProc input data ===============================================

# Load data with one header and one unit row from (tab-delimited) text file
EddyData.F <- read_eddy(input, sep = "\t")
# head(EddyData.F)
# str(EddyData.F)

# If not provided or if including gaps, calculate VPD from Tair and rH
if (!"VPD" %in% names(EddyData.F) || anyNA(EddyData.F$VPD)) {
  EddyData.F$VPD <- fCalcVPDfromRHandTair(EddyData.F$rH, EddyData.F$Tair)
}

# Add time stamp in POSIX time format
EddyDataWithPosix.F <- fConvertTimeToPosix(
  EddyData.F, 'YDH', Year = 'Year', Day = 'DoY', Hour = 'Hour')

# Initalize R5 reference class sEddyProc for processing of eddy data
# with all variables needed for processing later
variables <- c('NEE', 'LE', 'H', 'Rg','Tair', 'Tsoil', 'VPD', 'Ustar')
EddyProc.C <- sEddyProc$new(siteyear, EddyDataWithPosix.F, variables)
EddyProc.C$sSetLocationInfo(lat, long, tzh)  # site location info

# See the content
str(EddyProc.C)
EddyProc.C$sPrintFrames(NumRows.i = 6L)

# Plot individual months/years to console (of current R graphics device)
if (plot_to_console) {
  for (Var in variables) {
    EddyProc.C$sPlotFingerprintY(Var, Year = year)
    title(main = Var, line = 0.2)
  }
  for (Var in variables) {
    EddyProc.C$sPlotHHFluxesY(Var, Year = year)
    title(ylab = Var, line = 2)
  }
}

### Apply uStar-filtering ======================================================

# Seasons contained within one year (e.g. Dec 2014 is pooled with Jan, Feb 2014)
set.seed(0815)
season_factor <- usCreateSeasonFactorMonthWithinYear(
  EddyDataWithPosix.F$DateTime + shift.by)
table(season_factor)
(uStarRes <- EddyProc.C$sEstUstarThresholdDistribution(
  nSample = 200L, seasonFactor = season_factor,
  ctrlUstarEst = usControlUstarEst(isUsingCPTSeveralT = use_CPD)))

# Round and save the results
uStarRes <- round_df(uStarRes)
write.csv(uStarRes, row.names = FALSE, 
          file.path(
            paths$Ustar_filtering,
            paste0("Ustar_thresholds_", Tstamp, ".csv")))

# Plot saturation of NEE with UStar for available seasons
for (i in seq_along(levels(uStarRes$season))) {
  EddyProc.C$sPlotNEEVersusUStarForSeason(
    levels(uStarRes$season)[i], dir = paths$Ustar_filtering)
}

# Use annual or seasonal estimates
UstarThres.df <- if (seasonal_ustar) {
  usGetSeasonalSeasonUStarMap(uStarRes)
} else {
  usGetAnnualSeasonUStarMap(uStarRes)
}

# Save results to a file for fast reload if needed
# - ?readRDS
saveRDS(UstarThres.df, file.path(paths$Ustar_filtering, "UstarThres.df.rds"))

### Run gap-filling ============================================================

# Fill gaps in energy fluxes with MDS gap filling algorithm 
EddyProc.C$sMDSGapFill(c('H'), FillAll = TRUE)	
EddyProc.C$sMDSGapFill(c('LE'), FillAll = TRUE)

# NEE gapfilling: the maximum of all available seasons is taken to mark periods 
# with low uStar if seasonal_ustar == FALSE (higher exclusion fraction)
# NB: the ustar filtering is implemented here only for nighttime and if ustar
# value is missing the halfhour is NOT filtered (last checked in 2017)
(uStarThAgg <- EddyProc.C$sGetEstimatedUstarThresholdDistribution())
EddyProc.C$sSetUstarScenarios(usGetSeasonalSeasonUStarMap(uStarThAgg))
EddyProc.C$sGetUstarScenarios() # check the applied thresholds
EddyProc.C$sMDSGapFillUStarScens('NEE', FillAll = TRUE)

# Fill gaps in variables with MDS gap filling algorithm without prior ustar 
# filtering for comparison
EddyProc.C$sMDSGapFill('NEE', FillAll = TRUE, suffix = 'UNone')

# Meteo must be gap-filled even when without gaps to run the partitioning 
for (met_var in meteo) {
  EddyProc.C$sMDSGapFill(met_var, FillAll = TRUE)
}
saveRDS(EddyProc.C, file.path(paths$Gap_filling, "EddyProc.C_GF.rds"))

### Run flux partitioning ======================================================

# Perform flux partitioning for gap-filled product
suffixes <- c('uStar', 'U05', 'U50', 'U95', 'UNone')

# Reichstein et al. (2005) - further abbreviated as MR05
# - in case not sufficient amount of data, try to reduce TempRange 
for (i in suffixes) {
  EddyProc.C$sMRFluxPartition(
    TempVar = paste(FP_temp, "_f", sep = ""), 
    QFTempVar = paste(FP_temp, "_fqc", sep = ""),
    suffix = i,
    parsE0Regression = list(TempRange = 5))
}
# Save reference class and the column names at the end of MR05 FP
saveRDS(EddyProc.C, file.path(paths$Gap_filling, "EddyProc.C_MR05.rds"))
MR05_cols <- colnames(EddyProc.C$sExportResults())

# Lasslop et al. (2010) - further abbreviated as GL10
FP_GL10_out_list <- vector("list", length(suffixes))
for (i in suffixes) {
  rm(EddyProc.C)
  EddyProc.C <- readRDS(file.path(paths$Gap_filling, "EddyProc.C_MR05.rds"))
  EddyProc.C$sGLFluxPartition(
    TempVar = paste(FP_temp, "_f", sep = ""), 
    QFTempVar = paste(FP_temp, "_fqc", sep = ""),
    suffix = i)
  if (i == 'uStar') 
    saveRDS(EddyProc.C, 
            file.path(paths$Gap_filling, "EddyProc.C_GL10_uStar.rds"))
  out <- EddyProc.C$sExportResults()
  cols <- colnames(out)
  cols_out <- cols[!cols %in% MR05_cols]
  out <- out[cols_out]
  corr_filter <- !cols_out %in% grep(paste(suffixes, collapse = "|"), 
                                     cols_out, value = TRUE)
  names(out)[corr_filter] <- paste(cols_out[corr_filter], i, sep = '_')
  FP_GL10_out_list[[which(suffixes == i)]] <- out
}

# Column-bind the list with results to single data frame and save GL10 results
FP_GL10_out <- do.call(cbind, FP_GL10_out_list)
saveRDS(FP_GL10_out, file.path(paths$Gap_filling, "FP_GL10_out.rds"))

### Make following sections independent on previous processing =================

# Reload objects created during previous steps
# - uStar-filtering, gap-filling (GF) and flux partitioning (FP) are 
#   computationally expensive operations
# - sEddyProc reference class object or its outputs are saved during different 
#   processing stages
# - if above sections were run at least once for given setup, uStar-filtering,
#   GF and FP can be omitted if only amends in further steps are required 
#   (note that all sections preceding 'Apply uStar-filtering' must be run) 
# - saving column names from different stages of post-processing simplifies
#   column selection for different output files and keeps output structure

# Save snapshots of result colnames in steps to a list
colnames <- list()

# Save the column names at the end of Reichstein et al. (2005) FP
FP_MR05 <- readRDS(file.path(paths$Gap_filling, "EddyProc.C_MR05.rds"))
colnames$FP_MR05 <- colnames(FP_MR05$sExportResults())

# Save the column names at the end of Lasslop et al. (2005) FP
# - contains results of only one scenario (uStar)
FP_GL10_out <- readRDS(file.path(paths$Gap_filling, "FP_GL10_out.rds"))
colnames$FP_GL10 <- colnames(FP_GL10_out)

# Reload the reference class in the state after uStar FP scenario
EddyProc.C <- readRDS(file.path(paths$Gap_filling, "EddyProc.C_GL10_uStar.rds"))

### Plot the results using REddyProc ===========================================

# Plot daily sums of fluxes with their uncertainties 
for (Var in c("H_f", "LE_f", "NEE_uStar_f")) {
  EddyProc.C$sPlotDailySums(Var, paste(Var, "sd", sep = ""), 
                            Dir = paths$Plots, Format = plot_as)
}
# Plot fingerprints of relevant variables (with gaps and also gap-filled)
FP_vars <- c(paste(meteo, "_f", sep = ""), "H", "LE", "NEE", "H_f", "LE_f", 
             "NEE_uStar_f", "Reco_uStar", "GPP_uStar_f", "Reco_DT_uStar", 
             "GPP_DT_uStar")
for (Var in FP_vars) {
  EddyProc.C$sPlotFingerprint(Var, Dir = paths$Plots, Format = plot_as)
}
# Plot diurnal cycle of relevant variables (only gap-filled)
DC_vars <- c(paste(meteo, "_f", sep = ""), "H_f", "LE_f", "NEE_uStar_f", 
             "Reco_uStar", "GPP_uStar_f", "Reco_DT_uStar", "GPP_DT_uStar") 
for (Var in DC_vars) {
  EddyProc.C$sPlotDiurnalCycle(Var, Dir = paths$Plots, Format = plot_as)
}

### Convert LE to ET and combine ustar filter (UF) with qc_forGF_NEE ===========

# Load essential QC file
# - remove columns with "_orig" suffices to prevent duplication
ess <- read_eddy(ess_in)
ess[grep("_orig$", names(ess), value = TRUE)] <- NULL

# Export input data, gap-filling and flux partitioning results
all_out <- cbind(ess["timestamp"], EddyData.F, FP_MR05$sExportResults(), 
                 FP_GL10_out)

# Convert LE to ET [mm hour-1]
# - columns other than LE_f are required for agg_fsd() uncertainty evaluation
# - LE_fqc is included in the conversion only to keep proper ordering 
LE_vars <- c("LE_orig", "LE_f", "LE_fqc", "LE_fall", "LE_fsd")
ET_vars <- gsub("LE", "ET", LE_vars)
all_out[, ET_vars] <- 
  lapply(LE_vars, 
         function(x) LE.to.ET(all_out[, x], all_out$Tair_f) * 3600)
openeddy::units(all_out[ET_vars]) <- rep("mm hour-1", length(ET_vars))

# Overwrite ET_fqc with proper values
all_out$ET_fqc <- all_out$LE_fqc
openeddy::units(all_out$ET_fqc) <- "-" 

# Data are excluded (flag 2) by ustar filter (UF) when Ustar_uStar_fqc > 0
qc_uStar <- as.data.frame(ifelse(all_out["Ustar_uStar_fqc"] > 0, 2, 0))
names(qc_uStar) <- "qc_uStar"

# Create resulting QC flag for NEE
all_out$qc_NEE_forGF_UF <- combn_QC(cbind(ess["qc_NEE_forGF"], qc_uStar), 
                                    c("qc_NEE_forGF", "qc_uStar"), 
                                    "qc_NEE_forGF_UF")

### Save outputs and documentation =============================================

# Round and save all results together with input data into CSV file
all_out <- round_df(all_out)
write_eddy(all_out,
           file.path(
             paths$Gap_filling,
             paste0(siteyear, "_GF_full_output_", Tstamp, ".csv")))

# Select the most important variables obtained during gap-filling
# - add also ET columns created by conversion from LE
gf_ess <- grep(paste(c("_orig$", "_f$", "_fqc$", "_fall$", "_fsd$", "_Thres$"), 
                     collapse = "|"), colnames(all_out), value = TRUE)
gf_ess <- gf_ess[!grepl("U05|U50|U95|UNone|GPP", gf_ess)]

# Select the most important variables obtained during flux partitioning
fp_ess <- grep("DT_uStar", colnames$FP_GL10, value = TRUE)

ess_vars <- c("qc_NEE_forGF_UF", gf_ess, "Reco_uStar", "GPP_uStar_f", fp_ess)

# Create object with essential variables including gap-filling (GF) 
# and flux partitioning results
ess_vars <- choose_avail(names(all_out), ess_vars)
ess_out <- cbind(ess, all_out[ess_vars])

# Save the essential outputs
# - rounding of numerical columns is not needed as 'all_out' is already rounded
write_eddy(ess_out, 
           file.path(
             paths$Gap_filling,
             paste0(siteyear, "_GF_essentials_", Tstamp, ".csv")))

# Save documentation about executed gap-filling and flux partitioning 
document_GF(all_out, Tstamp, name, mail, siteyear, lat, long, tzh,
            FP_temp, seasonal_ustar, use_CPD, paths$Gap_filling)

### Plot the results using openeddy ============================================

# Correct timestamp for proper display
ess_out$timestamp <- strptime_eddy(ess_out$timestamp, "%Y-%m-%d %H:%M",
                                   shift.by = shift.by)

# Saving plots with gap-filled fluxes 
# -two kinds of plots for each flux: 
# 1) _f: filled flux with original measurements incorporated
# 2) _fall: filled flux with original measurements excluded
pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_H_f_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(ess_out, "H", "qc_H_forGF", "qc_H_forGF",
          flux_gf = "H_f")
dev.off()
pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_H_fall_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(ess_out, "H", "qc_H_forGF", "qc_H_forGF", flux_gf = "H_fall")
dev.off()

pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_LE_f_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(ess_out, "LE", "qc_LE_forGF", "qc_LE_forGF",
          flux_gf = "LE_f")
dev.off()
pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_LE_fall_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(ess_out, "LE", "qc_LE_forGF", "qc_LE_forGF",
          flux_gf = "LE_fall")
dev.off()

# Saving plots with gap-filled and flux partitioned data 
# -two kinds of flux partitioning for NEE:
# 1) _MR05: Reichstein et al. (2005)
# 2) _GL10: Lasslop et al. (2010)
# -two kinds of plots for each flux: 
# 1) _f: filled flux with original measurements incorporated
# 2) _fall: filled flux with original measurements excluded
MR05 <- ess_out
MR05_FP_names_filter <- names(MR05) %in% c("Reco_uStar", "GPP_uStar_f")
names(MR05)[MR05_FP_names_filter] <- c("Reco", "GPP")
pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_NEE_uStar_f_MR05_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(MR05, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_f", NEE_sep = TRUE)
dev.off()
pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_NEE_uStar_fall_MR05_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(MR05, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_fall", NEE_sep = TRUE)
dev.off()

GL10 <- ess_out
GL10_FP_names_filter <- names(GL10) %in% c("Reco_DT_uStar", "GPP_DT_uStar")
names(GL10)[GL10_FP_names_filter] <- c("Reco", "GPP")
pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_NEE_uStar_f_GL10_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(GL10, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_f", NEE_sep = TRUE)
dev.off()
pdf(file.path(
  paths$Gap_filling,
  paste0(siteyear, "_NEE_uStar_fall_GL10_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(GL10, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_fall", NEE_sep = TRUE)
dev.off()

# EOF
