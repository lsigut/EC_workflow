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
# You can find example data set at https://doi.org/10.5281/zenodo.6631498
#
# Code developed by Ladislav Sigut (sigut.l@czechglobe.cz).

### Set working directory to the folder where this document is saved ===========

# This expects you are working in RStudio and this document is saved in the root
# of already existing folder structure
# - see structure_eddy()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Install and load required packages and functions ===========================

# Load eddy covariance workflow utility functions
utilities_file <- list.files(pattern = "utilities", full.names = TRUE)
source(utilities_file)

# Attach packages from GitHub
# - you might need to have RTools for Windows machine to install openeddy:
#   https://cran.r-project.org/bin/windows/Rtools
# - uses attach_pkg() function saved in utilities.R
attach_pkg("openeddy", github = "lsigut/openeddy")

# Attach packages from CRAN
# - see https://github.com/bgctw/REddyProc if mlegp not available
packages <- c("REddyProc", "bigleaf", "mlegp", "tibble")
invisible(lapply(packages, attach_pkg))

# Check if openeddy version conforms to requirements
if (packageVersion("openeddy") < package_version("0.0.0.9009"))
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9009'")

# Check if REddyProc version conforms to requirements
if (packageVersion("REddyProc") < package_version("1.3.0"))
  warning("this version of workflow works reliably with REddyProc ",
          "version '1.3.0'")

# REddyProc documentation:
# https://github.com/bgctw/REddyProc
# https://bgc.iwww.mpg.de/5624551/REddyProc-Rpackage

# Data formatting description:
# https://bgc.iwww.mpg.de/5624884/Data-Formats

# Article:
# https://www.biogeosciences.net/15/5015/2018/

### Provide metadata and set file paths and arguments ==========================

# Load the site-year settings file
settings_file <- list.files(pattern = "settings", full.names = TRUE)
source(settings_file)

# Load the list of folder structure paths
# - automated, no input required if proposed folder structure is followed
paths <- make_paths()

# Meteo variables that will be plotted, gap-filled and exported
# - FP expects Meteo columns produced during GF (even if no gaps in Meteo)
# - minimal set must be: c('Rg', 'Tair', 'Tsoil', 'VPD')
meteo <- c('Rg', 'Tair', 'Tsoil', 'VPD')

# The path where input for gap-filling is located (automated)
input <- list.files(paths$input_for_gf, pattern = paste0(siteyear, ".*txt"),
                    full.names = TRUE)[1]

# The path where file with essential variables is located (automated)
ess_in <- list.files(paths$input_for_gf, 
                     pattern = paste0(siteyear, ".*essentials.*csv"),
                     full.names = TRUE)[1]

# Timestamp of the computation
# - automated, will be included in file names
Tstamp <- format(Sys.time(), "%Y-%m-%d") 

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

### Apply uStar-filtering ======================================================

# Seasons contained within one year (e.g. Dec 2014 is pooled with Jan, Feb 2014)
# - skip estimation if fixed ustar threshold was provided
if (is.na(fixed_UT)) {
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
              paths$ustar_filtering,
              paste0("Ustar_thresholds_", Tstamp, ".csv")))
  
  # Plot saturation of NEE with UStar for available seasons
  for (i in seq_along(levels(uStarRes$season))) {
    EddyProc.C$sPlotNEEVersusUStarForSeason(
      levels(uStarRes$season)[i], dir = paths$ustar_filtering)
  }
  
  # Use annual or seasonal estimates
  UstarThres.df <- if (seasonal_ustar) {
    usGetSeasonalSeasonUStarMap(uStarRes)
  } else {
    usGetAnnualSeasonUStarMap(uStarRes)
  }
  
  # Save results to a file for fast reload if needed
  # - ?readRDS
  saveRDS(UstarThres.df, file.path(paths$ustar_filtering, "UstarThres.df.rds"))
}

### Run gap-filling ============================================================

# Fill gaps in energy fluxes with MDS gap filling algorithm 
EddyProc.C$sMDSGapFill(c('H'), FillAll = TRUE)	
EddyProc.C$sMDSGapFill(c('LE'), FillAll = TRUE)

# NEE gap filling
# - the maximum of all available seasons is taken to mark periods with low uStar
#   if seasonal_ustar == FALSE (higher exclusion fraction)
# - Note: the ustar filtering is implemented here only for nighttime and if 
#   uStar value is missing, the half hour is not filtered; thus respective NEE 
#   values are removed in quality control step (see ?set_OT_input)
if (is.na(fixed_UT)) {
  (uStarThAgg <- EddyProc.C$sGetEstimatedUstarThresholdDistribution())
  EddyProc.C$sSetUstarScenarios(usGetSeasonalSeasonUStarMap(uStarThAgg))
  EddyProc.C$sGetUstarScenarios() # check the applied thresholds
  EddyProc.C$sMDSGapFillUStarScens('NEE', FillAll = TRUE)
} else {
  # the alternative if using fixed_UT
  EddyProc.C$sMDSGapFillAfterUstar('NEE', uStarTh = fixed_UT, FillAll = TRUE)
}

# Fill gaps in variables with MDS gap filling algorithm without prior ustar 
# filtering for comparison
EddyProc.C$sMDSGapFill('NEE', FillAll = TRUE, suffix = 'UNone')

# Meteo must be gap-filled even when without gaps to run the partitioning 
for (met_var in meteo) {
  EddyProc.C$sMDSGapFill(met_var, FillAll = TRUE)
}
saveRDS(EddyProc.C, file.path(paths$gap_filling, "EddyProc.C_GF.rds"))

### Run flux partitioning ======================================================

# Perform flux partitioning for gap-filled product
suffixes <- if (is.na(fixed_UT)) {
  c('uStar', 'U05', 'U50', 'U95', 'UNone')
} else {
  c('uStar', 'UNone')
}

# Reichstein et al. (2005) - further abbreviated as MR05
# - TempRange is automatically reduced in case of not sufficient amount of data
# - if solution was not found for TempRange > 0, consider removing given suffix
for (i in suffixes) {
  TempRange <- 5
  repeat {
    suppressWarnings(
      EddyProc.C$sMRFluxPartition(
        TempVar = paste0(FP_temp, "_f"), 
        QFTempVar = paste0(FP_temp, "_fqc"),
        suffix = i,
        parsE0Regression = list(TempRange = TempRange))
    )
    TempRange <- TempRange - 1
    if (TempRange <= 0) stop("MR05 flux partitioning failed for suffix ", i)
    if (is.null(EddyProc.C$sExportResults()$E_0_NEW[1])) break
  }
}
# Save reference class and the column names at the end of MR05 FP
saveRDS(EddyProc.C, file.path(paths$gap_filling, "EddyProc.C_MR05.rds"))
MR05_cols <- colnames(EddyProc.C$sExportResults())

# Lasslop et al. (2010) - further abbreviated as GL10
FP_GL10_out_list <- vector("list", length(suffixes))
for (i in suffixes) {
  rm(EddyProc.C)
  EddyProc.C <- readRDS(file.path(paths$gap_filling, "EddyProc.C_MR05.rds"))
  EddyProc.C$sGLFluxPartition(
    TempVar = paste0(FP_temp, "_f"), 
    QFTempVar = paste0(FP_temp, "_fqc"),
    suffix = i)
  if (i == 'uStar') 
    saveRDS(EddyProc.C, 
            file.path(paths$gap_filling, "EddyProc.C_GL10_uStar.rds"))
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
saveRDS(FP_GL10_out, file.path(paths$gap_filling, "FP_GL10_out.rds"))

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
FP_MR05 <- readRDS(file.path(paths$gap_filling, "EddyProc.C_MR05.rds"))
colnames$FP_MR05 <- colnames(FP_MR05$sExportResults())

# Save the column names at the end of Lasslop et al. (2005) FP
# - contains results of only one scenario (uStar)
FP_GL10_out <- readRDS(file.path(paths$gap_filling, "FP_GL10_out.rds"))
colnames$FP_GL10 <- colnames(FP_GL10_out)

# Reload the reference class in the state after uStar FP scenario
EddyProc.C <- readRDS(file.path(paths$gap_filling, "EddyProc.C_GL10_uStar.rds"))

### Plot the results using REddyProc ===========================================

# Plot daily sums of fluxes with their uncertainties 
for (Var in c("H_f", "LE_f", "NEE_uStar_f")) {
  EddyProc.C$sPlotDailySums(Var, paste0(Var, "sd"), 
                            Dir = paths$plots, Format = plot_as)
}
# Plot fingerprints of relevant variables (with gaps and also gap-filled)
FP_vars <- c(paste0(meteo, "_f"), "H", "LE", "NEE", "H_f", "LE_f", 
             "NEE_uStar_f", "Reco_uStar", "GPP_uStar_f", "Reco_DT_uStar", 
             "GPP_DT_uStar")
for (Var in FP_vars) {
  EddyProc.C$sPlotFingerprint(Var, Dir = paths$plots, Format = plot_as)
}
# Plot diurnal cycle of relevant variables (only gap-filled)
DC_vars <- c(paste0(meteo, "_f"), "H_f", "LE_f", "NEE_uStar_f", 
             "Reco_uStar", "GPP_uStar_f", "Reco_DT_uStar", "GPP_DT_uStar") 
for (Var in DC_vars) {
  EddyProc.C$sPlotDiurnalCycle(Var, Dir = paths$plots, Format = plot_as)
}

### Convert LE to ET and combine ustar filter (UF) with qc_forGF_NEE ===========

# Load essential QC file
# - remove columns with "_orig" suffixes to prevent duplication
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
             paths$gap_filling,
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
ess_vars <- choose_avail(ess_vars, names(all_out))
ess_out <- cbind(ess, all_out[ess_vars])

# Save the essential outputs
# - rounding of numerical columns is not needed as 'all_out' is already rounded
write_eddy(ess_out, 
           file.path(
             paths$gap_filling,
             paste0(siteyear, "_GF_essentials_", Tstamp, ".csv")))

# Save documentation about executed gap-filling and flux partitioning 
document_GF(all_out, Tstamp, name, mail, siteyear, lat, long, tzh,
            FP_temp, fixed_UT, seasonal_ustar, use_CPD, paths$gap_filling)

### Plot the results using openeddy ============================================

# Correct timestamp for proper display
ess_out$timestamp <- strptime_eddy(ess_out$timestamp, "%Y-%m-%d %H:%M",
                                   shift.by = shift.by)

# Saving plots with gap-filled fluxes 
# -two kinds of plots for each flux: 
# 1) _f: filled flux with original measurements incorporated
# 2) _fall: filled flux with original measurements excluded
pdf(file.path(
  paths$gap_filling,
  paste0(siteyear, "_H_f_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(ess_out, "H", "qc_H_forGF", "qc_H_forGF",
          flux_gf = "H_f")
dev.off()
pdf(file.path(
  paths$gap_filling,
  paste0(siteyear, "_H_fall_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(ess_out, "H", "qc_H_forGF", "qc_H_forGF", flux_gf = "H_fall")
dev.off()

pdf(file.path(
  paths$gap_filling,
  paste0(siteyear, "_LE_f_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(ess_out, "LE", "qc_LE_forGF", "qc_LE_forGF",
          flux_gf = "LE_f")
dev.off()
pdf(file.path(
  paths$gap_filling,
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
  paths$gap_filling,
  paste0(siteyear, "_NEE_uStar_f_MR05_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(MR05, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_f", NEE_sep = TRUE)
dev.off()
pdf(file.path(
  paths$gap_filling,
  paste0(siteyear, "_NEE_uStar_fall_MR05_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(MR05, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_fall", NEE_sep = TRUE)
dev.off()

GL10 <- ess_out
GL10_FP_names_filter <- names(GL10) %in% c("Reco_DT_uStar", "GPP_DT_uStar")
names(GL10)[GL10_FP_names_filter] <- c("Reco", "GPP")
pdf(file.path(
  paths$gap_filling,
  paste0(siteyear, "_NEE_uStar_f_GL10_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(GL10, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_f", NEE_sep = TRUE)
dev.off()
pdf(file.path(
  paths$gap_filling,
  paste0(siteyear, "_NEE_uStar_fall_GL10_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
plot_eddy(GL10, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_fall", NEE_sep = TRUE)
dev.off()

# EOF
