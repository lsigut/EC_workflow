### Description ================================================================

# This code primarily aims at u*-filtering, gap-filling (GF), 
# flux partitioning (FP) and preparation of input data for summaries. 
# Original and final data are plotted and further statistics are computed.
# Computation of bootstrap u* thresholds and estimation of standard deviation 
# based on look-up tables provides further information about measurement 
# uncertainty.
# Code developed by Ladislav Šigut (sigut.l@czechglobe.cz).  

### Loading the required packages ==============================================

# You might need to install the packages first. In that case use:
# install.packages("packagename")
library(openeddy)
library(REddyProc)
# Documentation:
# https://www.bgc-jena.mpg.de/bgi/index.php/Services/REddyProcWebRPackage
# https://github.com/bgctw/REddyProc
# Data formatting description:
# https://www.bgc-jena.mpg.de/bgi/index.php/Services/REddyProcWebFormats
# Article:
# https://www.biogeosciences.net/15/5015/2018/

### User defined variables =====================================================

# Firstly manually set the desired working directory in RStudio interface
# (typically it is the folder specifying the siteyear, e.g. BKF16)
# This assumes that the subdirectory structure is already present
# (i.e. the !Site_structure_template from the server)

# Set filenames, paths and arguments
siteyear <- "KRP16" # used as output filename and as 'ID.s' in sEddyProc$new()
Tstamp <- format(Sys.time(), "%Y-%m-%d") # Timestamp of the computation
# (will be included in file names)
path_out <- "./Level 3/Gap-filling/REddyProc/"
(path_plots <- paste(path_out, "Plots/", sep = ""))
(path_UF <- paste(path_out, "Ustar filtering/", sep = ""))
pkg_version <- packageVersion("REddyProc")
input <- "./Level 2/Input for gap-filling/KRP16_2018-08-10.txt"
lat <- 49.6 # edit site latitude
long <- 15.1 # edit site longtitude
tz <- 1 # timezone
# Include meteo variables that will be plotted, gap-filled and exported
meteo <- c('Rg', 'Tair', 'Tsoil', 'VPD')
# Temperature used for flux partitioning ('Tsoil' or 'Tair'). 
# Default: FP_temp <- 'Tsoil'
# NB: all other procedures are based on 'Tair'
FP_temp <- 'Tsoil'
# Show precheck plots in the console? Default: plot_to_console <- TRUE
plot_to_console <- FALSE
year <- 2016 # used when plotting to console (single year allowed)
# Save plots as "png" (default) or "pdf"; NEEvsUStar plots are fixed to "pdf"
plot_as <- "png" 
# Choose between seasonal ustar threshold (seasonal_ustar <- TRUE; default) 
# or annual thresholds (seasonal_ustar <- FALSE)
# Seasonal UT is recommended since it keeps more data (see manuscript).
seasonal_ustar <- TRUE
# Should change-point detection (CPD) method (Barr et al. 2013) be used 
# (use_CPD <- TRUE) instead of classical approach (Reichstein et al. 2005, 
# Papale et al. 2006) using binned classes of ustar and temperature? 
# (use_CPD <- FALSE; default)
# The changepoint is estimated based on the entire subset within one season 
# and one temperature class; currently using argument 'isUsingCPTSeveralT' 
# of function usControlUstarEst()
use_CPD <- FALSE

### Data preparation =============================================================

# Load data with one header and one unit row from (tab-delimited) text file
EddyData.F <- read_eddy(input, sep = "\t")
# head(EddyData.F)
# str(EddyData.F)

# If not provided, calculate VPD from Tair and rH
if (!"VPD" %in% names(EddyData.F)) {
  EddyData.F$VPD <- fCalcVPDfromRHandTair(EddyData.F$rH, EddyData.F$Tair)
}

# Add time stamp in POSIX time format
EddyDataWithPosix.F <- fConvertTimeToPosix(
  EddyData.F, 'YDH', Year.s = 'Year', Day.s = 'DoY', Hour.s = 'Hour')

# Initalize R5 reference class sEddyProc for processing of eddy data
# with all variables needed for processing later
variables <- c('NEE', 'LE', 'H', 'Rg','Tair', 'Tsoil', 'VPD', 'Ustar')
EddyProc.C <- sEddyProc$new(siteyear, EddyDataWithPosix.F, variables)
EddyProc.C$sSetLocationInfo(lat, long, tz)  # site location info

# See the content
str(EddyProc.C)
EddyProc.C$sPrintFrames(NumRows.i = 6L)

# Plot individual months/years to console (of current R graphics device)
if (plot_to_console) {
  for (Var in variables) {
    EddyProc.C$sPlotFingerprintY(Var, Year.i = year)
    title(main = Var, line = 0.2)
  }
  for (Var in variables) {
    EddyProc.C$sPlotHHFluxesY(Var, Year.i = year)
    title(ylab = Var, line = 2)
  }
}

### Ustar filtering ============================================================

# Seasons contained within one year (e.g. Dec 2014 is pooled with Jan, Feb 2014)
set.seed(0815)
season_factor <- usCreateSeasonFactorMonthWithinYear(
  EddyDataWithPosix.F$DateTime - 15 * 60)
table(season_factor)
(uStarRes <- EddyProc.C$sEstUstarThresholdDistribution(
  nSample = 200L, seasonFactor.v = season_factor,
  ctrlUstarEst.l = usControlUstarEst(isUsingCPTSeveralT = use_CPD)))
# Save the results
write.csv(uStarRes, row.names = FALSE, 
          paste(path_UF, "Ustar_thresholds_", Tstamp, ".csv", sep = ""))
# Plot saturation of NEE with UStar for available seasons
for (i in seq_along(levels(uStarRes$season))) {
  EddyProc.C$sPlotNEEVersusUStarForSeason(
    levels(uStarRes$season)[i], dir = path_UF)
}
# Use annual or seasonal estimates
UstarThres.df <- if (seasonal_ustar) {
  usGetSeasonalSeasonUStarMap(uStarRes)
} else {
  usGetAnnualSeasonUStarMap(uStarRes)
}
saveRDS(UstarThres.df, paste0(path_UF, "UstarThres.df.rds"))

### Gap-filling ================================================================

# Save snapshots of result colnames in steps to a list
colnames <- list()

# Fill gaps in energy fluxes with MDS gap filling algorithm 
EddyProc.C$sMDSGapFill(c('H'), FillAll.b = TRUE)
EddyProc.C$sMDSGapFill(c('LE'), FillAll.b = TRUE)
# Gapfilling: the maximum of all available seasons is taken to mark periods 
# with low uStar if seasonal_ustar == FALSE (higher exclusion fraction)
# NB: the ustar filtering is implemented here only for nighttime and if ustar
# value is missing the halfhour is NOT filtered (last checked in 2017)
EddyProc.C$sMDSGapFillAfterUStarDistr('NEE', FillAll.b = TRUE, 
                                      UstarThres.df = UstarThres.df)
# Fill gaps in variables with MDS gap filling algorithm 
# without prior ustar filtering for comparison
EddyProc.C$sMDSGapFill('NEE', FillAll.b = TRUE, Suffix.s = 'UNone')
# Meteo must be gap-filled even when without gaps to run the partitioning 
for (met_var in meteo) {
  EddyProc.C$sMDSGapFill(met_var, FillAll.b = TRUE)
}
colnames$GF <- colnames(EddyProc.C$sExportResults())
saveRDS(EddyProc.C, file = paste0(path_out, "EddyProc.C_GF.rds"))

### Flux partitioning ==========================================================

# Perform flux partitioning for gap-filled product
suffixes <- c('uStar', 'U05', 'U50', 'U95', 'UNone')

# Reichstein et al. (2005)
for (i in suffixes) {
  EddyProc.C$sMRFluxPartition(
    TempVar.s = paste(FP_temp, "_f", sep = ""), 
    QFTempVar.s = paste(FP_temp, "_fqc", sep = ""),
    Suffix.s = i)
}
colnames$FP_MR05 <- colnames(EddyProc.C$sExportResults())
saveRDS(EddyProc.C, file = paste0(path_out, "EddyProc.C_MR05.rds"))

# Lasslop et al. (2010)
FP_GL10_out <- EddyProc.C$sExportResults()[, 0]
for (i in suffixes) {
  rm(EddyProc.C)
  EddyProc.C <- readRDS(paste0(path_out, "EddyProc.C_MR05.rds"))
  EddyProc.C$sGLFluxPartition(
    TempVar.s = paste(FP_temp, "_f", sep = ""), 
    QFTempVar.s = paste(FP_temp, "_fqc", sep = ""),
    Suffix.s = i)
  if (i == 'uStar') saveRDS(EddyProc.C, 
                            paste0(path_out, "EddyProc.C_GL10_uStar.rds"))
  out <- EddyProc.C$sExportResults()
  cols <- colnames(out)
  cols_out <- cols[!cols %in% colnames$FP_MR05]
  out <- out[cols_out]
  corr_filter <- !cols_out %in% grep(paste(suffixes, collapse = "|"), 
                                     cols_out, value = TRUE)
  names(out)[corr_filter] <- paste(cols_out[corr_filter], i, sep = '_')
  FP_GL10_out <- cbind(FP_GL10_out, out)
}
rm(EddyProc.C)
EddyProc.C <- readRDS(paste0(path_out, "EddyProc.C_GL10_uStar.rds"))
colnames$FP_GL10 <- colnames(FP_GL10_out)

### Plotting the results =======================================================

# Plot daily sums of fluxes with their uncertainties 
for (Var in c("H_f", "LE_f", "NEE_uStar_f")) {
  EddyProc.C$sPlotDailySums(Var, paste(Var, "sd", sep = ""), 
                            Dir.s = path_plots, Format.s = plot_as)
}
# Plot fingerprints of relevant variables (with gaps and also gap-filled)
FP_vars <- c(paste(meteo, "_f", sep = ""), "H", "LE", "NEE", "H_f", "LE_f", 
             "NEE_uStar_f", "Reco_uStar", "GPP_uStar_f", "Reco_DT_uStar", 
             "GPP_DT_uStar")
for (Var in FP_vars) {
  EddyProc.C$sPlotFingerprint(Var, Dir.s = path_plots, Format.s = plot_as)
}
# Plot diurnal cycle of relevant variables (only gap-filled)
DC_vars <- c(paste(meteo, "_f", sep = ""), "H_f", "LE_f", "NEE_uStar_f", 
             "Reco_uStar", "GPP_uStar_f", "Reco_DT_uStar", "GPP_DT_uStar") 
for (Var in DC_vars) {
  EddyProc.C$sPlotDiurnalCycle(Var, Dir.s = path_plots, Format.s = plot_as)
}

### Combine ustar filter (UF) with qc_forGF_NEE ================================

# load essential QC file
ess_path <- paste("./Level 2/Input for gap-filling/", 
                  "KRP16_forGF_QC_essentials_2018-08-10.csv", sep = "")
ess_in <- read_eddy(ess_path)
# Export input data, gap-filling and flux partitioning results
GF_MR05 <- readRDS(paste0(path_out, "EddyProc.C_MR05.rds"))
all_out <- cbind(EddyData.F, GF_MR05$sExportResults(), FP_GL10_out)
# Data are excluded (flag 2) by ustar filter (UF) when Ustar_uStar_fqc > 0
qc_uStar <- as.data.frame(ifelse(all_out["Ustar_uStar_fqc"] > 0, 2, 0))
names(qc_uStar) <- "qc_uStar"
# Create resulting QC flag for NEE
all_out$qc_NEE_forGF_UF <- combn_QC(cbind(ess_in["qc_NEE_forGF"], qc_uStar), 
                                   c("qc_NEE_forGF", "qc_uStar"), 
                                   "qc_NEE_forGF_UF")

### Saving outputs, settings and plots =========================================

# Save all results together with input data into CSV file 
write_eddy(all_out, 
           paste(path_out, siteyear, "_GF_full_output_", Tstamp, ".csv", 
                 sep = ""))
# Create essentials file with gap-filling (GF) results
cols_MR05 <- colnames$FP_MR05[!colnames$FP_MR05 %in% colnames$GF]
ess_vars <- c(
  grep(paste(c("_f$", "_fqc$", "_fall$", "_fsd$", "_Thres$"), 
             collapse = "|"), colnames$GF, value = TRUE),
  "qc_NEE_forGF_UF",
  grep(paste(c("^E_0", "^R_ref", "^Reco", "^GPP[[:print:]]+f$"), 
             collapse = "|"), cols_MR05, value = TRUE),
  grep(paste(c("GPP2000", "_k_", "beta", "alpha", "RRef_[^N]", "E0", "Reco", 
               "GPP"), collapse = "|"), names(FP_GL10_out), value = TRUE))

ess_vars[!ess_vars %in% names(all_out)]
ess_vars <- ess_vars[ess_vars %in% names(all_out)]
ess_out <- cbind(ess_in, all_out[ess_vars])
write_eddy(ess_out, paste(path_out, siteyear, "_GF_essentials_", Tstamp, ".csv", 
                          sep = ""))

# Save settings and other info
perc_records <- nrow(all_out) / 100
flux_avail_names <- c("H_orig", "LE_orig", "NEE_uStar_orig")
avail_rec <- lapply(flux_avail_names, function(x) {
  temp <- table(!is.na(all_out[x]))
  unname(temp["TRUE"] / perc_records)
})
names(avail_rec) <- flux_avail_names
settings <- list(REddyProc_pkg_version = pkg_version,
                 siteyear = siteyear, 
                 path_out = path_out,
                 path_plots = path_plots,
                 input = input,
                 latitude = lat,
                 longitude = long, 
                 timezone = tz, 
                 meteo_variables = meteo,
                 flux_partitioning_temperature = FP_temp,
                 seasonal_ustar = seasonal_ustar,
                 use_changepoint_detection = use_CPD,
                 available_H_perc = avail_rec$H_orig,
                 available_LE_perc = avail_rec$LE_orig,
                 available_NEE_perc_after_UF = avail_rec$NEE_uStar_orig)
sink(paste(path_out, siteyear, '_settings.txt', sep = ""))
settings
sink()

# Saving plots with gap-filled and flux partitioned data 
ess_out$timestamp <- strptime_eddy(ess_out$timestamp, "%Y-%m-%d %H:%M",
                                      center.by = -900)
head(ess_out$timestamp)
names(ess_out)

# -two kinds of plots for each flux: 
# 1) _f: filled flux with original measurements incorporated
# 2) _fall: filled flux with original measurements excluded
pdf(paste(path_out, siteyear, "_H_f_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(ess_out, "H_stc", "qc_H_stc_forGF", "qc_H_stc_forGF",
          flux_gf = "H_f")
dev.off()
pdf(paste(path_out, siteyear, "_H_fall_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(ess_out, "H_stc", "qc_H_stc_forGF", "qc_H_stc_forGF",
          flux_gf = "H_fall")
dev.off()

pdf(paste(path_out, siteyear, "_LE_f_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(ess_out, "LE_stc", "qc_LE_stc_forGF", "qc_LE_stc_forGF",
          flux_gf = "LE_f")
dev.off()
pdf(paste(path_out, siteyear, "_LE_fall_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(ess_out, "LE_stc", "qc_LE_stc_forGF", "qc_LE_stc_forGF",
          flux_gf = "LE_fall")
dev.off()

# -two kinds of flux partitioning for NEE:
# 1) _MR05: Reichstein et al. (2005)
# 2) _GL10: Lasslop et al. (2010)
MR05 <- ess_out
MR05_FP_names_filter <- names(MR05) %in% c("Reco_uStar", "GPP_uStar_f")
names(MR05)[MR05_FP_names_filter] <- c("Reco", "GPP")
pdf(paste(path_out, siteyear, "_NEE_uStar_f_MR05_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(MR05, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_f", NEE_sep = TRUE)
dev.off()
pdf(paste(path_out, siteyear, "_NEE_uStar_fall_MR05_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(MR05, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_fall", NEE_sep = TRUE)
dev.off()

GL10 <- ess_out
GL10_FP_names_filter <- names(GL10) %in% c("Reco_DT_uStar", "GPP_DT_uStar")
names(GL10)[GL10_FP_names_filter] <- c("Reco", "GPP")
pdf(paste(path_out, siteyear, "_NEE_uStar_f_GL10_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(GL10, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_f", NEE_sep = TRUE)
dev.off()
pdf(paste(path_out, siteyear, "_NEE_uStar_fall_GL10_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
plot_eddy(GL10, "NEE", "qc_NEE_forGF_UF", "qc_NEE_forGF_UF",
          flux_gf = "NEE_uStar_fall", NEE_sep = TRUE)
dev.off()

# EOF
