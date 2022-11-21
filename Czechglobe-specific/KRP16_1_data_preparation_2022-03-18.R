### Description ================================================================

# Eddy covariance workflow part 1/4 (https://github.com/lsigut/EC_workflow).
# This code primarily aims at data preparation for quality control (QC). Meteo
# data, EddyPro full output files and instruments QC are validated, accordingly
# formatted, merged and saved. All numeric values are rounded to a reasonable
# precision. Meteo variable names are remapped according to Czechglobe
# standards. If Meteo contains replicates of the same variable, mean of all
# replicates and their QC is returned. Radek Czerny´s quality control of
# instruments includes "SA_instrum" and "GA_instrum" columns from SmartFlux
# files in "./Level 1/Post-processing/Preliminary" folder.
#
# For documentation of EddyPro variable names see:
# https://www.licor.com/env/support/EddyPro/topics/output-files-full-output.html
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
attach_pkg("tidyverse") # required for tribble()

# Workflow is currently aligned only with specific package version
# - package version should be neither higher or lower
if (packageVersion("openeddy") < "0.0.0.9006")
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9006' and above")

### Provide metadata and set file paths and arguments ==========================

# Contact information
name <- "Ladislav Sigut" # person that performed processing
mail <- "sigut.l@czechglobe.cz" # mail of the person that performed processing

# Edit the siteyear
# - included in folder and file names
siteyear <- "KRP16"

# Edit the year
# - used to define the time extent of used data, i.e. full year
year <- 2016 

# Required input files are meteo data, processed eddy covariance data and a file
# with QC flags marking instrument malfunctions

# Edit the folder name including meteo data
# - this folder is created by the user and can contain multiple MeteoDBS files
#   that will be merged along regular timestamp - check merge_eddy()
# - the path is difficult to standardize because it might depend on processing
#   versioning (e.g. .../EddyProOutput/Run2/Meteo data/)
# - MeteoDBS groups checklist (30 min): 
#   1) Radiation: GRin, PARin, Rn
#   2) Temperature: Ta (from EC height), Ts (from the top soil layer) 
#   3) Moisture: Ma (from EC height), Ms (from the top soil layer) 
#   4) Precipitation: sumP
#   5) Heat flux: HF
Meteo_path <- "./Level 1/Post-processing/EddyPro Output/Meteo data/"

# Edit the folder name including EddyPro data
# - this folder is created by the user and can contain multiple EddyPro files
#   that will be merged along regular timestamp - check merge_eddy()
# - the path is difficult to standardize because it might depend on processing
#   versioning (e.g. .../EddyProOutput/Run2/EddyPro data/)
EP_path <- "./Level 1/Post-processing/EddyPro Output/EddyPro data/"

# Edit the file name including flags marking instrument malfunctions 
instrum_input <- "./Level 1/Post-processing/Preliminary/RAJ_2020_SmartFlux.csv"

# Edit the folder name for output files (data)
# - the path is difficult to standardize because it might depend on processing
#   versioning (e.g. .../EddyProOutput/Run2)
out_path <- "./Level 1/Post-processing/EddyPro Output/"

# Timestamp of the computation
# - automated, will be included in file names
Tstamp <- format(Sys.time(), "%Y-%m-%d") 

# Set MeteoDBS mapping =========================================================

# Edit MeteoDBS mapping according to available variables at given site
# - regular expressions (?regex) can be used to select replicates of given 
#   variable (e.g. Ts) that will be averaged by remap_vars()
# - typical set of variables: global radiation (GR), photosynthetically active
#   radiation (PAR), Net radiation (Rn), air temperature (Tair), soil 
#   temperature (Tsoil), relative humidity (RH), vapor pressure dificit (VPD),
#   soil water content (SWC), precipitation (P), soil heat flux (G)
# - variables required by openeddy: PAR, Rn, Tair, Tsoil, VPD, P (only PAR
#   really needed for optimal setup, the rest can be initialized with NA values)
# - variables required by REddyProc: GR, Tair, Tsoil, VPD (and/or RH)	
Met_mapping <- tribble(
  ~MeteoDBS_varname, ~workflow_varname,
  "date/time",       "timestamp",
  "GRin",            "GR",
  "PARin",           "PAR",
  "^Rn",             "Rn",
  "TaKA02.0",        "Tair",
  "TsKO0.05",        "Tsoil",
  "RHKA02.0",        "RH",
  "VPDKA02.0",       "VPD",
  "sumP",            "P"
)

### Load and format Meteo data =================================================

# read_MeteoDBS() reads all meteo files at Czechglobe MeteoDBS format at given
# path and merges them together. The expectation is that files represent meteo
# variables for given site and different periods. Function merges them
# vertically (along generated complete timestamp). Original column names are
# retained for reliable variable remapping. See details in utilities.R file.
M <- read_MeteoDBS(Meteo_path, start = year, end = year)

# Rename Meteo data variables to FLUXNET standard
# - other columns than those required are removed
# - not available columns are reported and automatically initialized with NAs
M <- remap_vars(M, Met_mapping$workflow_varname, Met_mapping$MeteoDBS_varname,
                regexp = TRUE, qc = "_qcode")

# Correct units
# - not included in correct() as it is Czech(globe) specific
openeddy::units(M) <- gsub("st. ", "deg", openeddy::units(M))

### Load and format EddyPro full output ========================================

# Load the EddyPro files (untouched originals) and bind them together

# read_EddyPro() reads all EddyPro files at given path and merges them together.
# The expectation is that files represent variables for given site and different
# periods. Function merges them vertically (along generated complete timestamp).
# Original column names are retained for reliable variable remapping.

# - set correct skip parameter (number of lines above header in the input)
# - set correct file encoding (e.g. fileEncoding = "UTF-8") 
# - make sure that all files have the same encoding
# - EddyPro output contains duplicated column names "co2_mean" and "h2o_mean"
#   that are respective concentrations in different units (their correct match
#   across files is not automated)
# - older versions of EddyPro had more duplicated column names in the output
# - old EddyPro column name "max_speed" is corrected to "max_wind_speed" 
EP <- read_EddyPro(EP_path, start = year, end = year, skip = 1, 
                   fileEncoding = "UTF-8")

# Correct column names
names(EP) <- correct(names(EP))

# Correct units
openeddy::units(EP) <- correct(openeddy::units(EP), attr = "units")
# - check whether µ was correctly interpreted due to the encoding issues

# Rename EddyPro variables RH and VPD to avoid duplication with meteo variables
# measured by slow response sensors (slow response sensors are more reliable)
# - not included in correct() as it is context specific 
names(EP)[names(EP) %in% c("RH", "VPD")] <- c("RH_EddyPro", "VPD_EddyPro")

# Merge Meteo, EddyPro data and instrum flags ==================================

# Simple merge() works because both M and EP have continuous equidistant
# timestamp covering whole year. Assumption is that M and EP do not have other
# shared column names than timestamp. In case user decides to shorten timestamp
# of one of the data frames by "start" and "end" arguments, the merge() output
# will be also shortened.
data <- merge(M, EP)
nrow(data) # if 0 rows, M & EP shared more columns than just timestamp

# Load Instrum flags from SmartFlux file in "Preliminary" folder
# - QC columns currently provided by Radek in Excel file
# - file needs to be resaved to CSV with corrected column names (sometimes 
#   missing) 
instrum <- read_eddy(instrum_input)[c("timestamp", "SA_instrum", "GA_instrum")]

# Rename them by adding "qc_" prefix to follow the QC naming strategy
names(instrum)[match(c("SA_instrum", "GA_instrum"), names(instrum))] <- 
  c("qc_SA_instrum", "qc_GA_instrum")

# Convert timestamp to POSIXct format
instrum$timestamp <- strptime_eddy(instrum$timestamp, allow_gaps = TRUE)

# Merge with data to assure correct lines and to fill timestamp gaps if present
data <- merge(data, instrum, by = "timestamp", all.x = TRUE)

# Correct Instrum flags that are in range 2-3
# - flag 0 is not reported (NA values instead)
# - flag 1 is not used (hard flag does not resolve this level)
# - flag 2 values are hard flags that will be used for screening
# - flag 3 in qc_SA_instrum for ts_var>2 and in qc_GA_instrum for co2_var>50
#   - it is corrected to flag 0 as variance tests can be applied separately
data$qc_SA_instrum[!(data$qc_SA_instrum %in% 2)] <- 0L
data$qc_GA_instrum[!(data$qc_GA_instrum %in% 2)] <- 0L

# Change the timestamp formatting to its original state
data$timestamp <- format(data$timestamp, format = "%Y-%m-%d %H:%M", tz = "GMT")

# Round the columns of numeric mode type double to 6 significant digits
data <- round_df(data)

# Reset units lost during merging
# - varnames are used only for documentation below
openeddy::units(data) <- c(openeddy::units(M), openeddy::units(EP[-1]), 
                           openeddy::units(instrum)[-1])

# Save the merged data with documentation ======================================

# Names of merged output files
# - data in CSV files and documentation in TXT files
data_name_out <- list.files(EP_path)
data_name_out <- grep("[.][Cc][Ss][Vv]$", data_name_out, value = TRUE)
if (length(data_name_out) == 1) {
  data_name_out <- gsub("[.][Cc][Ss][Vv]$", "_met_instrum.csv", data_name_out)
} else {
  data_name_out <- paste0("eddypro_", siteyear, 
                      "_full_output_merged_adv_met_instrum.csv")
}
docu_name_out <- gsub("[.][Cc][Ss][Vv]$", "\\.txt", data_name_out)

# Save the merged EddyPro, Meteo and instrum flags
write_eddy(data, paste0(out_path, "/", data_name_out))

# Combine documentation
EP_names <- list.files(EP_path, full.names = TRUE)
EP_names <- grep("[.][Cc][Ss][Vv]$", EP_names, value = TRUE)
M_names <- list.files(Meteo_path, full.names = TRUE)
M_names <- grep("[.][Cc][Ss][Vv]$", M_names, value = TRUE)

docu_name_in <- list.files(c(EP_path, Meteo_path), full.names = TRUE)
docu_name_in <- grep("[.][Tt][Xx][Tt]$", docu_name_in, value = TRUE)

# The documentation file will not be overwritten if it already exists
# - this is to avoid overwriting manually edited documentation
# - to overwrite it, check file content and delete it manually if safe
if (docu_name_out %in% list.files(out_path)) {
  message("Combined documentation already exists")
} else {
  writeLines(c(paste0(Tstamp, ":"),
               paste0("Files merged by ", name, " (", mail, ")"),
               "",
               "Merged files:", 
               M_names,
               EP_names,
               instrum_input,
               "",
               "Variables from meteo database remapped to:",
               paste(names(M), varnames(M), sep = " = ", collapse = "\n"), 
               "", 
               combine_docu(docu_name_in),
               "Information about the R session:",
               capture.output(sessionInfo())), 
             paste0(out_path, "/", docu_name_out), sep = "\n")
}

# EOF
