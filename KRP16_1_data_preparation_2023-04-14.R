### Description ================================================================

# Eddy covariance workflow part 1/4 (https://github.com/lsigut/EC_workflow) This
# code primarily aims at data preparation for quality control (QC). Meteo data
# and EddyPro full output files are validated, accordingly formatted, merged and
# saved with documentation. All numeric values are rounded to a reasonable
# precision. Meteo variable names are remapped according to the requirements of
# openeddy and REddyProc packages. If Meteo contains replicates of the same
# variable, mean of all replicates and their QC is returned. It is expected that
# Meteo data underwent separate quality control and gap-filling (not in the
# scope of openeddy). Notice that especially gaps in incoming radiation (PAR or
# shortwave radiation) can have negative impact on the reliability and quality
# of the final products.
#
# You can find example data set at https://doi.org/10.5281/zenodo.6631498
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

# Check if openeddy version conforms to requirements
if (packageVersion("openeddy") == package_version("0.0.0.9006"))
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9006'")

### Provide metadata and set file paths and arguments ==========================

# Contact information
name <- "Ladislav Sigut" # person that performed processing
mail <- "sigut.l@czechglobe.cz" # mail of the person that performed processing

# Edit the siteyear
# - included in folder and file names
siteyear <- "KRP16"

# Edit the start and end of time period to post-process
# - start <- 2016; end <- 2016: to assure complete year
# - start <- NULL; end <- NULL: use timestamp extent from input files
# - start <- "2016-02-01 10:00:00"; end <- "2017-02-01 10:00:00": arbitrary period
start <- 2016
end <- 2016

# Load the list of folder structure paths
# - automated, no input required if proposed folder structure is followed
paths <- structure_eddy()

# Required input files are meteo data and processed eddy covariance data

# Meteo data
# - automated, no input required if proposed folder structure is followed
# - folder can contain multiple MeteoDBS files to be merged by merge_eddy()
# - only CSV (csv) files found in the folder will be used
# - MeteoDBS groups checklist (30 min): 
#   1) Radiation: GRin, PARin, Rn
#   2) Temperature: Ta (from EC height), Ts (from the top soil layer) 
#   3) Moisture: Ma (from EC height), Ms (from the top soil layer) 
#   4) Precipitation: sumP
#   5) Heat flux: HF
Meteo_path <- paths$qc_input_meteo

# EddyPro data
# - automated, no input required if proposed folder structure is followed
# - folder can contain multiple EddyPro files to be merged by merge_eddy()
# - folder is expected to contain CSV (csv) files with EddyPro data to merge and
#   optionally also TXT (txt) files with their documentation; their file names 
#   should differ only in their file extension
EP_path <- paths$qc_input_eddypro

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
# retained for reliable variable remapping.
M <- read_MeteoDBS(Meteo_path, start = start, end = end)

# Rename Meteo data variables as required by openeddy and REddyProc packages
# - other columns than those included in Met_mapping table are dropped
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
# - included EddyPro output contains duplicated column names "co2_mean" and 
#   "h2o_mean" that are respective concentrations in different units (their
#   correct match across files is not automated)
# - older versions of EddyPro had more duplicated column names in the output
# - old EddyPro column name "max_speed" is corrected to "max_wind_speed" 
EP <- read_EddyPro(EP_path, start = start, end = end, skip = 1, 
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

# Merge Meteo and EddyPro data =================================================

# Simple merge() works because both M and EP have continuous equidistant
# timestamp covering whole year. Assumption is that M and EP do not have other
# shared column names than timestamp. In case user decides to shorten timestamp
# of one of the data frames by "start" and "end" arguments, the merge() output
# will be also shortened.
data <- merge(M, EP)
nrow(data) # if 0 rows, M & EP shared more columns than just timestamp

# Change the timestamp formatting to its original state
data$timestamp <- format(data$timestamp, format = "%Y-%m-%d %H:%M", tz = "GMT")

# Round the columns of numeric mode type double to 6 significant digits
data <- round_df(data)

# Reset units lost during merging
# - varnames are used only for documentation below
openeddy::units(data) <- c(openeddy::units(M), openeddy::units(EP[-1]))

# Save the merged data with documentation ======================================

# Set the name of merged output file
data_name_out <- name_merged(EP_path, siteyear)

# Save the merged Meteo and EddyPro data
write_eddy(data, file.path(paths$input_for_qc, data_name_out))

# Documentation of merged files
# - TXT files with the same name as the merged input CSV files are used if 
#   present in respective folders with input data
# - information about the Meteo remapping and session info are also included
# - the documentation file will not be overwritten if it already exists, this is
#   to avoid overwriting manually edited documentation; to overwrite it, check 
#   file content and delete it manually if safe
document_merged(data_name_out, EP_path, Meteo_path, paths$input_for_qc, Tstamp, 
                name, mail, M)

# EOF
