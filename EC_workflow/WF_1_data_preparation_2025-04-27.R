### Description ================================================================

# Eddy covariance workflow part 1/4 (https://github.com/lsigut/EC_workflow) This
# code primarily aims at data preparation for quality control (QC). Meteo data
# and EddyPro full output files are validated, accordingly formatted, merged and
# saved with documentation. All numeric values are rounded to a reasonable
# precision. Meteo variable names are remapped according to the requirements of
# openeddy and REddyProc packages. If Meteo contains replicates of the same
# variable, mean of all replicates and their QC is returned. It is expected that
# Meteo data underwent separate quality control and gap-filling (not in the
# scope of openeddy). Notice that especially gaps in incoming radiation (GR -
# global radiation) can have negative impact on the reliability and quality
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
attach_pkg("tibble") # required for tribble()

# Check if openeddy version conforms to requirements
if (packageVersion("openeddy") < package_version("0.0.0.9009"))
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9009'")

### Provide metadata and set file paths and arguments ==========================

# Load the site-year settings file
settings_file <- list.files(pattern = "settings", full.names = TRUE)
source(settings_file)

# Load the list of folder structure paths
# - automated, no input required if proposed folder structure is followed
paths <- make_paths()

# Required input files are meteo data and processed eddy covariance data

# Meteo data
# - automated, no input required if proposed folder structure is followed
# - Meteo groups checklist (30 min): 
#   1) Radiation: GRin, PARin, Rn
#   2) Temperature: Ta (from EC height), Ts (from the top soil layer) 
#   3) Moisture: Ma (from EC height), Ms (from the top soil layer) 
#   4) Precipitation: sumP
#   5) Heat flux: HF

# EddyPro data
# - automated, no input required if proposed folder structure is followed
# - folder can contain multiple EddyPro files to be merged by merge_eddy()
# - folder is expected to contain CSV (csv) files with EddyPro data to merge and
#   optionally also TXT (txt) files with their documentation; their file names 
#   should differ only in their file extension

# Timestamp of the computation
# - automated, will be included in file names
Tstamp <- format(Sys.time(), "%Y-%m-%d") 

### Load and format Meteo data =================================================

# read_eddy() reads single meteo CSV file including units (placed on the second 
# row below header) at a standardized path. 
mf <- list.files(paths$qc_input_meteo, pattern = "\\.[Cc][Ss][Vv]$",
                 full.names = TRUE)[1]
M <- read_eddy(mf, check.names = FALSE)

# Rename Meteo data variables as required by openeddy and REddyProc packages
# - other columns than those included in Met_mapping table are dropped
# - not available columns are reported and automatically initialized with NAs
M <- remap_vars(M, Met_mapping$workflow_varname, Met_mapping$Meteo_varname,
                regexp = TRUE, qc = "_qcode")

# strptime_eddy() rewrites the original varname of the timestamp column 
# - retain the original varname for documentation purposes
vars <- openeddy::varnames(M)

# timestamp requires conversion to POSIXct for validation
M$timestamp <- strptime_eddy(M$timestamp, format = meteo_format, 
                             allow_gaps = TRUE)

# reset original varnames 
openeddy::varnames(M) <- vars

# merge_eddy() assures that timestamp is complete and has defined range
M <- merge_eddy(list(M), start = start, end = end)

# Correct units
# - not included in correct() as it is Czechglobe specific formatting
# - it should not impact sites with other formatting
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
EP <- read_EddyPro(paths$qc_input_eddypro, start = start, end = end, skip = 1, 
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

# Change the timestamp formatting to default timestamp format
# - format = "%Y-%m-%d %H:%M" is enforced to simplify further processing steps
data$timestamp <- format(data$timestamp, format = "%Y-%m-%d %H:%M", tz = "GMT")

# Round the columns of numeric mode type double to 6 significant digits
data <- round_df(data)

# Reset units lost during merging
# - varnames are used only for documentation below
openeddy::units(data) <- c(openeddy::units(M), openeddy::units(EP[-1]))

# Save the merged data with documentation ======================================

# Set the name of merged output file
data_name_out <- name_merged(paths$qc_input_eddypro, siteyear)

# Save the merged Meteo and EddyPro data
write_eddy(data, file.path(paths$input_for_qc, data_name_out))

# Documentation of merged files
# - TXT files with the same name as the merged input CSV files are used if 
#   present in respective folders with input data
# - information about the Meteo remapping and session info are also included
# - the documentation file will not be overwritten if it already exists, this is
#   to avoid overwriting manually edited documentation; to overwrite it, check 
#   file content and delete it manually if safe
document_merged(data_name_out, paths$qc_input_eddypro, paths$qc_input_meteo, 
                paths$input_for_qc, Tstamp, 
                name, mail, M)

# EOF
