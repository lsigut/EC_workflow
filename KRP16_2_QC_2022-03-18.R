### Description ================================================================

# Eddy covariance workflow part 2/4 (https://github.com/lsigut/EC_workflow) This
# code primarily aims at data quality checking (QC) and preparation of input
# data for gap-filling. Data before, during and after QC are plotted and further
# statistics are computed. Intermediate results are plotted or reported in
# console to assist the user. Storage correction and correction of rotated
# vertical wind speed (w_rot) is optional and overwrites the original values of
# respective variables.
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
packages <- c("tidyverse", "ggplot2", "gridExtra", "reshape2")
invisible(lapply(packages, attach_pkg))

# Workflow is currently aligned only with specific package version
# - package version should be neither higher or lower
if (packageVersion("openeddy") < "0.0.0.9008")
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9008' and above")

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

# Do you want to perform manual data quality check? 
# - default interactive_session <- TRUE is currently recommended
# - it mainly affects manual QC step >check_manually()<
# - set FALSE if you want to simply rerun this workflow part  
# - if FALSE, informative plotting will be skipped
# - if FALSE, manual QC will be still used if created previously
interactive_session <- TRUE

# Do you want to apply storage correction to H, LE and CO2 flux?
# - currently only storage correction estimated by EddyPro using discrete 
#   (one point) approach is implemented
# - if TRUE, storage flux is added to the respective original flux
# - recommended for sites with short canopy, e.g. grasslands, wetlands, 
#   croplands
apply_storage <- TRUE

# Specify the type of coordinate rotation applied during EddyPro processing
# - supported rotation types: 
#   - "double": for double (2D) rotation
#   - "planar fit": for planar fit rotation
rotation_type <- "double"

# Specify the type of IRGA used in eddy covariance system
# - supported IRGA types: 
#   - "en_closed": both for closed and enclosed path systems 
#   - "open": for open path systems
IRGA_type <- "en_closed"

# Load the list of folder structure paths
# - automated, no input required if proposed folder structure is followed
paths <- structure_eddy()

# Specify the time shift (in seconds) to be applied to the date-time information
# in order to represent the center of averaging period
shift.by <- -900

# Timestamp of the computation
# - automated, will be included in file names
Tstamp <- format(Sys.time(), "%Y-%m-%d") 

# Set Fetch filter boundary ====================================================

# Edit the region of interest boundary for your site
# - here KRP boundary version 20160206
# - see ROI boundary concept description at https://github.com/lsigut/EC_workflow
boundary <- 
  c(453, 489, 469, 455, 444, 410, 375, 348, 86, 82, 78, 76, 74, 73, 72, 72, 73, 
    74, 76, 78, 81, 85, 91, 97, 106, 116, 114, 113, 131, 372, 496, 500, 507, 
    519, 531, 541, 555, 562, 565, 572, 584, 605, 633, 749, 863, 1012, 1128, 
    1098, 802, 863, 871, 903, 403, 360, 328, 303, 283, 486, 466, 451, 441, 412, 
    390, 373, 360, 350, 349, 356, 367, 381, 399, 422)

### Load and format data =======================================================

# Edit the path to required inputs - these are meteo and EddyPro data merged in
# the previous >data preparation< processing step to a single validated CSV file
data <- read_eddy("./Level 1/Post-processing/EddyProOutput/eddypro_KrP_full_output_2017-05-20T180735_adv_met.csv")

# Convert timestamp to POSIXct and shift the date-time information to represent 
# the center of averaging period which is required for reliable processing
data$timestamp <- strptime_eddy(data$timestamp, shift.by = shift.by)

### Plot data for visual precheck ==============================================

# Save plots of variables that can help identify problems with instruments
# - precheck_vars is an object defined within openeddy
# - choose only names available in data
precheck <- choose_avail(precheck_vars, names(data))

# Save plots of precheck variable to single pdf at paths$Precheck
save_precheck_plots(data, precheck, siteyear, Tstamp, paths$Precheck)

# Show dependence of w_unrot on wind direction with additional statistics
# - to reduce sensitivity of computed stats to outliers, median and mad is used
# - qrange is quantile range for ylim to reduce impact of outliers
#   - qrange has only visual effect, it does not affect computed statistics 
#   - if you do not want to limit y-axis, set qrange = NULL
ggsave(file.path(
  paths$WD_dependency,
  paste0(siteyear, "_w_unrot_WD_stats_", Tstamp, ".png")),
  ggplot_stats(data, "wind_dir", "w_unrot", circular = TRUE),
  type = "cairo-png", width = 297, height = 210, units = "mm")

# Show dependence of w_rot on wind direction with additional statistics
ggsave(file.path(
  paths$WD_dependency,
  paste0(siteyear, "_orig_w_rot_WD_stats_", Tstamp, ".png")),
  ggplot_stats(data, "wind_dir", "w_rot", circular = TRUE),
  type = "cairo-png", width = 297, height = 210, units = "mm")

# Check the overall w_rot median
# - this section is mainly relevant if rotation_type == "planar fit" and
#   evaluation of w wind component residual will be applied in the QC scheme
# - ideally, w_rot median should be close to zero as any substantial deviations
#   could lead to excessive flagging by wresid filter
# - this can be corrected by forcing the w_rot median to zero or by correcting
#   it by B0 planar fit coefficient derived from the data

# If applied, document here B0 used for correction from the folder:
# none
# file: none

# Edit w_rot_correction value that will be substracted from w_rot
# - if w_rot does not require correction, set "none"
# - if rotation_type == "double", w_rot_correction is not considered
w_rot_correction <- "none"

# Store original value of w_rot_correction for documentation purposes
# - once correction is applied, w_rot_correction is set to "none"
applied_w_rot_correction <- w_rot_correction

# Skip if w_rot_correction = "none"
if (!w_rot_correction == "none" && rotation_type == "planar fit") {
  data$w_rot <- data$w_rot - w_rot_correction 
  # Confirm the succesful w_rot correction (skip if w_rot_correction = "none")
  ggsave(file.path(
    paths$WD_dependency,
    paste0(siteyear, "_corrected_w_rot_WD_stats_", Tstamp, ".png")),
    ggplot_stats(data, "wind_dir", "w_rot", circular = TRUE),
    type = "cairo-png", width = 297, height = 210, units = "mm")
  w_rot_correction <- "none" # avoid rerunning of the correction by mistake
}

# Specify the vector of supported flux names used repeatedly within workflow 
fluxes <- c("Tau", "H", "LE", "NEE")

# Save flux time series precheck plots with distinguished QC along basic meteo 
# - SSITC is the standard "Foken flag" (e.g. qc_H) from EddyPro renamed by 
#   correct() within data_preparation workflow
save_flux_plots(data, "SSITC", siteyear, "%s_precheck", Tstamp, paths$Precheck, 
                fluxes)

### Extract flags of predefined tests/filters ==================================

# A set of filters is extracted that may or may not be useful to apply at a
# given site
# - it is suggested to experiment with the setup to get optimal results
# - applied filters in this workflow seem to be a practical combination but 
#   they should be further tested considering their flagging efficiency 

# Notes for wresid filter:
# - recommended only for short canopy sites with setting rotation == "double"
# - typically excessive flagging for sites with rotation == "planar fit"
# - needs further testing
QC <- extract_QC(data, rotation = rotation_type)
summary_QC(QC, names(QC))

# Save the results to the main data frame
data[names(QC)] <- QC 

### Combine flags representing mainly technical issues with fluxes =============

# Create a tibble with QC filter names to combine for each flux 
prelim <- tribble(
  ~Tau,              ~H,               ~LE,                ~NEE,
  "qc_Tau_SSITC",    "qc_H_SSITC",     "qc_LE_SSITC",      "qc_NEE_SSITC",
  "qc_SA_abslim",    "qc_SA_abslim",   "qc_SAGA_abslim",   "qc_SAGA_abslim",
  "qc_SA_spikesHF",  "qc_SA_spikesHF", "qc_SAGA_spikesHF", "qc_SAGA_spikesHF",  
  "qc_Tau_missfrac", "qc_H_missfrac",  "qc_LE_missfrac",   "qc_NEE_missfrac",
  "qc_Tau_scf",      "qc_H_scf",       "qc_LE_scf",        "qc_NEE_scf",
  NA,                "qc_H_var",       NA,                 NA,
  "qc_Tau_runs",     "qc_H_runs",      "qc_LE_runs",       "qc_NEE_runs",
  NA,                "qc_H_lowcov",    "qc_LE_lowcov",     "qc_NEE_lowcov",
  NA,                NA,               "qc_GA_LI7200",     "qc_GA_LI7200", 
  "qc_ALL_wresid",   "qc_ALL_wresid",  "qc_ALL_wresid",    "qc_ALL_wresid"
)

# Combine specified flags for given flux to produce preliminary flags
pre_res <- combn_prelim_QC(data, prelim)

### Apply flux interdependency =================================================

# Evaluate flux interdependency based on the preliminary flags
# - preliminary H flag is needed only if IRGA = "open", otherwise not used
# - preliminary Tau flag is not used
interdep <- interdep(pre_res$qc_LE_prelim, pre_res$qc_H_prelim, IRGA_type)
summary_QC(interdep, names(interdep))

# Save the results to the main data frame
data[names(interdep)] <- interdep

# In case of IRGA_type == "en_closed" LE is not affected by flux interdependency
qc_LE_interdep <- ifelse(IRGA_type == "open", "qc_LE_interdep", NA)

# Include flux interdependency among QC filter names to combine
# - flux interdependency does not affect Tau
prelim_ad <- tribble(
  ~Tau, ~H,              ~LE,            ~NEE,
  NA,   "qc_H_interdep", qc_LE_interdep, "qc_NEE_interdep"
)

prelim2 <- rbind(prelim, prelim_ad)

# Combine specified flags for given flux to produce preliminary flags
pre2_res <- combn_prelim_QC(data, prelim2)

if (interactive_session) {
  # Show effect of flux interdependency test on QC flags
  # - no effect on Tau flags (flux interdependency not defined for Tau) 
  # - no effect on LE flags if IRGA_type == "en_closed"
  plot_QC_summary(data, prelim2, cumul = FALSE)
  plot_QC_summary(data, prelim2, cumul = TRUE)
  
  # Update flux time series precheck plots according to prelim2 QC flags
  save_flux_plots(cbind(data, pre2_res), "prelim2", siteyear, "%s_precheck",
                  Tstamp, paths$Precheck, fluxes)
}

### Apply storage correction ===================================================

# Store original value of apply_storage for documentation purposes
# - once correction is applied, apply_storage is set to FALSE
strg_applied <- apply_storage

if (apply_storage) {
  # Storage flux is added to the respective original flux
  # - correction  overwrites the original values of respective variables 
  data$H <- add_st(data$H, data$H_strg)
  data$LE <- add_st(data$LE, data$LE_strg)
  data$NEE <- add_st(data$NEE, data$co2_strg)
  apply_storage <- FALSE # avoid rerunning of the correction by mistake
}

### Apply prelim2 filters to fluxes ============================================

# Produce flux columns with suffix "_orig" only with fluxes passing previous QC 
orig_fluxes <- sapply(fluxes, function(x) paste0(x, "_orig")) 
for (i in seq_along(fluxes)) {
  data[orig_fluxes[i]] <- apply_QC(data[, fluxes[i]], pre2_res[, i])
}

### Extract flags based on low frequency data despiking ========================
desp <- data[0]

# Low frequency data despiking is not applied for Tau
# - change qrange if needed
plot_precheck(data, "H_orig", qrange = c(0, 1), pch = 19)
desp$qc_H_spikesLF <- 
  despikeLF(cbind(data, pre2_res), var = "H", qc_flag = "qc_H_prelim2",
            name_out = "qc_H_spikesLF", var_thr = c(-200, 800))

plot_precheck(data, "LE_orig", qrange = c(0, 1), pch = 19)
desp$qc_LE_spikesLF <- 
  despikeLF(cbind(data, pre2_res), var = "LE", qc_flag = "qc_LE_prelim2",
            name_out = "qc_LE_spikesLF", var_thr = c(-200, 800))

plot_precheck(data, "NEE_orig", qrange = c(0.005, 0.995), pch = 19)
desp$qc_NEE_spikesLF <- 
  despikeLF(cbind(data, pre2_res), var = "NEE", qc_flag = "qc_NEE_prelim2",
            name_out = "qc_NEE_spikesLF", var_thr = c(-100, 100))

# Check the results
# - the test has 3 outcomes: 0 - OK; 2 - spike; NA - excluded from despiking
# - since NA means not checked, NA is thus interpreted as flag 0
summary_QC(desp, names(desp))

# Save the results to the main data frame
data[names(desp)] <- desp 

### Extract fetch filter =======================================================

# Fetch filter is not applied for Tau
data$qc_ALL_fetch70 <- fetch_filter(
  data, "x_70perc", "wind_dir", boundary, "qc_ALL_fetch70")
summary_QC(data, "qc_ALL_fetch70")

### Combine QC flags before manual QC ==========================================

# Include despiking and fetch filter among QC filter names to combine
# - despiking and fetch filter is not applied for Tau
prelim2_ad <- tribble(
  ~Tau, ~H,               ~LE,                 ~NEE,
  NA,   "qc_H_spikesLF",  "qc_LE_spikesLF",    "qc_NEE_spikesLF",
  NA,   "qc_ALL_fetch70", "qc_ALL_fetch70",    "qc_ALL_fetch70"
)

prelim3 <- rbind(prelim2, prelim2_ad)

# Combine specified flags for given flux to produce preliminary flags
pre3_res <- combn_prelim_QC(data, prelim3)

if (interactive_session) {
  # Show effect of all filters before manual QC
  plot_QC_summary(data, prelim3, cumul = FALSE)
  plot_QC_summary(data, prelim3, cumul = TRUE)

  # Update flux time series precheck plots according to prelim3 QC flags
  save_flux_plots(cbind(data, pre3_res), qc_suffix = "prelim3", siteyear, 
                  "%s_precheck", Tstamp, paths$Precheck, fluxes)
}

### Apply prelim3 filters to fluxes ============================================

# Produce flux columns with suffix "_orig" only with fluxes passing previous QC 
for (i in seq_along(fluxes)) {
  data[orig_fluxes[i]] <- apply_QC(data[, fluxes[i]], pre3_res[, i])
}

### Run manual quality control =================================================

# Check all fluxes specified by 'vars' already screened by QC scheme above
# - at any point, intermediate progress can be saved using option '6. finalize'
#   and inserting 'y' to following dialog: 
#   "save current progress to file at 'path'?"
# - progress after the flagging of last flux should be saved in order to fully 
#   reproduce the manual flags when rerunning the code next time
# - if saving to a file is omitted, results are still saved to 'man' object
#   and later also to data frame 'data' but it does not allow easy rerunning
# - returned timestamp is removed in 'man' but kept in saved CSV file
# - if not interactive and no manual QC found: NULL returned 
man <- check_manually(cbind(data, pre3_res), paths$Quality_checking, 
                      vars = data.frame(
                        x = fluxes,
                        y = c("PAR", "Rn", "Rn", "PAR"),
                        z = c("wind_speed", "LE_orig" , "H_orig", "Tair")
                      ), 
                      qc_prefix = "qc_", qc_suffix = "_prelim3", 
                      interactive_session, siteyear)[-1]
summary_QC(man, names(man))

# Save the results to the main data frame
data[names(man)] <- man 

### Combine QC flags for last checkout =========================================

# Note: whole section can be run iteratively

# Consider only existing manual flags
# - there might be variables without manual QC or 'man' can be NULL
man_names <- set_man_names(fluxes, man)

# Include manual QC among QC filter names to combine
prelim3_ad <- tribble(
  ~Tau,          ~H,          ~LE,          ~NEE,
  man_names$Tau, man_names$H, man_names$LE, man_names$NEE
)

prelim4 <- rbind(prelim3, prelim3_ad)

# Combine specified flags for given flux to produce preliminary flags
pre4_res <- combn_prelim_QC(data, prelim4)

if (interactive_session) {
  # Update flux time series precheck plots according to prelim4 QC flags
  save_flux_plots(cbind(data, pre4_res), qc_suffix = "prelim4", siteyear, 
                  "%s_precheck", Tstamp, paths$Precheck, fluxes)
}

# Apply prelim4 filters to fluxes
# - produce flux columns with suffix "_orig" only with fluxes passing previous QC 
for (i in seq_along(fluxes)) {
  data[orig_fluxes[i]] <- apply_QC(data[, fluxes[i]], pre4_res[, i])
}

man <- check_manually(cbind(data, pre4_res), paths$Quality_checking, 
                      vars = data.frame(
                        x = fluxes,
                        y = c("PAR", "Rn", "Rn", "PAR"),
                        z = c("wind_speed", "LE_orig" , "H_orig", "Tair")
                      ), 
                      qc_prefix = "qc_", qc_suffix = "_prelim4", 
                      interactive_session, siteyear)[-1]
summary_QC(man, names(man))

# Save the results to the main data frame
data[names(man)] <- man 

### Combine QC flags to final flag for gap-filling =============================

# Consider only existing manual flags
# - there might be variables without manual QC or 'man' can be NULL
man_names <- set_man_names(fluxes, man)

# Include manual QC among QC filter names to combine
prelim3_ad <- tribble(
  ~Tau,          ~H,          ~LE,          ~NEE,
  man_names$Tau, man_names$H, man_names$LE, man_names$NEE
  )

forGF <- rbind(prelim3, prelim3_ad)

# Combine specified flags for given flux to produce final forGF flags
forGF_res <- combn_prelim_QC(data, forGF)

# Save the results to the main data frame
data[names(forGF_res)] <- forGF_res 

### Apply forGF filters to fluxes ==============================================

# Produce flux columns with suffix "_orig" only with fluxes passing final QC 
for (i in seq_along(fluxes)) {
  data[orig_fluxes[i]] <- apply_QC(data[, fluxes[i]], forGF_res[, i])
}

### Produce QC summary and save and plot the results ===========================

# Show the percentage of fluxes flagged by individual filters and their 
# cumulative effect
list_QC_ind <- lapply(names(forGF), function(x) 
  print(summary_QC(data, na.omit(pull(forGF, x)))))
list_QC_cum <- lapply(names(forGF), function(x) 
  print(summary_QC(data, na.omit(pull(forGF, x)), cumul = TRUE)))
list_QC <- c(list_QC_ind, list_QC_cum)
names(list_QC) <- paste0(names(forGF), rep(c("", "_cumulative"), each = 4))

# Save each element of list_QC to separate file in QC_summary folder
for (i in names(list_QC)) {
  write.csv(list_QC[[i]], 
            file.path(
              paths$QC_summary,
              paste0(siteyear, "_QC_summary_", i, "_", Tstamp, ".csv")))
}

# Save the QC summary results as plots
save_QC_summary_plots(data, forGF, paths$QC_summary, siteyear, Tstamp)

### Plot the quality checked data ==============================================

# Save the plots that show the data with QC flag used for gap-filling
save_flux_plots(data, qc_suffix = "forGF", siteyear, "forGF_QC_%s", Tstamp, 
                paths$Quality_checking, fluxes)

### Write QC results to a file =================================================

# Allow data to be usable also after saving (no timestamp reformatting)
save_data <- data

# Correct timestamp to its original state
save_data$timestamp <- save_data$timestamp - shift.by 

# Remap variable names and filter columns needed for online gap-filling tool:
OT_in <- set_OT_input(save_data,
                      c("qc_NEE_forGF", "NEE", "qc_LE_forGF", "LE",
                        "qc_H_forGF", "H", "GR", "Tair", "Tsoil", "RH",
                        "VPD", "qc_Tau_forGF", "ustar"))

# Save the standardized Online Tool/REddyProc input to '.txt' file
# - rounding not required as QC flags are integers and the rest was rounded
#   in data preparation step
write_eddy(OT_in, 
           file.path(
             paths$Input_for_GF,
             paste0(siteyear, "_", Tstamp, ".txt")),
           quote = FALSE, sep = "\t")

# Change the timestamp formatting to its original state:
save_data$timestamp <- format(save_data$timestamp, format = "%Y-%m-%d %H:%M", 
                              tz = "GMT")

# Save the full output
# - rounding not required as QC flags are integers and the rest was rounded
#   in data preparation step
write_eddy(save_data, 
           file.path(
             paths$Quality_checking,
             paste0(siteyear, "_forGF_QC_full_output_", Tstamp, ".csv")))

# Choose the most important variables to later combine with gap-filling results
# - essential_vars_QC is an object defined within openeddy
# - use only those that are available
essentials <- choose_avail(essential_vars_QC, names(save_data))

# Save essential output
write_eddy(save_data[essentials],
           file.path(
             paths$Input_for_GF,
             paste0(siteyear, "_forGF_QC_essentials_", Tstamp, ".csv")))

# Save documentation information about executed quality control 
names(forGF) <- names(forGF_res)
document_QC(Tstamp, name, mail, strg_applied, forGF,
            paths$Quality_checking, siteyear)
#EOF
