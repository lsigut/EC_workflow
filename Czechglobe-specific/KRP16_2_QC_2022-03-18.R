### Description ================================================================

# Eddy covariance workflow part 2/4 (https://github.com/lsigut/EC_workflow)
# This code primarily aims at data quality checking (QC) and preparation of
# input data for gap-filling. EddyPro and Meteo data inputs are merged and
# reformatted. If Meteo contains replicates of the same variable, mean of all
# replicates and their QC is returned. Data before and after QC are plotted and
# further statistics are computed. Storage correction and correction of rotated
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
if (packageVersion("openeddy") < "0.0.0.9007")
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9007' and above")

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

# BKF boundary version 20160206
# - see ROI boundary concept description at https://github.com/lsigut/EC_workflow
boundary <- 
  c(453, 489, 469, 455, 444, 410, 375, 348, 86, 82, 78, 76, 74, 73, 72, 72, 73, 
    74, 76, 78, 81, 85, 91, 97, 106, 116, 114, 113, 131, 372, 496, 500, 507, 
    519, 531, 541, 555, 562, 565, 572, 584, 605, 633, 749, 863, 1012, 1128, 
    1098, 802, 863, 871, 903, 403, 360, 328, 303, 283, 486, 466, 451, 441, 412, 
    390, 373, 360, 350, 349, 356, 367, 381, 399, 422)

### Load and format data =======================================================

# Required inputs are meteo, EddyPro and instrument QC data merged in the 
# previous >data preparation< processing step to a single validated CSV file
data <- read_eddy("./Level 1/Post-processing/EddyPro Output/eddypro_KrP_full_output_2017-05-20T180735_adv_met.csv")

# Convert timestamp to POSIXct and shift the date-time information to represent 
# the center of averaging period which is required for reliable processing
data$timestamp <- strptime_eddy(data$timestamp, shift.by = shift.by)

### Plot data for visual precheck ==============================================

# Display variables that can help identify problems with instruments
precheck <- c("u_rot", "v_rot", "w_unrot", "w_rot", "sonic_temperature", 
              "max_wind_speed",
              "Tau", "ustar", "H", "LE", "NEE", 
              "u_var", "v_var", "w_var", "ts_var", "h2o_var", "co2_var",
              "rand_err_Tau", "rand_err_H", "rand_err_LE", "rand_err_NEE",
              "Tau_scf", "H_scf", "LE_scf", "co2_scf",
              "u_spikes", "v_spikes", "w_spikes", "ts_spikes", "co2_spikes", 
              "h2o_spikes",
              "H_strg", "LE_strg", "co2_strg",
              "h2o_v_adv", "co2_v_adv",
              "co2_mixing_ratio", "h2o_mixing_ratio", 
              "co2_time_lag", "h2o_time_lag",
              "x_peak", "x_70perc",
              "mean_value_RSSI_LI_7200", "co2_signal_strength_7200_mean",
              "h2o_signal_strength_7200_mean", "flowrate_mean")

# Make sure that precheck names are available in data
# - message listing such names will be printed if any not available
precheck <- choose_avail(names(data), precheck)

# Save precheck plots to single pdf
# - qrange is quantile range for ylim to reduce impact of outliers 
pdf(file.path(
  paths$Precheck, 
  paste0(siteyear, "_auxiliary_precheck_", Tstamp, ".pdf")), 
  width = 11.00, height = 8.27)
invisible(lapply(precheck, plot_precheck, x = data, qrange = c(0.005, 0.995)))
dev.off()

# Show dependence of w_unrot on wind direction with additional statistics
# - to reduce sensitivity of computed stats to outliers, median and mad is used
# - qrange is quantile range for ylim to reduce impact of outliers
#   - qrange has only visual effect, it does not affect computed statistics 
#   - if you do not want to limit y-axis, set qrange = NULL or qrange = c(0, 1)
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

# Save flux time series precheck plots with distinguished QC along basic meteo 
for (i in c("Tau", "H", "LE", "NEE")) {
  pdf(file.path(
    paths$Precheck,
    paste0(siteyear, "_", i, "_precheck_", Tstamp, ".pdf")), 
    width = 11.00, height = 8.27)
  plot_eddy(data, i, paste0("qc_", i, "_SSITC"), paste0("qc_", i, "_SSITC"))
  dev.off()
}

### Extract flags of predefined tests/filters ==================================

### User should decide how wresid filter is handled 
#     - only for short canopy sites with setting rotation == "double"?
#     - excessive flagging for sites with rotation == "planar fit"?
#     - needs further testing, possibly rotation == "double" for all sites?
#     - is the excessive flagging for "planar fit" due to insufficient rotation?
QC <- extract_QC(data, rotation = rotation_type)
summary_QC(QC, names(QC))

# Save the results to the main data frame
data[names(QC)] <- QC 

### Extract user-specified flags ===============================================

# H fluxes are not reliable if sonic temperature variance is above 2
data$qc_H_thr_tsvar <- apply_thr(data$ts_var, c(2, 2), flag = "higher")

# Fluxes that form runs of equal values are likely spurious
data$qc_Tau_runs <- flag_runs(data$Tau)
data$qc_H_runs <- flag_runs(data$H)
data$qc_LE_runs <- flag_runs(data$LE)
data$qc_NEE_runs <- flag_runs(data$NEE)

# Fluxes with too low covariance are likely spurious
data$qc_H_lowcov <- apply_thr(data$H, c(-0.005, 0.005), flag = "between")
data$qc_LE_lowcov <- apply_thr(data$LE, c(-0.005, 0.005), flag = "between")
data$qc_NEE_lowcov <- apply_thr(data$NEE, c(-0.005, 0.005), flag = "between")

### Combine flags representing mainly technical issues with fluxes =============

# Instrum flags represent Radek´s manual check of instrument malfunctions
# - included during previous data preparation step

# Radek: you can remove safely all flags 2 as it covers periods of frozen sonic,
# maintenance, spikes >1% and low signal strength. Note that if SA_instrum is 2
# all fluxes are affected but concentrations are OK and if GA_instrum is 2 you
# can’t use information from gas analyzer but you can use sonic output (Tau, H,
# ustar, wind_speed, wind_dir)

# Create a tibble with QC filter names to combine for each flux 
prelim <- tribble(
  ~Tau,              ~H,               ~LE,                ~NEE,
  "qc_Tau_SSITC",    "qc_H_SSITC",     "qc_LE_SSITC",      "qc_NEE_SSITC",
  "qc_SA_abslim",    "qc_SA_abslim",   "qc_SAGA_abslim",   "qc_SAGA_abslim",
  "qc_SA_spikesHF",  "qc_SA_spikesHF", "qc_SAGA_spikesHF", "qc_SAGA_spikesHF",  
  "qc_Tau_missfrac", "qc_H_missfrac",  "qc_LE_missfrac",   "qc_NEE_missfrac",
  "qc_Tau_scf",      "qc_H_scf",       "qc_LE_scf",        "qc_NEE_scf",
  "qc_SA_instrum",   "qc_SA_instrum",  "qc_SA_instrum",    "qc_SA_instrum",
  NA,                NA,               "qc_GA_instrum",    "qc_GA_instrum",
  NA,                "qc_H_thr_tsvar", NA,                 NA,
  "qc_Tau_runs",     "qc_H_runs",      "qc_LE_runs",       "qc_NEE_runs",
  NA,                "qc_H_lowcov",    "qc_LE_lowcov",     "qc_NEE_lowcov",
  "qc_ALL_wresid",   "qc_ALL_wresid",  "qc_ALL_wresid",    "qc_ALL_wresid"
)

# Combine specified flags for given flux to produce preliminary flags
pre_res <- sapply(names(prelim), 
                  function(x) combn_QC(data, na.omit(pull(prelim, x))))
pre_res <- as.data.frame(pre_res)

### Apply flux interdependency =================================================

# Evaluate flux interdependency based on the preliminary flags
# - preliminary H flag is needed only if IRGA = "open", otherwise not used
# - preliminary Tau flag is not used
interdep <- interdep(pre_res$LE, pre_res$H, IRGA_type)
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
# - informative naming is useful during later despiking stage
pre2_res <- sapply(names(prelim2), 
                   function(x) combn_QC(data, na.omit(pull(prelim2, x))))
pre2_res <- as.data.frame(pre2_res)
names(pre2_res) <- paste0("qc_", names(pre2_res), "_prelim2")

# Show effect of flux interdependency test on QC flags
# - no effect on Tau flags (flux interdependency not defined for Tau) 
# - no effect on LE flags if IRGA_type == "en_closed"
if (interactive_session) {
  gridExtra::grid.arrange(grobs = lapply(names(prelim2), function(x) 
    summary_QC(data, na.omit(pull(prelim2, x)), plot = TRUE, flux = x)), 
    nrow = 2)
  gridExtra::grid.arrange(grobs = lapply(names(prelim2), function(x) 
    summary_QC(data, na.omit(pull(prelim2, x)), 
               plot = TRUE, cumul = TRUE, flux = x)), nrow = 2)
}

# Update flux time series precheck plots according to prelim2 QC flags
if (interactive_session) {
  for (i in c("Tau", "H", "LE", "NEE")) {
    pdf(file.path(
      paths$Precheck,
      paste0(siteyear, "_", i, "_precheck_", Tstamp, ".pdf")), 
      width = 11.00, height = 8.27)
    plot_eddy(cbind(data, pre2_res), i, 
              paste0("qc_", i, "_prelim2"), paste0("qc_", i, "_prelim2"))
    dev.off()
  }
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

### Extract flags based on low frequency data despiking ========================
desp <- list(qc_H_spikesLF = NULL, qc_LE_spikesLF = NULL, qc_NEE_spikesLF = NULL)

# Low frequency data despiking is not applied for Tau
plot_precheck(data, "H", qrange = c(0.006, 0.99), pch = 19)
desp$qc_H_spikesLF <- 
  despikeLF(cbind(data, pre2_res), var = "H", qc_flag = "qc_H_prelim2",
            name_out = "qc_H_spikesLF", var_thr = c(-200, 800))

plot_precheck(data, "LE", pch = 19)
desp$qc_LE_spikesLF <- 
  despikeLF(cbind(data, pre2_res), var = "LE", qc_flag = "qc_LE_prelim2",
            name_out = "qc_LE_spikesLF", var_thr = c(-200, 800))

plot_precheck(data, "NEE", pch = 19)
desp$qc_NEE_spikesLF <- 
  despikeLF(cbind(data, pre2_res), var = "NEE", qc_flag = "qc_NEE_prelim2",
            name_out = "qc_NEE_spikesLF", var_thr = c(-100, 100))
desp <- as.data.frame(desp)

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
# - informative naming is useful during later manual QC stage
pre3_res <- sapply(names(prelim3), 
                   function(x) combn_QC(data, na.omit(pull(prelim3, x))))
pre3_res <- as.data.frame(pre3_res)
names(pre3_res) <- paste0("qc_", names(pre3_res), "_prelim3")

# Show effect of all filters before manual QC
if (interactive_session) {
  gridExtra::grid.arrange(grobs = lapply(names(prelim3), function(x) 
    summary_QC(data, na.omit(pull(prelim3, x)), plot = TRUE, flux = x)), 
    nrow = 2)
  gridExtra::grid.arrange(grobs = lapply(names(prelim3), function(x) 
    summary_QC(data, na.omit(pull(prelim3, x)), 
               plot = TRUE, cumul = TRUE, flux = x)), nrow = 2)
}

# Update flux time series precheck plots according to prelim3 QC flags
if (interactive_session) {
  for (i in c("Tau", "H", "LE", "NEE")) {
    pdf(file.path(
      paths$Precheck,
      paste0(siteyear, "_", i, "_precheck_", Tstamp, ".pdf")), 
        width = 11.00, height = 8.27)
    plot_eddy(cbind(data, pre3_res), i, 
              paste0("qc_", i, "_prelim3"), paste0("qc_", i, "_prelim3"))
    dev.off()
  }
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
                        x = c("Tau", "H", "LE", "NEE"),
                        y = c("PAR", "Rn", "Rn", "PAR"),
                        z = c("wind_speed", "LE" , "H", "Tair")
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
mnames <- paste0("qc_", c("Tau", "H", "LE", "NEE"), "_man")
names(mnames) <- c("Tau", "H", "LE", "NEE")
is.na(mnames) <- !(mnames %in% names(man))
mnames <- as.list(mnames)

# Include manual QC among QC filter names to combine
prelim3_ad <- tribble(
  ~Tau,       ~H,       ~LE,       ~NEE,
  mnames$Tau, mnames$H, mnames$LE, mnames$NEE
)

prelim4 <- rbind(prelim3, prelim3_ad)

# Combine specified flags for given flux to produce preliminary flags
# - informative naming is useful during possible additional manual QC
pre4_res <- sapply(names(prelim4), 
                   function(x) combn_QC(data, na.omit(pull(prelim4, x))))
pre4_res <- as.data.frame(pre4_res)
names(pre4_res) <- paste0("qc_", names(pre4_res), "_prelim4")

# Update flux time series precheck plots according to prelim4 QC flags
if (interactive_session) {
  for (i in c("Tau", "H", "LE", "NEE")) {
    pdf(file.path(
      paths$Precheck,
      paste0(siteyear, "_", i, "_precheck_", Tstamp, ".pdf")), 
      width = 11.00, height = 8.27)
    plot_eddy(cbind(data, pre4_res), i, 
              paste0("qc_", i, "_prelim4"), paste0("qc_", i, "_prelim4"))
    dev.off()
  }
}

# If additional round of manual checking is not needed, skip to next chapter
man <- check_manually(cbind(data, pre4_res), paths$Quality_checking, 
                      vars = data.frame(
                        x = c("Tau", "H", "LE", "NEE"),
                        y = c("PAR", "Rn", "Rn", "PAR"),
                        z = c("wind_speed", "LE" , "H", "Tair")
                      ), 
                      qc_prefix = "qc_", qc_suffix = "_prelim4", 
                      interactive_session, siteyear)[-1]
summary_QC(man, names(man))

# Save the results to the main data frame
data[names(man)] <- man 

### Combine QC flags to final flag for gap-filling =============================

# Consider only existing manual flags
# - there might be variables without manual QC or 'man' can be NULL
mnames <- paste0("qc_", c("Tau", "H", "LE", "NEE"), "_man")
names(mnames) <- c("Tau", "H", "LE", "NEE")
is.na(mnames) <- !(mnames %in% names(man))
mnames <- as.list(mnames)

# Include manual QC among QC filter names to combine
prelim3_ad <- tribble(
  ~Tau,       ~H,       ~LE,       ~NEE,
  mnames$Tau, mnames$H, mnames$LE, mnames$NEE
  )

forGF <- rbind(prelim3, prelim3_ad)

# Combine specified flags for given flux to produce final forGF flags
forGF_res <- sapply(names(forGF), 
                   function(x) combn_QC(data, na.omit(pull(forGF, x))))
forGF_res <- as.data.frame(forGF_res)
names(forGF_res) <- paste0("qc_", names(forGF_res), "_forGF")

# Save the results to the main data frame
data[names(forGF_res)] <- forGF_res 

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
plot_QC_ind <- lapply(names(forGF), function(x) 
  summary_QC(data, na.omit(pull(forGF, x)), plot = TRUE, flux = x))
plot_QC_cum <- lapply(names(forGF), function(x) 
  summary_QC(data, na.omit(pull(forGF, x)), plot = TRUE, cumul = TRUE, flux = x))

ggsave(file.path(
  paths$QC_summary,
  paste0(siteyear, "_QC_summary_", Tstamp, ".png")),
  gridExtra::grid.arrange(grobs = plot_QC_ind, nrow = 2),
  type = "cairo-png", width = 297, height = 210, units = "mm")
ggsave(file.path(
  paths$QC_summary,
  paste0(siteyear, "_QC_summary_cumulative_", Tstamp, ".png")),
  gridExtra::grid.arrange(grobs = plot_QC_cum, nrow = 2),
  type = "cairo-png", width = 297, height = 210, units = "mm")

### Plot the quality checked data ==============================================

# Select the desired QC flags for plotting of respective fluxes
QC_plt <- names(forGF_res)
names(QC_plt) <- names(forGF)

# Save the plots that show the data with QC flag used for gap-filling
for (i in names(QC_plt)) {
  pdf(file.path(
    paths$Quality_checking,
    paste0(siteyear, "_forGF_QC_", i, "_", Tstamp, ".pdf")),
    width = 11.00, height = 8.27)
  plot_eddy(data, i, QC_plt[[i]], QC_plt[[i]])
  dev.off()
}

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
           
# Save essential output
essentials <- c(
  "timestamp", "GR", "qc_GR", "PAR", "qc_PAR", "Rn", "qc_Rn", "Tair",
  "qc_Tair", "Tsoil", "qc_Tsoil", "RH", "qc_RH", "VPD", "qc_VPD", "SWC",
  "qc_SWC", "P", "qc_P", "G", "qc_G", "Tau", "qc_Tau_forGF", "qc_Tau_SSITC",
  "rand_err_Tau", "H", "qc_H_forGF", "qc_H_SSITC", "rand_err_H", "LE",
  "qc_LE_forGF", "qc_LE_SSITC", "rand_err_LE", "NEE", "qc_NEE_forGF",
  "qc_NEE_SSITC", "rand_err_NEE", "H_strg", "LE_strg", "co2_strg", "wind_speed",
  "wind_dir", "ustar", "L", "zeta", "model", "x_peak", "x_70perc")

# Show which names are not available in the object save_data
# - use only those that are available
essentials <- choose_avail(names(save_data), essentials)

# Save the most important variables to later combine with gap-filling results
write_eddy(save_data[essentials],
           file.path(
             paths$Input_for_GF,
             paste0(siteyear, "_forGF_QC_essentials_", Tstamp, ".csv")))

# Save documentation information about executed quality control 
names(forGF) <- names(forGF_res)
writeLines(c(paste0(Tstamp, ":"),
             paste0("Quality controlled by ", name, " (", mail, ")"),
             "",
             paste0("Storage corrected fluxes: ", strg_applied),
             "",
             "Applied quality control scheme:", 
             capture.output(as.data.frame(forGF)),
             "",
             "Information about the R session:",
             capture.output(sessionInfo())), 
           file.path(
             paths$Quality_checking,
             paste0(siteyear, '_QC_info_', Tstamp, '.txt')), 
           sep = "\n")

#EOF
