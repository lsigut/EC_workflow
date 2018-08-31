### Description ================================================================

#' This code primarily aims at data quality checking (QC) and preparation of
#' input data for gap-filling. EddyPro and Meteo data inputs are merged and
#' reformatted. Data before and after QC are plotted and further statistics are
#' computed. Storage correction and correction of rotated vertical wind speed
#' (w_rot) is optional and overwrites the original values of respective
#' variables. Code developed by Ladislav Šigut (sigut.l@czechglobe.cz).

### Loading the required packages ==============================================

# You might need to install the packages first. In that case use:
# install.packages("packagename")
library("openeddy")
library("ggplot2")
library("gridExtra")
library("reshape2")

### Setting common paths, arguments and loading fetch filter boundaries ========

# Firstly manually set the desired working directory in RStudio interface
# (typically it is the folder specifying the siteyear, e.g. BKF16)
# This assumes that the subdirectory structure is already present
# (i.e. the !Site_structure_template from the server)

# Required inputs are post-processed data, meteo data 
# (./Level 1/Post-processing/EddyPro Output/) and
# Fetch_filter_boundaries.csv file (./Level 2/Quality checking)
pkg_version <- packageVersion("openeddy")
siteyear <- "KRP16" # Edit the site-year (included in folder and file names)
Tstamp <- format(Sys.time(), "%Y-%m-%d") # Timestamp of the computation
# (will be included in file names)

# Choose whether storage correction based on discrete (one point) approach 
# should be applied 
# (recommended for sites with short canopy, e.g. BKG, KRP, TRE)
apply_storage <- TRUE

# Output path for quality checking (Level 2 data):
(path <- "./Level 2/Quality checking/")

# Set precheck path (it will contain plots useful before quality checking)
(path_precheck <- paste(path, "Precheck/", sep = ""))

# Set wind direction dependency path
(path_wd_dep <- paste(path_precheck, "WD_dependency/", sep = ""))

# Set path for saving results as a gap-filling input 
GF_input <- "./Level 2/Input for gap-filling/"

# Set QC_summary path
(path_QC <- paste0(path, "QC_summary/"))

# Set current version of fetch filter boundaries
boundary_version <- "Fetch_filter_boundaries_20160206.csv"

# Edit the site abbreviation to load the correct boundary for given site
site_abr <- "KRP"

# Load the boundary used for fetch filter
boundary <- read.csv(paste(path, boundary_version, sep = ""))[, site_abr]

### Load EddyPro input and format it ===========================================

# Edit the filename:
EddyPro_input <- "./Level 1/Post-processing/EddyPro Output/eddypro_KrP_full_output_2017-05-20T180735_adv.csv"

# Load the EddyPro file
# -set correct skip parameter (number of lines above header in the input)
# -set correct file encoding (e.g. fileEncoding = "UTF-8") 
EP <- read_eddy(EddyPro_input, skip = 1, fileEncoding = "UTF-8")
str(EP) # check if loaded properly

# Load the original column names and check duplicates 
# -set correct skip parameter (number of lines above header in the input)
names <- read.csv(EddyPro_input, header = FALSE, nrows = 1, 
                  colClasses = "character", skip = 1)
names <- as.character(names) # compare with names(EP)
str(names) # check if loaded properly

# Correct column names
names <- gsub("co2_flux", "NEE", names) # assumption: co2_flux = NEE
names[names == "(z-d)/L"] <- "zeta"
names <- gsub("\\*", "star", names) # ustar, Tstar
names <- gsub("\\%", "perc", names) # signal contribution percentages
(names <- gsub("\\-|\\/", "_", names)) # dash and slash changed to underscore

# Update object EP
names(EP) <- names
varnames(EP) <- names
str(EP) # Are the names and attribute varnames edited correctly?

# Units corrections
(units <- units(EP, names = TRUE)) # Check the units
(units <- gsub(c("\\[|\\]"), "", units)) # Remove the square brackets
# -check whether µ was correctly interpreted due to the encoding issues

# Update object EP
units(EP) <- units
units(EP, names = TRUE) # Are units edited correctly?

# Remove RH and VPD from EP to avoid duplication 
# - eddy covariance estimates are less reliable than those of slow sensors
EP[c("RH", "VPD")] <- NULL

### Combine EddyPro timestamp from date and time ===============================

# Check the loaded input
head(EP[, 1:10]) # Starting at the beginning of year? (e.g. 1.1.2016 0:30)
tail(EP[, 1:10]) # Ending correctly? (e.g. 1.1.2017 0:00)

# Create timestamp, convert it to POSIXct format and center it
# Note: use GMT which just means "don’t modify times and dates – just use them
# exactly as they are irrespective of timezone information of your computer 
# system", i.e. use as is
timestamp <- paste(EP$date, EP$time)
timestamp <- strptime_eddy(timestamp, "%Y-%m-%d %H:%M", shift.by = -900)

# Combine timestamp with EP
EP <- cbind(timestamp, EP)
class(EP$timestamp) # should be POSIXt class now
head(EP$timestamp) # check if converted and centered correctly

### Load and format Meteo data =================================================

# Edit the filename
Meteo_input <- "./Level 1/Post-processing/EddyPro Output/Meteo_KRP16.csv"

# Load the Meteo file 
M <- read_eddy(Meteo_input, skip = 10, sep = ";")
str(M) # check if units loaded properly, names are incorrect

# Correct column names
M_names <- read.csv2(Meteo_input, header = FALSE, nrows = 1, 
                    colClasses = "character", skip = 9)
(M_names <- as.character(M_names)) # compare with names(M)
M_names_orig <- M_names # keep for documentation
(M_names <- gsub("\\/", "_", M_names)) # slash changed to underscore

# - link variable names from Meteo database with those required by scripts
DTB_varnames <- data.frame(timestamp = "date", P = "sumP", GR = "GR", Rn = "Rn", 
                           PAR = "PAR", Tair = "Ta", Tsoil = "Ts", RH = "RH", 
                           VPD = "VPD", stringsAsFactors = FALSE)
qc_vars <- grep("qcode", M_names)
for (i in names(DTB_varnames)) {
  M_names[grep(DTB_varnames[i], M_names)] <- i
}

# - remove other columns than those required
M_names_filter <- grep(paste(names(DTB_varnames), collapse = "|"), M_names)
M_names <- M_names[M_names_filter] 
M_names_orig <- M_names_orig[M_names_filter] 
M <- M[M_names_filter]
(M_names[qc_vars] <- paste0("qc_", M_names[qc_vars]))

# Update Meteo (var)names
names(M) <- M_names
varnames(M) <- M_names

# Units corrections
M_units <- units(M)
M_units <- gsub("\\*", "\\+1", M_units)
M_units <- gsub("st. ", "deg", M_units)

# Update Meteo units
units(M) <- M_units

# Convert timestamp to POSIXct format and center it
head(M$timestamp) # Starting at the beginning of year? (e.g. 1.1.2016 0:30:00)
tail(M$timestamp) # Ending correctly? (e.g. 1.1.2017 0:00:00)
M$timestamp <- strptime_eddy(M$timestamp, "%d.%m.%Y %H:%M", shift.by = -900)

# Merge Meteo and EddyPro inputs and save it with documentation ================

site <- merge(M, EP, by = "timestamp", all.x = TRUE)
units(site) <- c(units(M), units(EP[-1]))
varnames(site) <- c(varnames(M), varnames(EP[-1]))
save_site <- site

# Correct timestamp shift to its original state for saving
save_site$timestamp <- save_site$timestamp + 900 

# Change the timestamp formatting to its original state for saving
save_site$timestamp <- format(save_site$timestamp, format = "%Y-%m-%d %H:%M", 
                              tz = "GMT")

# Save the merged Meteo and EddyPro inputs with documentation
write_eddy(save_site, gsub("\\.csv", "_met.csv", EddyPro_input))
rm(save_site)
EP_doc <- readLines(gsub("\\.csv", ".txt", EddyPro_input), warn = FALSE)
docu <- c(paste0(Tstamp, ":"), 
          paste("Merged", basename(EddyPro_input), "and", basename(Meteo_input)),
          "Duplicated EddyPro columns RH and VPD were removed",
          "",
          "Variables from meteo database remapped to:",
          paste(M_names_orig, M_names, sep = " = ", collapse = "\n"), 
          "")
writeLines(c(docu, EP_doc), 
           gsub("\\.csv", "_met.txt", EddyPro_input), sep = "\n")

### Visual precheck ============================================================

# Display variables that can help identify problems with instruments
names <- c("u_rot", "v_rot", "w_unrot", "w_rot", "sonic_temperature", 
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

# Show which names are not available in the object site
names[!names %in% names(site)]
names <- names[names %in% names(site)]

# Edit ylim_list argument for plotted variables if needed
plot_precheck <- function(x, name, ylim_list = list(
  Tau = c(0, 2), H = c(-200, 800), LE = c(-200, 800), NEE = c(-100, 50),
  LE_strg = c(-5, 5), x_peak = c(0, 100), x_70perc = c(0, 500), 
  ustar = c(0, 1))) {
  ylim <- ylim_list[[name]]
  plot(x$timestamp, x[, name], ylim = ylim, 
       pch = 20, cex = 0.5, col = adjustcolor('black', 0.5),
       xaxt = "n", xlab = "timestamp", ylab = name, main = name)
  r <- as.POSIXct(round(range(x$timestamp), "months"))
  axis.POSIXct(1, at = seq(r[1], r[2], by = "month"), format = "%b-%y")
}

# save plots to single pdf
pdf(paste(path_precheck, siteyear, "_precheck_plots_", Tstamp, ".pdf",
          sep = ""), width = 11.00, height = 8.27)
invisible(lapply(names, plot_precheck, x = site))
dev.off()

# Show wind direction dependence of w_unrot with additional statistics
split_WDs <- function(x, var) {
  l <- split(x[, var], cut(x$wind_dir, seq(0, 360, 10)))
  data.frame(wind_dir = seq(5, 356, 10),
             mean = sapply(l, mean, na.rm = TRUE),
             SD = sapply(l, sd, na.rm = TRUE))
}
plot_WD_stats <- function(site, w_name, w_df) {
  y_span <- diff(range(site[, w_name], na.rm = TRUE))
  mean_val <- round(mean(site[, w_name], na.rm = TRUE), 3)
  ggplot(site, aes_string("wind_dir", w_name)) +
    geom_point(size = 1, na.rm = TRUE, alpha = 0.5) +
    geom_ribbon(data = w_df, 
                aes(x = wind_dir, ymin = mean - SD, ymax = mean + SD), 
                fill = "steelblue", alpha = 1/2, inherit.aes = FALSE) +
    geom_point(data = w_df, aes(wind_dir, mean), color = 'red') +
    geom_line(data = w_df, aes(wind_dir, mean), color = 'red') +
    annotate("text", size = 6,
             x = 250, 
             y = max(site[, w_name], na.rm = TRUE) - 0.1 * y_span, 
             label = paste0("mean ", w_name, ": ", mean_val))
}

w_unrot_df <- split_WDs(site, "w_unrot")
png(paste0(path_wd_dep, siteyear, "_w_unrot_WD_stats_", Tstamp, ".png"),
    width = 3508, height = 2480, res = 300)
print(plot_WD_stats(site, "w_unrot", w_unrot_df))
dev.off()

# Show wind direction dependence of w_rot with additional statistics
w_rot_df <- split_WDs(site, "w_rot")
png(paste0(path_wd_dep, siteyear, "_orig_w_rot_WD_stats_", Tstamp, ".png"),
    width = 3508, height = 2480, res = 300)
print(plot_WD_stats(site, "w_rot", w_rot_df))
dev.off()

# NB: w_rot has to be corrected by B0 planar fit coefficient if mean w_rot
# is not ~ 0 m s-1

# B0 applied for correction from the folder:
# none
# file: none

# Edit w_rot_correction. If w_rot does not require correction, set "none".
w_rot_correction <- "none" # if 2D rotated, set "none"

# Skip if w_rot_correction = "none"
if (!w_rot_correction == "none") {
  site$w_rot <- site$w_rot - w_rot_correction 
  # Confirm the succesful w_rot correction (skip if w_rot_correction = "none")
  w_rot_df <- split_WDs(site, "w_rot")
  png(paste0(path_wd_dep, siteyear, "_corrected_w_rot_WD_stats_", Tstamp, ".png"),
      width = 3508, height = 2480, res = 300)
  print(plot_WD_stats(site, "w_rot", w_rot_df))
  dev.off()
}

# WD dependency for individual months
wd_dep_names <- c("w_rot", "w_unrot", "ustar")
site$month <- ordered(month.name[as.POSIXlt(site$timestamp)$mon + 1], month.name)
q <- ggplot(site, aes_string("wind_dir", "w_rot")) +
  geom_point(size = 1, na.rm = TRUE, alpha = 0.25) +
  facet_wrap(~ month, scales = "free_x")

pdf(paste(path_wd_dep, siteyear, "_WD_dependency_", Tstamp, ".pdf",
          sep = ""), width = 11.00, height = 8.27)
invisible(lapply(wd_dep_names, function(x) print(q + aes_string(y = x))))
dev.off()

# Create the precheck plot (PC_plt) vector for simplified plotting:
(PC_plt <- c(Tau = "qc_Tau",
             H   = "qc_H",
             LE  = "qc_LE",
             NEE = "qc_NEE"))

# Save the plots in Precheck folder
for (i in names(PC_plt)) {
  pdf(paste(path_precheck, siteyear, "_", i, "_precheck_", Tstamp, ".pdf",
            sep = ""), width = 11.00, height = 8.27)
  plot_eddy(site, i, PC_plt[[i]], PC_plt[[i]])
  dev.off()
}

### Quality checking - preliminary flags =======================================

# Extract the standardized tests/filters
QC <- extract_QC(site) 
str(QC) # Check the results
lapply(QC, table, useNA = "always")
site[names(QC)] <- QC # Save the results to the main data frame

# Create data frame with test names and specification of test additivity
fluxes <- c("Tau", "H", "LE", "NEE")
qc_names <- c("qc_flux", "qc_instrum_abslim", "qc_instrum_spikesHF", 
              "qc_flux_missfrac", "qc_flux_scf", "qc_ALL_wresid")
additive <- c(rep(FALSE, 5), TRUE) # Only the qc_ALL_wresid filter is additive
na.as <- NA
pattern <- data.frame(qc_names = qc_names, additive = additive, na.as = na.as)

# Function that renames tests according to used flux
assign_tests <- function(x, flux) {
  x$qc_names <- gsub("flux", flux, x$qc_names)
  x$qc_names <- 
    if (flux %in% c("Tau", "H")) gsub("instrum", "SA", x$qc_names) else
      gsub("instrum", "SA_IRGA", x$qc_names)
  return(x)
}

# Retrieve a list of data frames with test specifications for each flux
prelim <- lapply(fluxes, assign_tests, x = pattern)
names(prelim) <- fluxes
prelim # check the list

# Function to simplify combining flags based on a list of test names
simplify_combn_QC <- function(x, list, flux, name_suffix) {
  combn_QC(x, list[[flux]]$qc_names, paste("qc", flux, name_suffix, sep = "_"), 
           list[[flux]]$additive, list[[flux]]$na.as)
}

# Combine existing tests/filters to produce preliminary flags:
H_prelim <- simplify_combn_QC(site, prelim, "H", "prelim")
LE_prelim <- simplify_combn_QC(site, prelim, "LE", "prelim")
NEE_prelim <- simplify_combn_QC(site, prelim, "NEE", "prelim")

### Quality checking - flux interdependency and composite flags ================

# H_prelim is needed in interdep() only if IRGA = "open", otherwise not used
interdep_enclosed <- interdep(LE_prelim) # Interdependency for enclosed IRGA
str(interdep_enclosed)
lapply(interdep_enclosed, table, useNA = "always")
site[names(interdep_enclosed)] <- interdep_enclosed # Save the results

# Create data frame with test names and specification of test additivity
qc_names <- c("qc_flux", "qc_instrum_abslim", "qc_instrum_spikesHF", 
              "qc_flux_missfrac", "qc_flux_scf", "qc_ALL_wresid", 
              "qc_flux_interdep")
additive <- c(rep(FALSE, 5), rep(TRUE, 2)) # _wresid & _interdep are additive
na.as <- NA
pattern <- data.frame(qc_names = qc_names, additive = additive, na.as = na.as)

# Retrieve a list of data frames with test specifications for each flux
comp <- lapply(fluxes, assign_tests, x = pattern)
names(comp) <- fluxes

# Remove not applied tests and check the list
comp$Tau <- comp$Tau[-7, ] # qc_Tau_interdep not used
comp$LE <- comp$LE[-7, ] # !!! With open-path analyzer (IRGA = "open") comp$LE
                         # would include also "qc_LE_interdep"
comp # check the list

# Combine existing tests/filters to produce composite flags:
Tau_comp <- simplify_combn_QC(site, comp, "Tau", "composite")
H_comp <- simplify_combn_QC(site, comp, "H", "composite")
LE_comp <- simplify_combn_QC(site, comp, "LE", "composite")
NEE_comp <- simplify_combn_QC(site, comp, "NEE", "composite")

# Show effect of flux interdependency test on QC flags:
# Interpretation:
# 0 shows how many flag 0 data were corrected to be flag 1,
# 1 shows overall change in flag 1 data
# 2 shows how many flag 1 data were corrected to be flag 2,
abs(table(H_comp) - table(H_prelim))
abs(table(LE_comp) - table(LE_prelim)) # No difference for (en)closed IRGA
abs(table(NEE_comp) - table(NEE_prelim))

# Save the computed composite flags:
site$qc_Tau_composite <- Tau_comp
site$qc_H_composite <- H_comp
site$qc_LE_composite <- LE_comp
site$qc_NEE_composite <- NEE_comp

### Storage correction of fluxes ===============================================

if (apply_storage) {
  # Storage flux is added to the respective original flux
  # - correction  overwrites the original values of respective variables 
  site$H <- add_st(site$H, "H", site$H_strg)
  site$LE <- add_st(site$LE, "LE", site$LE_strg)
  site$NEE <- add_st(site$NEE, "NEE", site$co2_strg)
}

### Low frequency flux despiking ===============================================

# Not applied to Tau
plot(site$H, ylim = c(-200, 400))
desp_H <- despikeLF(site, var = "H", qc_flag = "qc_H_composite",
                    name_out = "qc_H_spikesLF", var_thr = c(-200, 800))
plot(site$LE, ylim = c(-200, 400))
desp_LE <- despikeLF(site, var = "LE", qc_flag = "qc_LE_composite",
                     name_out = "qc_LE_spikesLF", var_thr = c(-200, 800))
plot(site$NEE, ylim = c(-50, 50))
desp_NEE <- despikeLF(site, var = "NEE", qc_flag = "qc_NEE_composite",
                      "qc_NEE_spikesLF", var_thr = c(-100, 100))

# Check the results
# The test has 3 outcomes: 0 - OK; 2 - spike; NA - excluded from despiking
# NB: Since NA means not checked, NA should be thus interpreted as flag 0
lapply(list(desp_H = desp_H, desp_LE = desp_LE, desp_NEE = desp_NEE), 
       table, useNA = "a")

# Save the results
site$qc_H_spikesLF   <- desp_H
site$qc_LE_spikesLF  <- desp_LE
site$qc_NEE_spikesLF <- desp_NEE

### Filtering fluxes with outlying random error ================================

# plot(site$rand_err_H)
# err_H <- despikeLF(site, var = "rand_err_H", qc_flag = "qc_H_composite",
#                    name_out = "qc_H_rand_err_spikesLF", var_thr = c(0, 50))
# plot(site$rand_err_LE)
# err_LE <- despikeLF(site, var = "rand_err_LE", qc_flag = "qc_LE_composite",
#                     name_out = "qc_LE_rand_err_spikesLF", var_thr = c(0, 50))
# plot(site$rand_err_NEE)
# err_NEE <- despikeLF(site, var = "rand_err_NEE", 
#                       qc_flag = "qc_NEE_composite", 
#                       name_out = "qc_NEE_rand_err_spikesLF", 
#                       var_thr = c(0, 20))
# 
# # Check the results
# lapply(list(err_H = err_H, err_LE = err_LE, err_NEE = err_NEE), 
#        table, useNA = "a")
# 
# # Save the results
# site$qc_H_rand_err_spikesLF    <- err_H
# site$qc_LE_rand_err_spikesLF   <- err_LE
# site$qc_NEE_rand_err_spikesLF <- err_NEE

### Filtering fluxes with outlying variance ====================================

# plot(site$ts_var, ylim = c(0, 10))
# var_H <- despikeLF(site, var = "ts_var", qc_flag = "qc_H_composite",
#                    name_out = "qc_H_ts_var_spikesLF", var_thr = c(0, 5))
# plot(site$h2o_var, ylim = c(0, 10))
# var_LE <- despikeLF(site, var = "h2o_var", qc_flag = "qc_LE_composite",
#                     name_out = "qc_LE_h2o_var_spikesLF", var_thr = c(0, 5))
# plot(site$co2_var, ylim = c(0, 200))
# var_NEE <- despikeLF(site, var = "co2_var", qc_flag = "qc_NEE_composite", 
#                       name_out = "qc_NEE_co2_var_spikesLF", 
#                       var_thr = c(0, 200))
# 
# # Check the results
# lapply(list(var_H = var_H, var_LE = var_LE, var_NEE = var_NEE), 
#        table, useNA = "a")
# 
# # Save the results
# site$qc_H_ts_var_spikesLF  <- var_H
# site$qc_LE_h2o_var_spikesLF <- var_LE
# site$qc_NEE_co2_var_spikesLF <- var_NEE

### Fetch filter ===============================================================

# NB: qc_ALL_fetch70 is actually not applied to Tau
FF <- fetch_filter(site, "x_70perc", "wind_dir", boundary, "qc_ALL_fetch70")
table(FF, useNA = "always", dnn = "qc_ALL_fetch70")
site$qc_ALL_fetch70 <- FF

### Combine QC flags used for gap-filling ======================================

# Create a list for each flux with additional tests and their properties
# add <- data.frame(H = c("qc_H_spikesLF", "qc_H_rand_err_spikesLF", 
#                             "qc_H_ts_var_spikesLF", "qc_ALL_fetch70"),
#                   LE = c("qc_LE_spikesLF", "qc_LE_rand_err_spikesLF", 
#                              "qc_LE_h2o_var_spikesLF", "qc_ALL_fetch70"),
#                   NEE = c("qc_NEE_spikesLF", "qc_NEE_rand_err_spikesLF", 
#                           "qc_NEE_co2_var_spikesLF", "qc_ALL_fetch70"),
#                   stringsAsFactors = FALSE)
# properties <- data.frame(additive = FALSE, na.as = c(rep(0, 3), NA))
####### add and properties above include rand_err and var despiking
####### if not required rewrite add and properties below
add <- data.frame(H = c("qc_H_spikesLF", "qc_ALL_fetch70"),
                  LE = c("qc_LE_spikesLF", "qc_ALL_fetch70"),
                  NEE = c("qc_NEE_spikesLF", "qc_ALL_fetch70"),
                  stringsAsFactors = FALSE)
properties <- data.frame(additive = FALSE, na.as = c(0, NA))

add_list <- lapply(add, function(qc_names) cbind(qc_names, properties))

# Create a complete list of tests within forGF scheme (does not apply for Tau)
forGF <- lapply(1:3, function(x) rbind(comp[-1][[x]], add_list[[x]]))
names(forGF) <- c("H", "LE", "NEE")
forGF # Check the list

# Combine existing tests/filters to produce forGF flags:
H_forGF <- simplify_combn_QC(site, forGF, "H", "forGF")
LE_forGF <- simplify_combn_QC(site, forGF, "LE", "forGF")
NEE_forGF <- simplify_combn_QC(site, forGF, "NEE", "forGF")

# Show the difference between composite and forGF flags:
# Interpretation:
# 0 shows how many flag 0 data were corrected to be flag 1,
# 1 shows overall change in flag 1 data
# 2 shows how many flag 1 data were corrected to be flag 2,
abs(table(H_forGF) - table(H_comp))
abs(table(LE_forGF) - table(LE_comp))
abs(table(NEE_forGF) - table(NEE_comp))

# Save the computed forGF flags:
site$qc_H_forGF  <- H_forGF
site$qc_LE_forGF <- LE_forGF
site$qc_NEE_forGF    <- NEE_forGF

### QC summary and the graphical display of the results ========================

# Function to simplify combining flags based on a list of test names
simplify_summary_QC <- function(x, df, cumul = FALSE, plot = FALSE, 
                                perc = TRUE, flux = NULL) {
  with(df, summary_QC(x, qc_names, na.as, cumul, additive, plot, perc, flux))
}

# Show the percantage of fluxes flagged by individual tests/filter and
# their cumulative effect
list_QC <- list(Tau = simplify_summary_QC(site, comp$Tau),
                H = simplify_summary_QC(site, forGF$H),
                LE = simplify_summary_QC(site, forGF$LE),
                NEE = simplify_summary_QC(site, forGF$NEE),
                
                Tau_cumul = simplify_summary_QC(site, comp$Tau, TRUE),
                H_cumul = simplify_summary_QC(site, forGF$H, TRUE),
                LE_cumul = simplify_summary_QC(site, forGF$LE, TRUE),
                NEE_cumul = simplify_summary_QC(site, forGF$NEE, TRUE))
                
# Save each element of list_QC to separate file in path_QC
for (i in names(list_QC)) {
  write.csv(list_QC[[i]], paste(path_QC, siteyear, "_QC_summary_", i, "_",
                                Tstamp, ".csv", sep = ""))
}

# Save the QC summary results as plots
pdf(paste(path_QC, siteyear, "_QC_summary_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
gridExtra::grid.arrange(
  simplify_summary_QC(site, comp$Tau, plot = TRUE, flux = "Tau"),
  simplify_summary_QC(site, forGF$H, plot = TRUE, flux = "H"),
  simplify_summary_QC(site, forGF$LE, plot = TRUE, flux = "LE"),
  simplify_summary_QC(site, forGF$NEE, plot = TRUE, flux = "NEE"),
  nrow = 2)
gridExtra::grid.arrange(
  simplify_summary_QC(site, comp$Tau, TRUE, TRUE, flux = "Tau"),
  simplify_summary_QC(site, forGF$H, TRUE, TRUE, flux = "H"),
  simplify_summary_QC(site, forGF$LE, TRUE, TRUE, flux = "LE"),
  simplify_summary_QC(site, forGF$NEE, TRUE, TRUE, flux = "NEE"),
  nrow = 2)
dev.off()

### Quality checked data plotting ==============================================

# Select the desired QC flag for plotting
QC_plt <- c(Tau    = "qc_Tau_composite",
            H  = "qc_H_forGF",
            LE = "qc_LE_forGF",
            NEE    = "qc_NEE_forGF")

# Save the plots that show the data with QC flag used for gap-filling
for (i in names(QC_plt)) {
  pdf(paste(path, siteyear, "_forGF_QC_", i, "_", Tstamp, ".pdf", sep = ""),
      width = 11.00, height = 8.27)
  plot_eddy(site, i, QC_plt[[i]], QC_plt[[i]])
  dev.off()
}

### Write QC results to a file =================================================

# Correct timestamp to its original state
site$timestamp <- site$timestamp + 900 

# Remap variable names and filter columns needed for online gap-filling tool:
OT_in <- set_OT_input(site,
                      c("qc_NEE_forGF", "NEE", "qc_LE_forGF", "LE",
                        "qc_H_forGF", "H", "GR", "Tair", "Tsoil", "RH",
                        "VPD", "qc_Tau_composite", "ustar"))

# Save the standardized Online Tool/REddyProc input to '.txt' file:
write_eddy(OT_in, paste(GF_input, siteyear, "_", Tstamp, ".txt", sep = ""),
           quote = FALSE, sep = "\t")

# Change the timestamp formatting to its original state:
site$timestamp <- format(site$timestamp, format = "%Y-%m-%d %H:%M", tz = "GMT")

# Save the full output:
write_eddy(site, paste(path, siteyear, "_forGF_QC_full_output_", Tstamp, ".csv", 
                       sep = ""))

# Save essential output
essentials <- c("timestamp", "DOY", "date", "time", "P", "GR", "Rn", "PAR", 
                "Tair", "Tsoil", "RH", "VPD", "Tau", "qc_Tau", 
                "qc_Tau_composite", "rand_err_Tau", "w_var", "u_var", 
                "H", "qc_H", "qc_H_composite", "H_strg", "H", 
                "qc_H_forGF", "rand_err_H", "ts_var", "LE", "ET", 
                "h2o_flux", "rand_err_h2o_flux", "h2o_strg", "h2o_v_adv", 
                "qc_LE", "qc_LE_composite", "LE_strg", "LE", 
                "qc_LE_forGF", "rand_err_LE", "h2o_var", "NEE", 
                "qc_NEE", "qc_NEE_composite", "co2_strg", "NEE", 
                "qc_NEE_forGF", "rand_err_NEE", "co2_var", "co2_v_adv", 
                "wind_speed", "wind_dir", "ustar", "L", "zeta", "bowen_ratio")

# Show which names are not available in the object site
essentials[!essentials %in% names(site)]
# - use only those that are available
essentials <- essentials[essentials %in% names(site)]

write_eddy(site[essentials], paste(GF_input, siteyear, "_forGF_QC_essentials_", 
                                   Tstamp, ".csv", sep = ""))

# Save settings used for the computations
settings <- list(openeddy_pkg_version = pkg_version,
                 siteyear = siteyear, 
                 QC_input = EddyPro_input,
                 QC_output_path = path,
                 precheck_path = path_precheck,
                 QC_summary_path = path_QC,
                 gap_filling_input = GF_input,
                 wind_direction_dependency_path = path_wd_dep,
                 ROI_boundary_version = boundary_version,
                 boundary_used = site_abr,
                 w_rot_correction = w_rot_correction,
                 storage_correction = apply_storage)
sink(paste(path, siteyear, '_settings_', Tstamp, '.txt', sep = ""))
settings
sink()

# Save test names that would be needed to reconstruct composite and forGF flags
tests <- list(composite = comp,
              forGF = forGF)
sink(paste(path, siteyear, '_test_combining_', Tstamp, '.txt', sep = ""))
tests
sink()

#EOF
