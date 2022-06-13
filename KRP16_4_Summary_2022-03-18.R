### Description ================================================================

# Eddy covariance workflow part 4/4 (https://github.com/lsigut/EC_workflow)
# This code primarily aims to produce the summary of processed data at different
# timescales (daily, weekly, monthly, yearly) and to plot them.
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
attach_pkg("openair")

# Workflow is currently aligned only with specific package version
# - package version should be neither higher or lower
if (packageVersion("openeddy") != "0.0.0.9006")
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9006'")

### Provide metadata and set file paths and arguments ==========================

# Edit the siteyear
# - included in file names
siteyear <- "KRP16"

# Timestamp of the computation
# - automated, will be included in file names
Tstamp <- format(Sys.time(), "%Y-%m-%d") 

# Load the list of folder structure paths
# - automated, no input required if proposed folder structure is followed
paths <- structure_eddy()

# Input path for summary (automated)
# - gap-filled and partitioned data
path_in <- grep(paste0(siteyear, ".*GF_essentials.*csv"), 
                list.files(paths$Gap_filling, full.names = TRUE),
                value = TRUE)

# Specify the time shift (in seconds) to be applied to the date-time information
# in order to represent the center of averaging period
shift.by <- -900

### Load the input file and convert timestamp and date =========================

data <- read_eddy(path_in)
str(data) # check if loaded properly
data$timestamp <- strptime_eddy(data$timestamp, "%Y-%m-%d %H:%M", 
                                shift.by = shift.by)
head(data$timestamp)

### Plot half-hourly data before unit conversion ===============================

# Specify variables needed for different procedures later 
# - used for averaging, summation, uncertainty estimation and plotting
mean <- c("Tair", "Tsoil", "RH", "VPD", "GR", "Rn", "PAR", "H_f", "H_fqc", 
          "LE_f", "LE_fqc", "ET_f", "ET_fqc", "NEE_uStar_f", "NEE_uStar_fqc", 
          "GPP_uStar_f", "GPP_DT_uStar", "Reco_uStar", "Reco_DT_uStar")
mean <- choose_avail(names(data), mean)
sum <- c("P", "GR", "Rn", "PAR", "H_f", "LE_f", "ET_f", 
         grep("(NEE|GPP).+f$", names(data), value = TRUE),
         grep("GPP_DT.+[^D]$", names(data), value = TRUE),
         grep("Reco.+[^D]$", names(data), value = TRUE))
sum <- choose_avail(names(data), sum)
err_agg <- grep("(^H|^LE|^ET|^NEE|^Reco|^GPP).*(sd|SD)$", names(data), 
                value = TRUE)

# Specify variables for plots at half-hour resolution
hh_vars <- grep("[^c]$", unique(c(mean, sum, err_agg)), value = TRUE)

# Print half-hourly results to pdf and png
pdf(file.path(
  paths$Summary,
  paste0(siteyear, "_half-hourly_plots_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
invisible(lapply(hh_vars, plot_hh, x = data))
dev.off()

for (i in hh_vars) {
  png(file.path(
    paths$png,
    paste0(siteyear, "_hh_", i, "_", Tstamp, ".png")),
      width = 3508, height = 2480, res = 400)
  plot_hh(data, i)
  dev.off()
}

### Plot wind roses and 1D footprint results ===================================

# Format and prepare data for openair package
wrose_all <- data.frame(ws = data$wind_speed, wd = data$wind_dir)
wrose_all$time <- cut(data$PAR, c(-Inf, 10, Inf), 
                      labels = c("nighttime", "daytime"))
wrose_all$months <- ordered(month.name[as.POSIXlt(data$timestamp)$mon + 1], 
                            month.name)
# - number of groups to which zeta should be cut
ngroups <- 6
# - identify cut breakpoints
breakpoints <- quantile(data$zeta, seq(0, 1, len = ngroups + 1), na.rm = TRUE)
# - cut zeta so intervals are closed on the left and highest value is included
# - negative (positive) zeta represents unstable (stable) conditions
wrose_all$stability <- cut(data$zeta, breakpoints, right = FALSE,
                           include.lowest = TRUE)

# Print all to pdf
pdf(file.path(
  paths$Summary,
  paste0(siteyear, "_wind_roses_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd")]), ], 
         angle = 22.5, paddle = FALSE, breaks = 5)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "time")]), ], 
         type = "time", angle = 22.5, paddle = FALSE, breaks = 5)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "months")]), ], 
         type = "months", angle = 45, paddle = FALSE, breaks = 5, 
         grid.line = 10)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "stability")]), ], 
         type = "stability", angle = 22.5, paddle = FALSE, 
         breaks = 5, grid.line = 10, 
         main = "Zeta parameter based stability classes")
ggplot_stats(data, "wind_dir", "x_peak", circular = TRUE)
ggplot_stats(data, "wind_dir", "x_70perc", circular = TRUE)
ggplot_stats(data, "wind_dir", "wind_speed", circular = TRUE)
ggplot_stats(data, "wind_dir", "ustar", circular = TRUE)
ggplot_stats(data, "wind_dir", "zeta", circular = TRUE)
dev.off()

# Print separately to png
# - png_util() is defined in utilities.R as a helper function for saving plots
png_util("wind_rose_all")
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd")]), ], 
         angle = 22.5, paddle = FALSE, breaks = 5)
dev.off()

png_util("wind_rose_day-night")
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "time")]), ], 
         type = "time", angle = 22.5, paddle = FALSE, breaks = 5)
dev.off()

png_util("wind_rose_months")
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "months")]), ], 
         type = "months", angle = 45, paddle = FALSE, breaks = 5, 
         grid.line = 10)
dev.off()

png_util("wind_rose_stability")
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "stability")]), ], 
         type = "stability", angle = 22.5, paddle = FALSE, 
         breaks = 5, grid.line = 10, 
         main = "Zeta parameter based stability classes")
dev.off()

png_util("wind_dir_x_peak")
ggplot_stats(data, "wind_dir", "x_peak", circular = TRUE)
dev.off()

png_util("wind_dir_x_70perc")
ggplot_stats(data, "wind_dir", "x_70perc", circular = TRUE)
dev.off()

png_util("wind_dir_wind_speed")
ggplot_stats(data, "wind_dir", "wind_speed", circular = TRUE)
dev.off()

png_util("wind_dir_ustar")
ggplot_stats(data, "wind_dir", "ustar", circular = TRUE)
dev.off()

png_util("wind_dir_zeta")
ggplot_stats(data, "wind_dir", "zeta", circular = TRUE)
dev.off()

### Compute summaries for different intervals ==================================

# Daily, weekly, monthly and yearly means, sums and uncertainties
intervals <- c("%Y-%m-%d", "%W_%y", "%b-%y", "%Y")
agg_periods <- c("day-1", "week-1", "month-1", "year-1")

means <- lapply(intervals, function(x) agg_mean(data[c("timestamp", mean)], x))
sums <- mapply(function(x, y) agg_sum(data[c("timestamp", sum)], x, agg_per = y),
               x = intervals,
               y = agg_periods,
               SIMPLIFY = FALSE)
fsd <- mapply(function(x, y) agg_fsd(data, x, agg_per = y),
              x = intervals,
              y = agg_periods,
              SIMPLIFY = FALSE)
DT_SD <- mapply(function(x, y) agg_DT_SD(data, x, agg_per = y),
                x = intervals,
                y = agg_periods,
                SIMPLIFY = FALSE)

# Compute additional parameters
pars <- vector("list", length(means))

for (i in seq_along(means)) {
  pars[[i]] <- means[[i]][c("Intervals", "days")]
  pars[[i]]$bowen_ratio_f <- sums[[i]]$H_f_sum / sums[[i]]$LE_f_sum
  pars[[i]]$evaporative_fraction <- sums[[i]]$LE_f_sum / 
    (sums[[i]]$H_f_sum + sums[[i]]$LE_f_sum)
  pars[[i]]$closure_fraction <- (sums[[i]]$LE_f_sum + sums[[i]]$H_f_sum) / 
    sums[[i]]$Rn_sum
  openeddy::units(pars[[i]]) <- rep("-", ncol(pars[[i]]))
}

# Mark days with positive carbon uptake and Tair above 5 degC (CUP and GSL)
pars[[1]]$CUP <- 0L
pars[[1]]$CUP[sums[[1]]$NEP_uStar_f_sum > 0] <- 1L
pars[[1]]$GSL <- 0L
pars[[1]]$GSL[means[[1]]$Tair_mean > 5] <- 1L

# Count days with positive carbon uptake and Tair above 5 degC in intervals
# - requires daily resolution of factors representing aggregation intervals
# - number of days in interval is included
# - first element of list "days" is left NULL to simplify merging
resol <- c("%W_%y", "%b-%y", "%Y")
days <- vector("list", length(means))

for (i in seq_along(resol)) {
  grouping_str <- strftime(means[[1]]$Intervals, format = resol[i], tz = "GMT")
  grouping <- factor(grouping_str, levels = unique(grouping_str))
  days[[i+1]] <- aggregate(means[[1]]$Intervals, list(grouping), length)[2]
  days[[i+1]]$CUP <- aggregate(pars[[1]]$CUP, list(grouping), sum)[, 2]
  days[[i+1]]$GSL <- aggregate(pars[[1]]$GSL, list(grouping), sum)[, 2]
  names(days[[i+1]])[1] <- "days"
}

# Report percentage of original measured data in intervals
orig <- c("H_orig", "LE_orig", "NEE_uStar_orig")
avail <- vector("list", length(intervals))

for (i in seq_along(intervals)) {
  g_str <- strftime(data$timestamp, format = intervals[i], tz = "GMT")
  g <- factor(g_str, levels = unique(g_str))
  ilen <- aggregate(data$timestamp, list(g), length)
  av <- aggregate(data[orig], list(g), function(x) sum(!is.na(x)))
  avail[[i]] <- round(av[-1] / ilen$x * 100, 1)
  openeddy::units(avail[[i]]) <- rep("%", 3)
}

# Merge additional variables with pars data frame and add varnames and units
for (i in seq_along(means)) {
  pars[[i]][names(days[[i]][-1])] <- days[[i]][-1]
  pars[[i]][orig] <- avail[[i]]
  openeddy::varnames(pars[[i]][c("CUP", "GSL", orig)]) <- c(
    "Carbon Uptake Period", "Growing Season Length", 
    paste("percentage of", gsub("_orig", "", orig), "original records"))
  openeddy::units(pars[[i]][c("CUP", "GSL")]) <- rep(
    paste0("days ", agg_periods[i]), 2)
}


# Combine summaries to a single data frame per interval (all in one list)
# and round and save the results
summaries <- vector("list", length(means))
resol_names <- c("daily", "weekly", "monthly", "yearly")
names(summaries) <- resol_names

# Loop uses filenames() function saved in utilities.R
for (i in seq_along(means)) {
  summaries[[i]] <- cbind(means[[i]], fsd[[i]]$mean[-c(1:2)], 
                          DT_SD[[i]]$mean[-c(1:2)], sums[[i]][-c(1:2)], 
                          fsd[[i]]$sum[-c(1:2)], DT_SD[[i]]$sum[-c(1:2)],
                          pars[[i]][-c(1:2)])
  write_eddy(round_df(summaries[[i]]), filenames(resol_names[i]))
}

### Plot summaries for different intervals ===================================== 

vars <- names(summaries$daily)
vars <- vars[!(vars %in% c("Intervals", "days"))]

# Plot aggregated variables per given intervals to pdf 
for (interval in c("daily", "weekly", "monthly")) {
  pdf(file.path(
    paths$Summary,
    paste0(siteyear, "_", interval, "_plots_", Tstamp, ".pdf")),
      width = 11.00, height = 8.27)
  for (var in vars) {
    barplot_agg(summaries[[interval]], var = var, interval)
  }
  dev.off()
}

# Plot aggregated variables per given intervals to png 
for (interval in c("daily", "weekly", "monthly")) {
  for (var in vars) {
    png(file.path(
      paths$png,
      paste0(siteyear, "_", interval, "_", var, "_", Tstamp, ".png")),
      width = 3508, height = 2480, res = 350)
    barplot_agg(summaries[[interval]], var = var, interval)
    dev.off()
  }
}

### Compute Griebel et al. (2020) budgets ======================================

# Compute budgets with different consideration of space-time equity
G20 <- Griebel20_budgets(
  data, "timestamp", "wind_dir", "NEE_uStar_f", "NEE_uStar_fqc")

# Fix random number generator state before spti_boot() for reproducible results
set.seed(621)

# Estimate uncertainty of space-time-equitable budgets
spti_unc <- spti_boot(
  data, "timestamp", "wind_dir", "NEE_uStar_f", "NEE_uStar_fqc")

# Calculate spatio-temporal sampling coverage
spti_cov <- spti_coverage(
  data, "timestamp", "wind_dir", "NEE_uStar_f", "NEE_uStar_fqc")

# Combine, round and save results
# - traditional_budget and standardized_budget differ only if multiple years 
#   are included
G20_all <- round_df(cbind(G20[-6], spti_unc[-(1:2)], spti_cov[-1]))
write_eddy(G20_all, 
           file.path(
             paths$Summary,
             paste0(siteyear, "_Griebel_et_al_2020_budgets_", Tstamp, ".csv")))

### Plot spatio-temporal sampling coverage to pdf and png ======================

# Create a list of ggplot objects
spti_covp <- spti_coverage(
  data, "timestamp", "wind_dir", "NEE_uStar_f", "NEE_uStar_fqc",
  plot = TRUE)

# Save all ggplots to pdf
pdf(file.path(
  paths$Summary,
  paste0(siteyear, "_spatio-temporal_sampling_coverage_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
spti_covp
dev.off()

# Save specified ggplots to png
png_util("spatial_sampling_coverage")
spti_covp[[1]]$spatial_sampling_coverage
dev.off()

png_util("temporal_sampling_coverage")
spti_covp[[1]]$temporal_sampling_coverage
dev.off()

### Quantify and plot estimates of flux uncertainty ============================

# Select variables and compute errors in half-hourly resolution
# - only for visualization purposes (outputs not saved)
# - only original measured data after quality control are plotted
# - nighttime errors have different properties than daytime (not distinguished)
# - absolute relative random error (divided by actual flux value)
# - absolute relative error: based on gag-filling look-up table (divided by 
#   actual flux value); '_fsd' represents SD, thus fsd * 1.96 represents 
#   95% confidence interval estimate
vars <- c("H", "LE", "NEE")

# Open pdf device
pdf(file.path(
  paths$Summary,
  paste0(siteyear, "_flux_uncertainty_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
for (var in vars) {
  # - only a single match is expected for rel_err, var_fsd, rel_fsd
  rel_err <- paste0("abs_rel_rand_err_", var)
  var_fsd <- grep(paste0("^", var, ".*_fsd$"), names(data), value = TRUE)
  rel_fsd <- paste0("abs_rel_", var_fsd)
  
  data[rel_err] <- abs(data[paste0("rand_err_", var)]/data[var]) * 100
  data[rel_fsd] <- abs(data[var_fsd] * 1.96 /data[var]) * 100
  openeddy::units(data[c(rel_err, rel_fsd)]) <- c("%", "%")
  
  # find respective gap-filling quality flag '_fqc' to given var
  # - e.g. "H_fqc" for H, but "NEE_uStar_fqc" for NEE 
  var_fqc <- grep(paste0(var, ".*_fqc$"), names(data), value = TRUE)
  # choose only original data (not gap-filled)
  fqc_subset <- data[data[, var_fqc] == 0, ]
  # restore units (subsetting removes units)
  openeddy::units(fqc_subset) <- openeddy::units(data)
  
  # Plot to opend pdf device
  print(ggplot_stats(fqc_subset, "wind_dir", rel_err, 
                     circular = TRUE, qrange = c(0.005, 0.99)))
  print(ggplot_stats(fqc_subset, "wind_dir", rel_fsd, 
                     circular = TRUE, qrange = c(0.005, 0.99)))
  
  # Open png device, plot and close the device
  png_util(rel_err)
  print(ggplot_stats(fqc_subset, "wind_dir", rel_err, 
                     circular = TRUE, qrange = c(0.005, 0.99)))
  dev.off()
  
  # Open png device, plot and close the device
  png_util(rel_fsd)
  print(ggplot_stats(fqc_subset, "wind_dir", rel_fsd, 
                     circular = TRUE, qrange = c(0.005, 0.99)))
  dev.off()
}
# Close the pdf device with all plots from the loop
dev.off()

# EOF
