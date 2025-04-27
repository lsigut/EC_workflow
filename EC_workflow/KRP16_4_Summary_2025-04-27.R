### Description ================================================================

# Eddy covariance workflow part 4/4 (https://github.com/lsigut/EC_workflow)
# This code primarily aims to produce the summary of processed data at different
# timescales (daily, weekly, monthly, yearly) and to plot them.
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
source("utilities.R")

# Attach packages from GitHub
# - you might need to have RTools for Windows machine to install openeddy:
#   https://cran.r-project.org/bin/windows/Rtools
# - uses attach_pkg() function saved in utilities.R
attach_pkg("openeddy", github = "lsigut/openeddy")

# Attach packages from CRAN
attach_pkg("openair", "tibble")

# Check if openeddy version conforms to requirements
if (packageVersion("openeddy") < package_version("0.0.0.9009"))
  warning("this version of workflow works reliably only with openeddy version ",
          "'0.0.0.9009'")

### Provide metadata and set file paths and arguments ==========================

# Load the site-year settings file
settings_file <- list.files(pattern = "settings.R", full.names = TRUE)
source(settings_file)

# Timestamp of the computation
# - automated, will be included in file names
Tstamp <- format(Sys.time(), "%Y-%m-%d") 

# Load the list of folder structure paths
# - automated, no input required if proposed folder structure is followed
paths <- make_paths()

# Input path for summary (automated)
# - gap-filled and partitioned data
path_in <- list.files(paths$gap_filling, 
                      pattern = paste0(siteyear, ".*GF_essentials.*csv"),
                      full.names = TRUE)[1]

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

mean <- choose_avail(mean, names(data))
sum <- choose_avail(sum, names(data))
err_agg <- choose_avail(err_agg, names(data))

# Specify variables for plots at half-hour resolution
hh_vars <- grep("[^c]$", unique(c(mean, sum, err_agg)), value = TRUE)

# Print half-hourly results to pdf and png
pdf(file.path(
  paths$summary,
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
  paths$summary,
  paste0(siteyear, "_wind_roses_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd")]), ], 
         angle = 15, paddle = FALSE, breaks = 5)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "time")]), ], 
         type = "time", angle = 15, paddle = FALSE, breaks = 5)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "months")]), ], 
         type = "months", angle = 15, paddle = FALSE, breaks = 5, 
         grid.line = 10)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "stability")]), ], 
         type = "stability", angle = 15, paddle = FALSE, 
         breaks = 5, grid.line = 10, 
         main = "Zeta parameter based stability classes")
print(ggplot_stats(data, "wind_dir", "x_peak", circular = TRUE))
print(ggplot_stats(data, "wind_dir", "x_70perc", circular = TRUE))
print(ggplot_stats(data, "wind_dir", "wind_speed", circular = TRUE))
print(ggplot_stats(data, "wind_dir", "ustar", circular = TRUE))
print(ggplot_stats(data, "wind_dir", "zeta", circular = TRUE))
dev.off()

# Print separately to png
# - save_png() is defined in utilities.R as a helper function for saving plots
save_png("wind_rose_all", paths$png, siteyear, Tstamp)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd")]), ], 
         angle = 22.5, paddle = FALSE, breaks = 5)
dev.off()

save_png("wind_rose_day-night", paths$png, siteyear, Tstamp)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "time")]), ], 
         type = "time", angle = 22.5, paddle = FALSE, breaks = 5)
dev.off()

save_png("wind_rose_months", paths$png, siteyear, Tstamp)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "months")]), ], 
         type = "months", angle = 45, paddle = FALSE, breaks = 5, 
         grid.line = 10)
dev.off()

save_png("wind_rose_stability", paths$png, siteyear, Tstamp)
windRose(wrose_all[complete.cases(wrose_all[c("ws", "wd", "stability")]), ], 
         type = "stability", angle = 22.5, paddle = FALSE, 
         breaks = 5, grid.line = 10, 
         main = "Zeta parameter based stability classes")
dev.off()

save_png("wind_dir_x_peak", paths$png, siteyear, Tstamp)
print(ggplot_stats(data, "wind_dir", "x_peak", circular = TRUE))
dev.off()

save_png("wind_dir_x_70perc", paths$png, siteyear, Tstamp)
print(ggplot_stats(data, "wind_dir", "x_70perc", circular = TRUE))
dev.off()

save_png("wind_dir_wind_speed", paths$png, siteyear, Tstamp)
print(ggplot_stats(data, "wind_dir", "wind_speed", circular = TRUE))
dev.off()

save_png("wind_dir_ustar", paths$png, siteyear, Tstamp)
print(ggplot_stats(data, "wind_dir", "ustar", circular = TRUE))
dev.off()

save_png("wind_dir_zeta", paths$png, siteyear, Tstamp)
print(ggplot_stats(data, "wind_dir", "zeta", circular = TRUE))
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

# Save the summaries to CSV files
for (i in seq_along(means)) {
  summaries[[i]] <- cbind(means[[i]], fsd[[i]]$mean[-c(1:2)], 
                          DT_SD[[i]]$mean[-c(1:2)], sums[[i]][-c(1:2)], 
                          fsd[[i]]$sum[-c(1:2)], DT_SD[[i]]$sum[-c(1:2)],
                          pars[[i]][-c(1:2)])
  write_eddy(
    round_df(summaries[[i]]), 
    file.path(
      paths$summary,
      paste0(siteyear, "_", resol_names[i], "_summary_", Tstamp, ".csv"))
  )
}

### Plot summaries for different intervals ===================================== 

vars <- names(summaries$daily)
vars <- vars[!(vars %in% c("Intervals", "days"))]

# Plot aggregated variables per given intervals to pdf 
for (interval in c("daily", "weekly", "monthly")) {
  pdf(file.path(
    paths$summary,
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

### Plot spatio-temporal sampling coverage to pdf and png ======================

# - see Griebel et al. (2020) or ?spti_coverage for details

# Create a list of ggplot objects
spti_covp <- spti_coverage(
  data, "timestamp", "wind_dir", "NEE_uStar_f", "NEE_uStar_fqc",
  plot = TRUE)

# Save all ggplots to pdf
pdf(file.path(
  paths$summary,
  paste0(siteyear, "_spatio-temporal_sampling_coverage_", Tstamp, ".pdf")),
  width = 11.00, height = 8.27)
print(spti_covp)
dev.off()

# Save specified ggplots to png
save_png("spatial_sampling_coverage", paths$png, siteyear, Tstamp)
print(spti_covp[[1]]$spatial_sampling_coverage)
dev.off()

save_png("temporal_sampling_coverage", paths$png, siteyear, Tstamp)
print(spti_covp[[1]]$temporal_sampling_coverage)
dev.off()

# EOF
