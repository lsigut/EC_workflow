### Description ================================================================

# This code primarily aims to produce the summary of processed data
# at different timescales (daily, weekly, monthly, yearly) and to plot them. 
# Code developed by Ladislav Šigut (sigut.l@czechglobe.cz).  

### Setting paths and loading inputs ===========================================

# Firstly manually set the desired working directory in RStudio interface
# (typically it is the folder specifying the siteyear, e.g. BKF16)
# This assumes that the subdirectory structure is already present
# (i.e. the !Site_structure_template from the server)
# NB: the script does not support multiyear processing

# Load the package 
library(openeddy)

siteyear <- "KRP16" # Edit the site-year (included in folder and file names)
Tstamp <- format(Sys.time(), "%Y-%m-%d") # Timestamp of the computation
# (will be included in file names)

# Input and output path for summary (Level 3 data):
path_in <- "./Level 3/Gap-filling/REddyProc/KRP16_GF_essentials_2018-08-20.csv"
path_out <- "./Level 3/Summary/REddyProc/"
path_png <- "./Level 3/Summary/REddyProc/png/"

# Load the input file and convert timestamp and date:
# Only one year can be loaded; 
# e.g. May 2012-May 2013 could lead to incorrect computation
site <- read_eddy(path_in)
str(site) # check if loaded properly
site$timestamp <- strptime_eddy(site$timestamp, "%Y-%m-%d %H:%M", 
                                center.by = -900)
head(site$timestamp)
# NB: extraction from site$date would be incorrect at midnight:
site$date <- as.Date(site$timestamp) 

### Plot half-hourly data before unit conversion ===============================

# Specify variables for averaging and summation needed later (used for plotting)  
mean <- c("Tair", "Tair_fsd", "Tsoil", "Tsoil_fsd", "RH", "VPD", "VPD_fsd", 
          "ustar", "zeta", "H_fqc", "LE_fqc", "NEE_uStar_fqc")  
sum <- c("P", "GR", "Rn", "PAR", "H_f", "LE_f", 
         grep("(NEE|GPP).+f$", names(site), value = TRUE),
         grep("GPP_DT.+[^S][^D]$", names(site), value = TRUE),
         grep("Reco.+[^S][^D]$", names(site), value = TRUE))
err_agg <- grep("^[^TV].*(sd$|SD)$", names(site), value = TRUE)

hh_vars <- grep("[^c]$", c(mean, sum, err_agg), value = TRUE)

# Print half-hourly results to pdf and png
plot_hh <- function(x, name, pch = ".") {
  plot(x$timestamp, x[, name], pch = pch, xaxt = "n", xlab = "timestamp", 
       ylab = paste0(name, " [", openeddy::units(x[name]), "]"),
       main = name)
  r <- as.POSIXct(round(range(x$timestamp), "months"))
  axis.POSIXct(1, at = seq(r[1], r[2], by = "month"), format = "%b-%y")
}

pdf(paste(path_out, siteyear, "_half-hourly_plots_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
invisible(lapply(hh_vars, plot_hh, x = site))
dev.off()

for (i in hh_vars) {
  png(paste(path_png, siteyear, "_hh_", i, "_", Tstamp, ".png", sep = ""),
      width = 3508, height = 2480, res = 400)
  plot_hh(site, i)
  dev.off()
}

# Plot estimated parameters
pars <- grep("E_0_uStar|R_ref_uStar|^FP.*uStar$", names(site), value = TRUE)
pdf(paste(path_out, siteyear, "_parameter_estimates_", Tstamp, ".pdf", 
          sep = ""),
    width = 11.00, height = 8.27)
invisible(lapply(pars, plot_hh, x = site, pch = 20))
dev.off()

### Prepare variables for summaries ============================================

# Select carbon fluxes to convert units for summaries
carbon <- grep("NEE|GPP|Reco", c(sum, err_agg), value = TRUE)

# Perform sign correction in case the mean GPP is negative 
GPP <- grep("GPP", sum, value = TRUE)
sign_cor <- sapply(site[GPP], function(x) mean(x) < 0)
site[GPP][sign_cor] <- -site[GPP][sign_cor]

# Convert units
# http://alexis.czechglobe.cz/projects/vt04/knowledgebase/articles/10
# - from umol+1s-1m-2 to MJ+1hh-1m-2 (hh = half-hour)
site$PAR <- site$PAR * 1800e-6 / 4.57 # 4.57 Thimijan and Heins (1983)
units(site$PAR) <- "MJ+1hh-1m-2"
# - from W m-2 to MJ+1hh-1m-2 (hh = half-hour)
energy <- c("GR", "Rg_fsd", "Rn", "H_f", "H_fsd", "LE_f", "LE_fsd")
site[energy] <- lapply(site[energy], function(x) x * 1800e-6)
units(site[energy]) <- rep("MJ+1hh-1m-2", length(energy))
# - from umol+1s-1m-2 to g(C)+1hh-1m-2
site[carbon] <- lapply(site[carbon], function(x) x * 0.0216)
units(site[carbon]) <- rep("g(C)+1hh-1m-2", length(carbon))

### Daily summary ==============================================================

# Compute daily means
d_mean <- aggregate(site[mean], list(day = site$date), mean)
units(d_mean) <- c("-", units(site[mean]))

# Compute daily sums
d_sum <- aggregate(site[sum], list(day = site$date), sum)
d_sum_units <- units(site[sum])
d_sum_units <- gsub("hh", "day", d_sum_units)
units(d_sum) <- c("-", d_sum_units)

# Aggregate SD 'in quadrature' for variables gap-filled by look-up table
d_sd <- aggregate(site[err_agg], list(day = site$date), 
                  function(x) sqrt(sum(x^2)))
d_sd_units <- units(site[err_agg])
d_sd_units <- gsub("hh", "day", d_sd_units)
units(d_sd) <- c("-", d_sd_units)

# Compute relative error of bootstrapped variables
estim_err <- function(x, names) {
  fun <- function(x) round(abs((max(x) - min(x)) / median(x)) * 100, 2)
  l <- apply(names, 1, function(y) x[y])
  out <- as.data.frame(lapply(l, function(x) apply(x, 1, fun)))
  names(out) <- sapply(names[, 1], function(x) sub("_U05", "", x))
  names(out) <- paste0("rel_err_", names(out))
  units(out) <- rep("%", ncol(out))
  return(out)
}
xn <- data.frame(U05 = grep("(NEE|GPP|Reco).*U05", sum, value = TRUE),
                 U50 = grep("(NEE|GPP|Reco).*U50", sum, value = TRUE),
                 U95 = grep("(NEE|GPP|Reco).*U95", sum, value = TRUE),
                 stringsAsFactors = FALSE)
d_rel_err <- estim_err(d_sum, xn)

# Compute additional parameters
d_pars <- d_sum["day"]
d_pars$bowen_ratio_f <- d_sum$H_f / d_sum$LE_f
d_pars$evaporative_fraction <- d_sum$LE_f / (d_sum$LE_f + d_sum$H_f)
d_pars$closure_fraction <- (d_sum$LE_f + d_sum$H_f) / d_sum$Rn
d_pars$CUP <- FALSE
d_pars$CUP[d_sum$NEE_uStar_f < 0] <- TRUE
d_pars$GSL_method1 <- FALSE
d_pars$GSL_method1[d_mean$Tair > 5] <- TRUE

# Merge the daily means, sums and pars  
daily <- cbind(d_mean, d_sum[-1], d_sd[-1], d_rel_err, d_pars[-1])

### Weekly summary =============================================================

# - week starts with Monday
weeks <- strftime(site$date, format = "%W-%y")
weeks <- factor(weeks, levels = unique(weeks))

# Compute weekly means
w_mean <- aggregate(site[mean], list(week = weeks), mean)
units(w_mean) <- c("-", units(site[mean]))

# Compute weekly sums
w_sum <- aggregate(site[sum], list(week = weeks), sum)
w_sum_units <- units(site[sum])
w_sum_units <- gsub("hh", "week", w_sum_units)
units(w_sum) <- c("-", w_sum_units)

# Aggregate SD 'in quadrature' for variables gap-filled by look-up table
w_sd <- aggregate(site[err_agg], list(week = weeks), 
                  function(x) sqrt(sum(x^2)))
w_sd_units <- units(site[err_agg])
w_sd_units <- gsub("hh", "week", w_sd_units)
units(w_sd) <- c("-", w_sd_units)

# Compute relative error of bootstrapped variables
w_rel_err <- estim_err(w_sum, xn) 

# Compute additional parameters
w_pars <- w_sum["week"]
w_pars$bowen_ratio_f <- w_sum$H_f / w_sum$LE_f
w_pars$evaporative_fraction <- w_sum$LE_f / (w_sum$LE_f + w_sum$H_f)
w_pars$closure_fraction <- (w_sum$LE_f + w_sum$H_f) / w_sum$Rn
# - requires daily resolution
weeks_d <- strftime(daily$day, format = "%W-%y")
weeks_d <- factor(weeks_d, levels = unique(weeks))
w_pars$CUP <- aggregate(d_pars["CUP"], list(week = weeks_d), sum)$CUP
w_pars$GSL_method1 <- aggregate(d_pars["GSL_method1"], list(week = weeks_d), 
                                sum)$GSL_method1

# Compute number of days within interval
w_days <- aggregate(daily$day, list(week = weeks_d), length)
names(w_days)[2] <- "days_in_week"

# Merge the weekly means, sums and pars  
weekly <- cbind(w_days, w_mean[-1], w_sum[-1], w_sd[-1], w_rel_err, w_pars[-1])

### Monthly summary ============================================================

months <- strftime(site$date, format = "%m-%y")
months <- factor(months, levels = unique(months))

# Compute monthly means
m_mean <- aggregate(site[mean], list(month = months), mean)
units(m_mean) <- c("-", units(site[mean]))

# Compute monthly sums
m_sum <- aggregate(site[sum], list(month = months), sum)
m_sum_units <- units(site[sum])
m_sum_units <- gsub("hh", "month", m_sum_units)
units(m_sum) <- c("-", m_sum_units)

# Aggregate SD 'in quadrature' for variables gap-filled by look-up table
m_sd <- aggregate(site[err_agg], list(month = months), 
                  function(x) sqrt(sum(x^2)))
m_sd_units <- units(site[err_agg])
m_sd_units <- gsub("hh", "month", m_sd_units)
units(m_sd) <- c("-", m_sd_units)

# Compute relative error of bootstrapped variables
m_rel_err <- estim_err(m_sum, xn) 

# Compute additional parameters
m_pars <- m_sum["month"]
m_pars$bowen_ratio_f <- m_sum$H_f / m_sum$LE_f
m_pars$evaporative_fraction <- m_sum$LE_f / (m_sum$LE_f + m_sum$H_f)
m_pars$closure_fraction <- (m_sum$LE_f + m_sum$H_f) / m_sum$Rn
# - requires daily resolution
months_d <- strftime(daily$day, format = "%m-%y")
months_d <- factor(months_d, levels = unique(months))
m_pars$CUP <- aggregate(d_pars["CUP"], list(month = months_d), sum)$CUP
m_pars$GSL_method1 <- aggregate(d_pars["GSL_method1"], list(month = months_d), 
                                sum)$GSL_method1

# Compute number of days within interval
m_days <- aggregate(daily$day, list(month = months_d), length)
names(m_days)[2] <- "days_in_month"

# Merge the monthly means, sums and pars  
monthly <- cbind(m_days, m_mean[-1], m_sum[-1], m_sd[-1], m_rel_err, m_pars[-1])

### Yearly summary =============================================================

years <- as.numeric(strftime(site$date, format = "%Y"))

# Compute yearly means
y_mean <- aggregate(site[mean], list(year = years), mean)
units(y_mean) <- c("-", units(site[mean]))

# Compute yearly sums
y_sum <- aggregate(site[sum], list(year = years), sum)
y_sum_units <- units(site[sum])
y_sum_units <- gsub("hh", "year", y_sum_units)
units(y_sum) <- c("-", y_sum_units)

# Aggregate SD 'in quadrature' for variables gap-filled by look-up table
y_sd <- aggregate(site[err_agg], list(year = years), 
                  function(x) sqrt(sum(x^2)))
y_sd_units <- units(site[err_agg])
y_sd_units <- gsub("hh", "year", y_sd_units)
units(y_sd) <- c("-", y_sd_units)

# Compute relative error of bootstrapped variables
y_rel_err <- estim_err(y_sum, xn) 

# Compute additional parameters
y_pars <- y_sum["year"]
y_pars$bowen_ratio_f <- y_sum$H_f / y_sum$LE_f
y_pars$evaporative_fraction <- y_sum$LE_f / (y_sum$LE_f + y_sum$H_f)
y_pars$closure_fraction <- (y_sum$LE_f + y_sum$H_f) / y_sum$Rn
# - requires daily resolution
years_d <- as.numeric(strftime(daily$day, format = "%Y"))
y_pars$CUP <- aggregate(d_pars["CUP"], list(year = years_d), sum)$CUP
y_pars$GSL_method1 <- aggregate(d_pars["GSL_method1"], list(year = years_d), 
                                sum)$GSL_method1
# Method 2 for estimating start and end of growing season: 
# if Tair > 5 and Tsoil > 0 for 3 consecutive days;
# the first day out of this sequence is assigned GS_start. 
# The GS_end is assigned the first day that holds 
# this condition when checked backwards from the last day of year 
# (from personal experience works for Bily Kriz forest).
if (any(is.na(daily[c("Tair", "Tsoil")]))) {
  GS_start <- NA 
  GS_end <- NA
} else {
  GS_start <- daily$day[1]
  repeat {
    rows <- which(daily$day == GS_start)
    if (min(daily[c("Tair", "Tsoil")][rows:(rows + 2), ] > 
            matrix(ncol = 2, c(rep(5, 3), rep(0, 3)))) > 0) break
    GS_start <- GS_start + 1
  }
  GS_end <- daily$day[nrow(daily)]
  repeat {
    rows <- which(daily$day == GS_end)
    if (min(daily[c("Tair", "Tsoil")][rows:(rows - 2), ] > 
            matrix(ncol = 2, c(rep(0, 3), rep(5, 3)))) > 0) break
    GS_end <- GS_end - 1
  }
}
y_pars$GS_start_method2 <- GS_start 
y_pars$GS_end_method2 <- GS_end
y_pars$GSL_method2 <- as.numeric(GS_end - GS_start)

# Method 3 for estimating start and end of growing season: 
# if Tair > 5 for 5 consecutive days; the first day out of this sequence 
# is assigned GS_start. The GS_end is assigned the first day that holds this 
# condition when checked backwards from the last day of year 
# (Marek et al., 2011; for original citation ask Radek Pokorny).
if (any(is.na(daily["Tair"]))) {
  GS_start <- NA 
  GS_end <- NA
} else {
  GS_start <- daily$day[1]
  repeat {
    rows <- which(daily$day == GS_start)
    if (min(daily["Tair"][rows:(rows + 4),] > 5) > 0) break
    GS_start <- GS_start + 1
  }
  GS_end <- daily$day[nrow(daily)]
  repeat {
    rows <- which(daily$day == GS_end)
    if (min(daily["Tair"][rows:(rows - 4),] > 5) > 0) break
    GS_end <- GS_end - 1
  }
}  
y_pars$GS_start_method3 <- GS_start
y_pars$GS_end_method3 <- GS_end
y_pars$GSL_method3 <- as.numeric(GS_end-GS_start)

# Compute number of days within interval
y_days <- aggregate(daily$day, list(year = years_d), length)
names(y_days)[2] <- "days_in_year"

# Merge the yearly means, sums and pars  
yearly <- cbind(y_days, y_mean[-1], y_sum[-1], y_sd[-1], y_rel_err, y_pars[-1])

### Saving all computed data =================================================== 
summary_path <- function(x) {
  paste(path_out, siteyear, "_", x, "_summary_", Tstamp, ".csv", sep = "")
}
write_eddy(daily, summary_path("daily"))
write_eddy(weekly, summary_path("weekly"))
write_eddy(monthly, summary_path("monthly"))
write_eddy(yearly, summary_path("yearly"))

# Plotting computed data ======================================================= 
vars <- c(mean, sum, err_agg, names(y_rel_err), 
          "bowen_ratio_f", "evaporative_fraction", "closure_fraction")
plot_summary <- function(x, interval = c("daily", "weekly", "monthly")) {
  interval <- match.arg(interval)
  d <- interval == "daily"
  format <- if (d) "%F" else if (interval == "weekly") "%W-%y" else "%b-%y"
  p <- barplot(as.numeric(x[, i]), space = 0,
               width = if (interval == "daily") 1 else x[, 2], 
               names.arg = if (d) NULL else {
                 unique(strftime(site$date, format = format))
               },
               ylim = if (i == "closure_fraction") c(0, 1) else NULL,
               xlab = paste(interval, "timescale"),
               ylab = paste(i, " [", units(x[i]), "]", sep = ""),
               main = i, axis.lty = if (d) 0 else 1)
  if (d) {
    index <- round(seq(1, length(p), l = 8))
    axis(1, at = p[index],
         labels = unique(strftime(site$date, format = format))[index])
  }
}

# Plot to pdf
pdf(paste(path_out, siteyear, "_daily_plots_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
for (i in vars) {
  plot_summary(daily, "daily")
}
dev.off()

pdf(paste(path_out, siteyear, "_weekly_plots_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
for (i in vars) {
  plot_summary(weekly, "weekly")
}
dev.off()

pdf(paste(path_out, siteyear, "_monthly_plots_", Tstamp, ".pdf", sep = ""),
    width = 11.00, height = 8.27)
for (i in vars) {
  plot_summary(monthly, "monthly")
}
dev.off()

# Plot to png
for (i in vars) {
  png(paste(path_png, siteyear, "_daily_", i, "_", Tstamp, ".png", 
            sep = ""),
      width = 3508, height = 2480, res = 350)
  plot_summary(daily, "daily")
  dev.off()
}

for (i in vars) {
  png(paste(path_png, siteyear, "_weekly_", i, "_", Tstamp, ".png", 
            sep = ""),
      width = 3508, height = 2480, res = 350)
  plot_summary(weekly, "weekly")
  dev.off()
}

for (i in vars) {
  png(paste(path_png, siteyear, "_monthly_", i, "_", Tstamp, ".png", 
            sep = ""),
      width = 3508, height = 2480, res = 350)
  plot_summary(monthly, "monthly")
  dev.off()
}

# EOF
