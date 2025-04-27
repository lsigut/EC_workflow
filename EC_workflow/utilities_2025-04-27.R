### Description ================================================================

# Collection of utility functions that seem to have merit only in the given
# workflow and may or may not be moved to openeddy package. To simplify the
# workflow for 4 different fluxes high level functions were designed and saved
# in 'utilities.R'. This is mostly when multiple commands should be run without
# the need of user intervention. User can still adapt function arguments.
#
# Code developed by Ladislav Sigut (sigut.l@czechglobe.cz).

#' Supported Fluxes
#'
#' A complete set of supported fluxes in eddy covariance workflow.
#'
#' While users can utilize [openeddy] to create their own processing workflow
#' for any set of fluxes or variables, the original intent was to process
#' specified fluxes. The processing workflow is available here:
#' <https://github.com/lsigut/EC_workflow>.
fluxes <- c("Tau", "H", "LE", "NEE")

# Attach an R package
# package: A character string specifying single package
# github: A character string specifying github repository
# - used in all workflows to attach or install and attach packages if not present
attach_pkg <- function(package, github = NULL) {
  # load package if available
  avail <- require(package, character.only = TRUE)
  if (!avail) {
    # if not available and placed at github
    if (!is.null(github)) {
      # devtools package is required
      if (!require("devtools")) install.packages("devtools")
      devtools::install_github(github)
    } else {
      # if not available and placed at CRAN
      install.packages(package)
    }
    # load installed package
    require(package, character.only = TRUE)
  }
}

# Name merged output
# EP_path: A character string specifying folder name including EddyPro data
# siteyear: A character string specifying siteyear
# - used in data_preparation workflow
name_merged <- function(EP_path, siteyear) {
  # Names of merged output files
  # - data in CSV files and documentation in TXT files
  data_name_out <- list.files(EP_path)
  data_name_out <- grep("[.][Cc][Ss][Vv]$", data_name_out, value = TRUE)
  if (length(data_name_out) == 1) {
    data_name_out <- gsub("[.][Cc][Ss][Vv]$", "_met.csv", data_name_out)
  } else {
    data_name_out <- paste0("eddypro_", siteyear, 
                            "_full_output_merged_adv_met.csv")
  }
  return(data_name_out)
}

#' Combine Documentation
#'
#' Read documentation from single or multiple TXT files. In case of multiple
#' files, combine them together with one additional line separating them.
#'
#' @param path A character vector. The full paths to TXT files.
#'
#' @seealso \code{\link{readLines}}.
#'
#' @export
#' - used in document_merged() below
combine_docu <- function(path) {
  unlist(lapply(path, function(x) c(readLines(x, warn = FALSE), "")))
}

# Document merged files
# - this function is called for its side effect - writing TXT documentation file
# data_name_out: A character string. Name of merged output file.
# EP_path: A character string specifying folder name including EddyPro data
# Meteo_path: A character string specifying folder name including Meteo data
# out_path: A character string specifying folder name for output files
# Tstamp: A character string specifying timestamp of the computation
# name, mail: character string with contact information
# M: A data frame with merged Meteo data 
# - used in data_preparation workflow
document_merged <- function(data_name_out, EP_path, Meteo_path, out_path,
                            Tstamp, name, mail, M) {
  docu_name_out <- gsub("[.][Cc][Ss][Vv]$", "\\.txt", data_name_out)
  
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
    fp <- file.path(out_path, docu_name_out)
    message("saving file to ", fp)
    writeLines(c(paste0(Tstamp, ":"),
                 paste0("Files merged by ", name, " (", mail, ")"),
                 "",
                 "Merged files:", 
                 M_names,
                 EP_names,
                 "",
                 "Variables from meteo database remapped to:",
                 paste(names(M), varnames(M), sep = " = ", collapse = "\n"), 
                 "", 
                 combine_docu(docu_name_in),
                 "Information about the R session:",
                 capture.output(sessionInfo())), 
               fp, sep = "\n")
  }
}

# Save plots of precheck variables in a single pdf to specified path
# data: A data frame with column names and "timestamp" column in POSIXt format.
# precheck: A character vector of available precheck variables.
# siteyear: A character string specifying siteyear
# Tstamp: A character string specifying timestamp of the computation
# path: A character string specifying folder name for saving the pdf
# width, height: The width and height of the graphics region in inches.
# qrange: A numeric vector of length 2, giving the quantile range of y-axis.
# - used in QC workflow
save_precheck_plots <- function(data, precheck, siteyear, Tstamp, path, 
                                width = 11.00, height = 8.27, 
                                qrange = c(0.005, 0.995)) {
  fp <- file.path(path, 
                  paste0(siteyear, "_auxiliary_precheck_", Tstamp, ".pdf"))
  message("saving file to ", fp)
  pdf(fp, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  invisible(lapply(precheck, plot_precheck, x = data, qrange = qrange))
}

# Save plots of fluxes with meteo in separate pdfs to specified path
# data: A data frame with column names and "timestamp" column in POSIXt format.
# qc_suffix: A character string identifying respective QC flag included in data.
# siteyear: A character string specifying siteyear
# sname: A character string to be evaluated by sprintf and %s substituted for flux  
# Tstamp: A character string specifying timestamp of the computation
# path: A character string specifying folder name for saving the pdf
# fluxes: A character vector of supported flux names
# width, height: The width and height of the graphics region in inches.
# - used in QC workflow
save_flux_plots <- function(data, qc_suffix = "prelim", siteyear, sname,
                            Tstamp, path, fluxes, 
                            width = 11.00, height = 8.27) {
  for (i in fluxes) {
    fp <- file.path(
      path,
      paste0(siteyear, "_", sprintf(sname, i), "_", Tstamp, ".pdf"))
    message("saving file to ", fp)
    pdf(fp, width = width, height = height)
    qc <- paste("qc", i, qc_suffix, sep = "_")
    plot_eddy(data, i, qc, qc)
    dev.off()
  }
}

# Show independent or cumulative effect of all filters 
# data: A data frame with column names. 
# prelim: A tibble with names of quality control flags to combine
# cumul: A logical value that determines if cumulative (cumul = TRUE) or
#   individual (cumul = FALSE) effects of quality control flags should be shown.
# - used in QC workflow
plot_QC_summary <- function(data, prelim, cumul) {
  gridExtra::grid.arrange(grobs = lapply(names(prelim), function(x) 
    summary_QC(data, na.omit(pull(prelim, x)), cumul = cumul, plot = TRUE, 
               flux = x)), 
    nrow = 2)
}

# Save QC summary plots produced by plot_QC_summary()
# data: A data frame with column names. 
# prelim: A tibble with names of quality control flags to combine
# path: A character string specifying folder name for saving the png
# siteyear: A character string specifying siteyear
# Tstamp: A character string specifying timestamp of the computation
# width, height: The width and height of the graphics region in inches.
# - used in QC workflow
save_QC_summary_plots <- function(data, prelim, path, siteyear, Tstamp,
                                  width = 297, height = 210) {
  fp_ind <- file.path(
    path,
    paste0(siteyear, "_QC_summary_", Tstamp, ".png"))
  message("saving file to ", fp_ind)
  ggsave(fp_ind,
    plot_QC_summary(data, prelim, cumul = FALSE),
    type = "cairo-png", width = width, height = height, units = "mm")
  fp_cum <- file.path(
    path,
    paste0(siteyear, "_QC_summary_cumulative_", Tstamp, ".png"))
  message("saving file to ", fp_cum)
  ggsave(fp_cum,
    plot_QC_summary(data, prelim, cumul = TRUE),
    type = "cairo-png", width = width, height = height, units = "mm")
}

# Combine specified flags for given flux to produce preliminary flags
# - informative naming is useful e.g. for despiking and manual QC
# data: A data frame with quality control flags specified in prelim
# prelim: A tibble with names of quality control flags to combine
# - used in QC workflow
combn_prelim_QC <- function(data, prelim) {
  res <- sapply(names(prelim), 
                function(x) combn_QC(data, na.omit(pull(prelim, x))))
  res <- as.data.frame(res)
  names(res) <- paste("qc", names(res), substitute(prelim), sep = "_")
  return(res)
}

# Set names of existing manual quality control columns
# fluxes: A character vector containing names of supported fluxes
# man: A data frame with manual quality control flags
# - there might be variables without manual QC or 'man' can be NULL
# - used in QC workflow
set_man_names <- function(fluxes, man) {
  mnames <- paste0("qc_", fluxes, "_man")
  names(mnames) <- fluxes
  is.na(mnames) <- !(mnames %in% names(man))
  mnames <- as.list(mnames)
  return(mnames)
}

# Document quality control step
# - this function is called for its side effect - writing TXT documentation file
# Tstamp: A character string specifying timestamp of the computation
# name, mail: character string with contact information
# strg_applied: A logical value documenting whether storage correction was applied
# forGF: A tibble with names of quality control flags used for final QC
# path: A character string specifying folder name for saving the TXT file
# siteyear: A character string specifying siteyear
# - used in QC workflow
document_QC <- function(Tstamp, name, mail, strg_applied, forGF,
                        path, siteyear) {
  fp <- file.path(
    path,
    paste0(siteyear, '_QC_info_', Tstamp, '.txt'))
  message("saving file to ", fp)
  writeLines(c(paste0(Tstamp, ":"),
               paste0("Quality controlled by ", name, " (", mail, ")"),
               "",
               paste0("Storage corrected fluxes: ", strg_applied),
               "",
               paste0("Applied w_rot correction: ", applied_w_rot_correction),
               "",
               "Applied quality control scheme:", 
               capture.output(as.data.frame(forGF)),
               "",
               "Information about the R session:",
               capture.output(sessionInfo())), 
             fp, 
             sep = "\n")
}

# Document gap-filling and flux partitioning
# - this function is called for its side effect - writing TXT documentation file
# all_out: A data frame with column names containing REddyProc exported results
# Tstamp: A character string specifying timestamp of the computation
# name, mail: character string with contact information
# siteyear: A character string specifying siteyear
# lat, long, tzh: Numeric values specifying latitude, longtitude and timezone
# FP_temp: A character string. Temperature used for flux partitioning
# seasonal_ustar: A logical value. Was ustar threshold resolved seasonally? 
# use_CPD: A logical value. Was change-point detection method used?
# path: A character string specifying folder name for saving the TXT file
# - used in GF workflow
document_GF <- function(all_out, Tstamp, name, mail, siteyear, lat, long, tzh,
                        FP_temp, fixed_UT, seasonal_ustar, use_CPD, path) {
  # compute flux availability percentage
  perc_records <- nrow(all_out) / 100
  flux_avail_names <- c("H_orig", "LE_orig", "NEE_uStar_orig")
  avail_rec <- lapply(flux_avail_names, function(x) {
    temp <- table(!is.na(all_out[x]))
    round(unname(temp["TRUE"] / perc_records), 1)
  })
  names(avail_rec) <- flux_avail_names
  
  # document GF output
  fp <- file.path(
    path,
    paste0(siteyear, '_documentation_', Tstamp, '.txt'))
  message("saving file to ", fp)
  writeLines(c(paste0(Tstamp, ":"),
               paste0("Processed by ", name, " (", mail, ")"),
               "",
               paste0("Siteyear:"),
               siteyear,
               "",
               "Used site metadata:",
               paste0("latitude = ", lat, ", longitude = ", long, ", timezone = ", 
                      tzh),
               "",
               "Temperature used for flux partitioning:",
               FP_temp,
               "",
               "Ustar filtering settings:",
               if (is.na(fixed_UT)) {
                 c(paste0("seasonal_ustar = ", seasonal_ustar),
                   paste0("use_changepoint_detection = ", use_CPD))
               } else {
                 paste("fixed_ustar_threshold:", fixed_UT, "m s-1")
               },
                "",
               "Availability of original records for respective flux:",
               paste0("H = ", avail_rec$H_orig, "%"),
               paste0("LE = ", avail_rec$LE_orig, "%"),
               paste0("NEE = ", avail_rec$NEE_uStar_orig, "%"),
               "",
               "Information about the R session:",
               capture.output(sessionInfo())), 
             fp, 
             sep = "\n")
}

# Create utility function for saving plots to png
# x: A character string specifying naming of the given plot
# path: A character string specifying folder name for saving the png
# siteyear: A character string specifying siteyear
# Tstamp: A character string specifying timestamp of the computation
# width, height: The width and height of the graphics region in inches
# res: An integer specifying png resolution (see ?png)
# - used in Summary workflow
save_png <- function(x, path, siteyear, Tstamp, width = 3508, height = 2480, 
                     res = 400) {
  png(file.path(
    path,
    paste0(siteyear, "_", x, "_", Tstamp, ".png")),
    width = width, height = height, res = res)
}

# EOF
