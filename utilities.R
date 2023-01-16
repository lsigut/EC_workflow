### Description ================================================================

# Collection of utility functions that seem to have merit only in the given
# workflow and may or may not be moved to openeddy package.
#
# Code developed by Ladislav Sigut (sigut.l@czechglobe.cz).

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
               file.path(out_path, docu_name_out), sep = "\n")
  }
}

# Create utility function for saving plots to png
# - x: A character string specifying naming of the given plot
# - used in Summary workflow
png_util <- function(x, width = 3508, height = 2480, res = 400) {
  png(file.path(
    paths$png,
    paste0(siteyear, "_", x, "_", Tstamp, ".png")),
    width = width, height = height, res = res)
}

# Construct saving path for different aggregation intervals
# - x: a character string with interval description (e.g. "daily")
# - used in Summary workflow
filenames <- function(x) {
  file.path(
    paths$Summary,
    paste0(siteyear, "_", x, "_summary_", Tstamp, ".csv"))
}

# EOF
