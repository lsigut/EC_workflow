### Description ================================================================

# Collection of utility functions that seem to have merit only in the given
# workflow and may or may not be moved to openeddy package.
#
# Code developed by Ladislav Sigut (sigut.l@czechglobe.cz).

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
