### Loading the required packages ==================================

# You might need to install the packages first. In that case use:
# install.packages("packagename")
library("eddyczechr")
library("REddyProc")

### Setting paths and loading inputs ===========================================

# Firstly manually set the desired working directory in RStudio interface
# (typically it is the folder specifying the siteyear, e.g. BKF16)
# This assumes that the subdirectory structure is already present
# (i.e. the !Site_structure_template from the server)

siteyear <- "KRP16" # Edit the site-year (included in folder and file names)
Tstamp <- format(Sys.time(), "%Y-%m-%d") # Timestamp of the computation
# (will be included in file names)

# Output path for storage plotting (Level 2 data):
(path <- "./Level 2/Storage flux/")

# Load the input files. Edit the filenames:
# - timestamp (combines date and time) and meteo data are expected
input <- paste("./Level 2/Input for gap-filling/",
               "KRP16_forGF_QC_essentials_2017-09-06.csv", 
               sep = "")
# -set correct skip parameter (number of lines above header in the input)
# -set correct file encoding (e.g. fileEncoding = "UTF-8") 
site <- read_eddy(input)
str(site) # check if loaded properly

### Convert timestamp ==========================================================

# Convert timestamp to POSIXct format without centering:
# Note: use GMT which just means "don’t modify times and dates – just use them
# exactly as they are irrespective of timezone information of your computer 
# system", i.e. use as is
site$DateTime <- strptime_eddy(site$timestamp, "%Y-%m-%d %H:%M")
class(site$DateTime) # should be POSIXt class now
head(site$DateTime) # check if converted and centered correctly

### Initalize R5 reference class sEddyProc for storage plotting ================
strg <- c("H_strg", "LE_strg", "h2o_strg", "co2_strg") # storage flux names
site_R5 <- sEddyProc$new(siteyear, site, strg)

### Plot diurnal cycle of storage fluxes =======================================
for (Var in strg) {
  site_R5$sPlotDiurnalCycle(Var, Dir.s = path, Format.s = 'png')
}



