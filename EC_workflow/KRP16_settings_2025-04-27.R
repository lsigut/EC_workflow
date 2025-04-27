### Description ================================================================

# Settings for all EC workflow files are stored here for easier workflow updates 
# and setup editing. General settings are provided first, settings relevant to
# particular workflow can be found in the respective section.
#
# Code developed by Ladislav Sigut (sigut.l@czechglobe.cz).

# Contact information
name <- "Ladislav Sigut" # person that performed processing
mail <- "sigut.l@czechglobe.cz" # mail of the person that performed processing

# Edit the siteyear
# - included in folder and file names
siteyear <- "KRP16"

# Edit the start and end of time period to post-process
# - start <- 2016; end <- 2016: to assure complete year (recommended)
# - start <- NULL; end <- NULL: use timestamp extent from input files
# - start <- "2016-02-01 10:00:00"; end <- "2017-02-01 10:00:00": arbitrary period
start <- 2016
end <- 2016

# Specify the time shift (in seconds) to be applied to the date-time information
# in order to represent the center of averaging period
shift.by <- -900

# Specify site metadata
lat <- 49.5732575 # edit site latitude
long <- 15.0787731 # edit site longtitude
tzh <- 1 # timezone hour


# WF1: data_preparation ========================================================

# Set Meteo mapping
# - edit Meteo mapping according to available variables at given site
# - regular expressions (?regex) can be used to select replicates of given 
#   variable (e.g. Ts) that will be averaged by remap_vars()
# - typical set of variables: global radiation (GR), photosynthetically active
#   radiation (PAR), Net radiation (Rn), air temperature (Tair), soil 
#   temperature (Tsoil), relative humidity (RH), vapor pressure dificit (VPD),
#   soil water content (SWC), precipitation (P), soil heat flux (G)
# - variables required by openeddy: PAR, Rn, Tair, Tsoil, VPD, P (only PAR
#   really needed for optimal setup, the rest can be initialized with NA values)
# - variables required by REddyProc: GR, Tair, Tsoil, VPD (and/or RH)	
Met_mapping <- tribble(
  ~Meteo_varname, ~workflow_varname,
  "date/time",       "timestamp",
  "GRin",            "GR",
  "PARin",           "PAR",
  "^Rn",             "Rn",
  "TaKA02.0",        "Tair",
  "TsKO0.05",        "Tsoil",
  "RHKA02.0",        "RH",
  "VPDKA02.0",       "VPD",
  "sumP",            "P"
)

# Provide timestamp format of your Meteo data (see ?strptime_eddy)
# - default: "%Y-%m-%d %H:%M"
meteo_format <- "%Y-%m-%d %H:%M"

# WF2: QC (quality control) ====================================================

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

# Set Fetch filter boundary 
# - edit the region of interest boundary for your site
# - here KRP boundary version 20160206
# - see ROI boundary concept description at https://github.com/lsigut/EC_workflow
# - it can be determined manually (e.g. using Google Earth app)
# - R scripts for obtaining ROI boundary based on defined KML file can be found
#   here: https://github.com/lsigut/ROI_boundary 
boundary <- 
  c(453, 489, 469, 455, 444, 410, 375, 348, 86, 82, 78, 76, 74, 73, 72, 72, 73, 
    74, 76, 78, 81, 85, 91, 97, 106, 116, 114, 113, 131, 372, 496, 500, 507, 
    519, 531, 541, 555, 562, 565, 572, 584, 605, 633, 749, 863, 1012, 1128, 
    1098, 802, 863, 871, 903, 403, 360, 328, 303, 283, 486, 466, 451, 441, 412, 
    390, 373, 360, 350, 349, 356, 367, 381, 399, 422)

# Precheck variables
# - set of variables typically available in EddyPro full output that can be
#   useful for preliminary check before quality control procedure
# - the selected variables will be plotted (no further dependencies)
precheck_vars <- c(
  "u_rot", "v_rot", "w_unrot", "w_rot",
  "sonic_temperature", "max_wind_speed",
  "Tau", "ustar", "H", "LE", "NEE",
  "u_var", "v_var", "w_var", "ts_var", "h2o_var", "co2_var",
  "rand_err_Tau", "rand_err_H", "rand_err_LE", "rand_err_NEE",
  "Tau_scf", "H_scf", "LE_scf", "co2_scf",
  "u_spikes", "v_spikes", "w_spikes", "ts_spikes", "co2_spikes", "h2o_spikes",
  "H_strg", "LE_strg", "co2_strg",
  "h2o_v_adv", "co2_v_adv",
  "co2_mixing_ratio", "h2o_mixing_ratio",
  "co2_time_lag", "h2o_time_lag",
  "x_peak", "x_70perc",
  "mean_value_RSSI_LI_7200", "co2_signal_strength_7200_mean",
  "h2o_signal_strength_7200_mean", "flowrate_mean")

# Quality Control Essential Variables
# - a set of variables useful when working with quality controlled data
# - some variables might have dependencies in further processing steps
essential_vars_QC <- c(
  "timestamp", "GR", "qc_GR", "PAR", "qc_PAR", "Rn", "qc_Rn", "Tair",
  "qc_Tair", "Tsoil", "qc_Tsoil", "RH", "qc_RH", "VPD", "qc_VPD", "SWC",
  "qc_SWC", "WTD", "qc_WTD", "GWL", "qc_GWL", "P", "qc_P", "G", "qc_G", "Tau",
  "Tau_orig", "qc_Tau_forGF", "qc_Tau_SSITC", "rand_err_Tau", "H", "H_orig",
  "qc_H_forGF", "qc_H_SSITC", "rand_err_H", "LE", "LE_orig", "qc_LE_forGF",
  "qc_LE_SSITC", "rand_err_LE", "NEE", "NEE_orig", "qc_NEE_forGF",
  "qc_NEE_SSITC", "rand_err_NEE", "H_strg", "LE_strg", "co2_strg", "wind_speed",
  "wind_dir", "ustar", "L", "zeta", "model", "x_peak", "x_70perc")

# WF3: GF_&_FP (gap filling and flux partitioning) =============================

# Temperature used for flux partitioning ('Tsoil' or 'Tair') 
# - default: FP_temp <- 'Tsoil'
# - note that MDS gap-filling is based on 'Tair'
FP_temp <- 'Tsoil'

# Save figures as "png" (default) or "pdf"
# - NEEvsUStar plots are fixed to "pdf"
plot_as <- "png" 

# Set fixed ustar threshold if needed (skip ustar threshold estimation)
# - default is: fixed_UT <- NA; i.e. UT estimation is recommended
# - if provided, seasonal_ustar and use_CPD settings will be ignored
fixed_UT <- NA

# Choose ustar threshold estimation method resolution
# - either seasonal ustar threshold (seasonal_ustar <- TRUE; default)
# - or annual thresholds (seasonal_ustar <- FALSE)
# - seasonal UT is recommended as it keeps more data (see Wutzler et al., 2018)
seasonal_ustar <- TRUE

# Should change-point detection (CPD) method (Barr et al. 2013) be used 
# (use_CPD <- TRUE) instead of classical approach (Reichstein et al. 2005, 
# Papale et al. 2006) using binned classes of ustar and temperature? 
# (use_CPD <- FALSE; default)
# The changepoint is estimated based on the entire subset within one season 
# and one temperature class; currently using argument 'isUsingCPTSeveralT' 
# of function usControlUstarEst()
use_CPD <- FALSE

# WF4: Summary =================================================================

# Specify variables needed for different procedures during aggregation 
# - used for averaging, summation, uncertainty estimation and plotting
# - some variables have dependencies in further processing steps
# - consider adding instead of removing from the list
# - GWL (ground water level), SWC (soil water content) and G (soil heat flux)
#   are not available in KRP16 example data set (they will be reported and 
#   skipped)
mean <- c("Tair", "Tsoil", "RH", "VPD", "GR", "Rn", "PAR", "GWL", "SWC", "G", 
          "H_f", "H_fqc", "LE_f", "LE_fqc", "ET_f", "ET_fqc", 
          "NEE_uStar_f", "NEE_uStar_fqc", "GPP_uStar_f", "GPP_DT_uStar", 
          "Reco_uStar", "Reco_DT_uStar")
sum <- c("P", "GR", "Rn", "PAR", "G", "H_f", "LE_f", "ET_f", 
         grep("(NEE|GPP).+f$", names(data), value = TRUE),
         grep("GPP_DT.+[^D]$", names(data), value = TRUE),
         grep("Reco.+[^D]$", names(data), value = TRUE))
err_agg <- grep("(^H|^LE|^ET|^NEE|^Reco|^GPP).*(sd|SD)$", names(data), 
                value = TRUE)
