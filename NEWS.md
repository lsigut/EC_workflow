## EC_workflow 2025-04-27

-   bug fix: selection of variables for WF_4 aggregation in settings now
    follows literal list - previous definition had not available
    dependencies (2025-07-30)
-   all workflow files moved to common folder EC_workflow
-   site-year prefix is now kept only for settings file as workflow
    files should not need editing
-   improved efficiency of folder structure and command line support of
    folder names by the implementation of make_paths()
-   automated loading of inputs using make_paths()
-   Meteo data is expected as a single file with a header and units that
    can be read by read_eddy()
-   settings for all workflow files were extracted to a single file
    (siteyear_settings.R) for easier workflow updates and setup editing
-   README: updated folder structure description and graphics
-   README: new sections QC principles and Flagging scheme
-   README: new section Change of EC system or EddyPro settings
-   allow flexibility of setting processing period (start and end
    variables instead of year)
-   automated setting of optimal TempRange in nighttime partitioning
    (REddyProc package)

## EC_workflow 2023-04-14

-   automate w_rot correction
-   update Zenodo link
-   README: include link to openeddy tutorials

## EC_workflow 2022-03-18

-   update Zenodo files, workflows and description
-   add utilities file with functions related to workflow
-   README: describe Manual QC and how to load meteo
-   rename workflows to reflect their order
-   include flags for checking instruments in Czechglobe-specific
    section
-   increase angular resolution of wind roses
-   support G, SWC and GWL variables
-   assure that ggplots are printed also when sourcing
-   specify auxiliary vars in check_manually()
-   update requirements for package version
-   README: include new QC filters
-   Czechglobe specific branch removed after formalizing instrument
    filters
-   formalize documentation functions and move them to utilities
-   simplify QC workflow - formalize functions and move to utilities
-   document w_rot correction
-   inform about saving files
-   amend GF workflow and move GF documentation to utilities
-   add support for fixed Ustar threshold

## EC_workflow 2018-09-17

-   include workflow version date suffix
-   documentation improvements
-   gap-fill VPD if any record is NA
-   include workflow description in README
-   provide link to datasets through Zenodo
-   include ROI boundary description and figs in README

## EC_workflow (no date suffix)

### 2018-09-09

-   update QC workflow with manual QC and remove spikesHF & wresid
    filter
-   correct units in summary (fractions & ratios)

### 2018-09-03

-   convert NEE to NEP
-   implement support for Meteo vars replicates

### 2018-08-31

-   implement automated merging of EddyPro & Meteo

### 2018-08-28

-   FCO2 renamed to NEE for both storage corrected and uncorrected
    fluxes to unify workflows
-   unify QC for sites with/without storage correction
-   include wind rose plotting
