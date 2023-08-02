## EC_workflow 2023-04-14

-   Automate w_rot correction
-   Update Zenodo link
-   README: include link to openeddy tutorials

## EC_workflow 2022-03-18

-   Update Zenodo files, workflows and description
-   Add utilities file with functions related to workflow
-   README: describe Manual QC and how to load meteo
-   Rename workflows to reflect their order
-   Include flags for checking instruments in Czechglobe-specific
    section
-   Increase angular resolution of wind roses
-   Support G, SWC and GWL variables
-   Assure that ggplots are printed also when sourcing
-   Specify auxiliary vars in check_manually()
-   Update requirements for package version
-   README: include new QC filters
-   Czechglobe specific branch removed after formalizing instrument
    filters
-   Formalize documentation functions and move them to utilities
-   Simplify QC workflow - formalize functions and move to utilities
-   Document w_rot correction
-   Inform about saving files
-   Amend GF workflow and move GF documentation to utilities
-   Add support for fixed Ustar threshold

## EC_workflow 2018-09-17

-   Include workflow version date suffix
-   Documentation improvements
-   Gap-fill VPD if any record is NA
-   Include workflow description in README
-   Provide link to datasets through Zenodo
-   Include ROI boundary description and figs in README

## EC_workflow (no date suffix)

### 2018-09-09

-   Update QC workflow with manual QC and remove spikesHF & wresid
    filter
-   Correct units in summary (fractions & ratios)

### 2018-09-03

-   Convert NEE to NEP
-   Implement support for Meteo vars replicates

### 2018-08-31

-   Implement automated merging of EddyPro & Meteo

### 2018-08-28

-   FCO2 renamed to NEE for both storage corrected and uncorrected
    fluxes to unify workflows
-   Unify QC for sites with/without storage correction
-   Include wind rose plotting
