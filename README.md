<!-- README.md is generated from README.Rmd. Please edit that file -->
EC\_workflow
============

The proposed workflow allows to process eddy covariance data with single processing chain consisting of four stages:

1.  **Quality control (QC):** load the EddyPro output and gap-filled meteorological data and apply automated tests and filters implemented in [openeddy](https://github.com/lsigut/openeddy) to quality check fluxes of momentum (Tau), sensible (H) and latent heat (LE) and net ecosystem exchange (NEE). Export documentation of applied QC and produce the outputs needed in next steps. QC workflow is placed in `.\Level 2\Quality checking\`.
2.  **Storage:** under development. Currently it only plots monthly means of daily courses for all storage fluxes using [REddyProc](https://github.com/bgctw/REddyProc). The plan is to include processing of temperature, RH and CO2 profiles and computation of storage fluxes with support for 1 Hz measurements. Storage workflow is placed in `.\Level 2\Storage flux\`.
3.  **Gap-filling and flux partitioning (GF & FP):** combine utilities of [REddyProc](https://github.com/bgctw/REddyProc) and [openeddy](https://github.com/lsigut/openeddy) to gap-fill (H, LE, NEE), partition (NEE) and visualize (H, LE, NEE) fluxes. The setup allows to change and document some processing options in an organized way. Computation of bootstrapped friction velocity threshold is included. GF & FP workflow is placed in `.\Level 3\Gap-filling\REddyProc\`.
4.  **Summary:** visualize processed data, convert units and aggregate results to daily, weekly, monthly and yearly timescales. A limited amount of computed parameters is also produced, including different uncertainty estimates. Summary workflow is placed in `.\Level 3\Summary\REddyProc\`.

The example dataset is available at Zenodo due to GitHub file size limitation: <https://doi.org/10.5281/zenodo.1442531>

The workflows assume certain folder structure:

![](Folder%20structure.jpg)

-   Level 0: raw data, related metadata and configuration files
-   Level 1: half-hourly data processed by EddyPro and gap-filled meteorological data
-   Level 2: results and documentation of QC, storage corrected fluxes for GF & FP
-   Level 3: results of GF & FP and the dataset summaries

The complete processing chain in context of above folder structure can be summarized as:

![](Processing%20chain.jpg)

In order to run fetch filter, QC workflow also requires fetch filter boundary vector for given site, here loaded from `.\Level 2\Quality checking\Fetch_filter_boundaries_20160206.csv`. For your site you will have to define your own region of interest (ROI) boundary (see *ROI boundary* section below).

Example
-------

In order to run the example siteyear KRP16:

1.  Download `KRP16 - before processing.zip` from [Zenodo](https://doi.org/10.5281/zenodo.1442531) and unzip
2.  Install devtools package if not available yet

    ``` r
    install.packages("devtools")
    ```

3.  Install [openeddy](https://github.com/lsigut/openeddy)

    ``` r
    devtools::install_github("lsigut/openeddy")
    ```

4.  Install [REddyProc](https://github.com/bgctw/REddyProc)

    ``` r
    install.packages("REddyProc")
    ```

5.  Go through the processing chain by running commands in workflow files in the order and at the locations as described above
6.  You can check your results with those of `KRP16 - processed.zip` at [Zenodo](https://doi.org/10.5281/zenodo.1442531)

Note: Sourcing the QC workflow does not produce desired outcome because interactive functions are included.

ROI boundary
------------

The spatial extent of the studied ecosystem (region of interest; ROI) is specified by its ROI boundary that describes the distance from EC tower to the edge of the studied ecosystem. In order to work with [openeddy](https://github.com/lsigut/openeddy), ROI boundary has to be provided as a numeric vector with following properties:

-   The number of circular sectors is the same as the number of provided distances (length of the vector)
-   The angular resolution of the ROI boundary is given by `360° / number of angular sectors`
-   The ROI boundary distances are assigned to the centers of their respective circular sectors with first sector centered on 0°

### ROI boundary example

![](ROI%20boundary%20example.jpg)

In this simplified case ROI boundary would be specified as:

``` r
c(150, 200, 250, 300)
```

**Interpretation:**

-   There would be 4 circular sectors with 90° angular resolution
-   ROI boundary is specified for the whole first sector (315°, 45°\] at the distance 150 m from tower (center of the sector is 0°)
-   Boundary of the second sector (45°, 135°\] is at the distance 200 m
-   Third sector (135°, 225°\] is at the distance 250 m
-   Fourth sector (225°, 315°\] is at the distance 300 m

Realistic representation of ROI boundary can look e.g. like this:

![](ROI%20boundary%20RAJ.jpg)

References
----------

Publication describing openeddy is not yet available. When describing the proposed quality control scheme, please refer to is as similar to:

Mauder, M., Cuntz, M., Drüe, C., Graf, A., Rebmann, C., Schmid, H.P., Schmidt, M., Steinbrecher, R., 2013. A strategy for quality and uncertainty assessment of long-term eddy-covariance measurements. Agric. For. Meteorol. 169, 122-135, <https://doi.org/10.1016/j.agrformet.2012.09.006>

The methodology and benchmark of REddyProc 1.1.3 is descibed in the following paper:

Wutzler, T., Lucas-Moffat, A., Migliavacca, M., Knauer, J., Sickel, K., Šigut, L., Menzer, O., and Reichstein, M. (2018): Basic and extensible post-processing of eddy covariance flux data with REddyProc, Biogeosciences, 15, 5015-5030, <https://doi.org/10.5194/bg-15-5015-2018>.
