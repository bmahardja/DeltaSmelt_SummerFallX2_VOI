---
editor_options: 
  markdown: 
    wrap: 72
---

# DeltaSmelt_SummerFallX2_VOI

Repository for Delta Smelt Summer-Fall X2 Value of Information analysis

## CalSim_output

Script to read and process dss files/output from CalSim3 model runs

## Consequence_Table

Folder containing summarized data from all IBMR runs and scripts to run
multi-criteria decision analysis (MCDA) and value of information (VOI)
calculations

-   **Data**: Formatted data from BoRVOI_IBMR_final_results excel sheet
    files
-   **Output**: Figures produced from scripts

## Figures

Files needed to construct map and draft workflow figure

-   **Map**: Data and script to produce map figure
-   **Workflow**: Draft modeling workflow schematics not used in
    manuscript

## IBMR

Data and code for Delta Smelt IBMR (Individual Based Model - R)

-   **data**: raw data for IBMR
    -   **data_raw**: original model inputs
    -    **data_processed**: adjusted zooplankton and X2 data for IBMR
        runs
    -   **output**: folder for IBMR run results
-   **scripts:** scripts for running IBMR in parallel mode
    (control.model script to run dependencies)
    -   **orig_scripts**: original IBMR scripts provided by Will Smith
        (IBMR) without parallelization

## Salinity_Zooplankton_Model

Data and code for IBMR sub-models

-   **Data**: Data for salinity-zooplankton models
-   **Output**: Figures from salinity-zooplankton models
-   **R**: Scripts for converting CalSim3 data to expected salinity for
    each IBMR region
