FibeRtypeR Analysis Template — README

Purpose:
This script runs a standardized analysis to quantify the relationship between skeletal muscle fiber type (from FibeRtypeR) and insulin sensitivity/resistance outcomes, using pre-specified models (A–F) for meta-analysis. It is designed to run locally in each lab (GDPR-friendly) and return only model results, not raw data.
Author: Max Ullrich
Last updated: 2026-01-16

**Quick start (for collaborators)**


1. Place your files

  - Put your metadata CSV and FibeRtypeR results CSV somewhere accessible.
  - Copy this script and the renv.lock file into a local project folder.



2. Open R/RStudio in that folder.


3. Configure the script
   Open the script and edit the CONFIG (EDIT) block near the top:

    - metadata_path – path to your metadata CSV
    - fibertype_path – path to your FibeRtypeR results CSV
    - study_id – a short identifier (e.g., a PMID)
    - outcome_var – one of "HOMA_IR", "M_value", "Matsuda_index"
    - (Optional) baseline_filter_col + baseline_filter_val – set to NULL to disable
    - Column mapping (column_map) – map your column names to the standardized names; set any that you don’t have to NULL.

Important: If you have a disease grouping variable for Model E, map disease_state to your column name (otherwise Model E is skipped).





4. Run the analysis by running these two commands in your console:

       source("Analysis_template.R")

       run_all()


  The script will:

  - Restore packages from renv.lock (first run may take a few minutes)
  - Run eligible models
  - Save outputs in ./<study_id>_outputs/



5. Return only the results
Send back the full "outputs" folder
