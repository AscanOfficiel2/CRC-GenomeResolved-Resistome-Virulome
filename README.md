# CRC-GenomeResolved-Resistome-Virulome
Genome-resolved metagenomic pipeline and latent functional modeling framework for resistance–virulence organization in colorectal cancer microbiomes.


# Overview

This repository contains the complete computational workflow described
in:

**Ascandari A, Aminu S, Benhida R, Daoud R.**\
*Resistance-enriched latent functional organization of the gut
microbiome underlies colorectal cancer progression (under review).*

The pipeline performs genome-resolved metagenomic reconstruction,
resistance and virulence profiling, and latent functional modeling
associated with colorectal cancer (CRC) progression.

------------------------------------------------------------------------

# Repository Structure

    01_Bioninformatics/
    02_Preprocessing/
    03_ARG_VF_analysis/
    04_Species_analysis/

Each module corresponds to a major analytical stage. Below is a concise
description of the purpose of each module and the scripts it contains.The modules and their scripts have been numbered in sequencial order to facilitate stepwise execution for the user.

------------------------------------------------------------------------

# 01_Bioninformatics

## Purpose

This module performs genome-resolved metagenomic assembly, binning,
dereplication, taxonomic annotation, and functional profiling (ARGs and
VFs). It produces curated MAGs and species-level abundance matrices used
in downstream analyses.

## Scripts

-   `01_run.sh`\
-   `02_Assembly_Binning.sh`\
-   `03_collect_files.sh`\
-   `04_global_drep_taxonomy_abundance.sh`\
-   `05_global_ARG.sh`\
-   `06_global_VFDB.sh`\
-   `07_clean_abundance_files.sh`\
-   `08_merge_abundances.sh`\
-   `09_mapping_to_species.sh`

------------------------------------------------------------------------

# 02_Preprocessing

## Purpose

This module ensures MAG quality control, evaluates assembly metrics,
filters contamination, estimates coverage, and harmonizes metadata prior
to functional and statistical modeling.

## Scripts

-   `00_gather_script.sh`\
-   `01_mag_quality_crc_new.py`\
-   `02_quast_crc.py`\
-   `03_fastqscreen_crc_mags.py`\
-   `04_Nonpareil_CRC.R`\
-   `05_mags_vf_arg_cleanup.py`\
-   `06_metadata_analysis.py`

------------------------------------------------------------------------

# 03_ARG_VF_analysis

## Purpose

This module conducts ARG and VF preprocessing, compositional
normalization, batch correction, diversity analysis, enrichment testing,
and multivariate latent modeling to identify resistance-enriched
functional structures associated with CRC progression.

## Scripts

-   `01_arg_vf_preprocessing.py`\
-   `02_VF_ARG_genes_CLR_Batch_correction.R`\
-   `03_Beta_diversity_arg_vf_genes.R`\
-   `04_ARG_VF_prevalence.R`\
-   `05_ARG_VF_cross_talk.R`\
-   `06_pls_dml.py`\
-   `07_PLS_Components_stats_Analysis.R`\
-   `08_Genes_log2foldchange.R`\
-   `09_Mechanism_enrichment.R`\
-   `10_pls_high_risk_pathobiont.py`\
-   `11_Maaslin3.R`

------------------------------------------------------------------------

# 04_Species_analysis

## Purpose

This module encodes species-level signals of CRC-associated functional
organization and evaluates predictive performance, generalizability, and
feature stability using shared and unshared species frameworks.

## Scripts

-   `01_species_matrix_processing.py`\
-   `02_Overlap_with_reference_panel.R`\
-   `03_Benchmark_shared_unshared.sh`\
-   `04_Nested_cv_holdout_shared_unshared.sh`\
-   `05_Comparison_Delong.sh`\
-   `06_Feature_importance_shared.sh`\
-   `07_Feature_importance_unshared.sh`

------------------------------------------------------------------------

# Software Requirements

## Core Tools

-   MEGAHIT
-   MetaBAT2
-   dRep
-   GTDB-Tk
-   CheckM2
-   QUAST
-   BBTools
-   FastQ Screen
-   Nonpareil
-   RGI (CARD)
-   DIAMOND (VFDB)
-   CoverM

## Programming Environments

-   **R** (vegan, limma, MaAsLin3)\
-   **Python** (scikit-learn, pandas, numpy)

Exact software versions are specified in the manuscript Methods section.

------------------------------------------------------------------------

# Data Availability

All sequencing datasets analyzed in this study are publicly available
through the **NCBI Sequence Read Archive (SRA)**.\
Accession numbers are provided in Supplementary Data S1 (In manuscript) and also uploaded in this repository as "SRA_meta_Data".

------------------------------------------------------------------------

# Citation

If you use or adapt this workflow, please cite:

Ascandari A et al.\
*Resistance-enriched latent functional organization of the gut
microbiome underlies colorectal cancer progression.*

------------------------------------------------------------------------

# Contact

**AbdulAziz Ascandari**\
abdulaziz.ascandari@um6p.ma\
ryandari87@gmail.com

**Suleiman Aminu**\
suleiman.aminu@um6p.ma\
saminu83@gmail.com

