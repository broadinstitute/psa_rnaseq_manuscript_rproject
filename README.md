# Table of Contents
Scripts and data files used for PsA RNAseq manuscript

## R project file
Open psa_rnaseq_manuscript.Rproj file to easily access all scripts and data mentioned below.

## Scripts
* 01_qc_and_vst.Rmd: Script for QC reads and calcuating vst</li>
* 02a_calculate_mod_zscore.Rmd: Calculate moderated z-scores (if you have variety of drug treatments)</li>
* 02a_calculate_mod_nzscore.Rmd: Calculate moderated nz-scores (z-scores relative to negative control)</li>
* 03_target_predictions.R: Predict targets of query compounds</li>

## Functions
Functions are R functions utilized by the scripts above.

## Data Tables
The data_tables folder includes the z-scores or nz-scores of the various data sets in the manuscript
* antimicrobial_reference_set: z-score values of treatments using drugs with known mechanism. Here, batch is defined as 384-well plate, strain, and time point
* internal_test_compounds: z-score or nz-score values of treatments of compounds with unknown mechanism. Tested internally.
* external_data: z-score and nz-score values of publicly available data PRJNA291292
* hypomorph_zscore_across_strains: z-score values of a subset of the reference set to assess hypomorph transcriptional signatures, where batch is re-defined as 384 well plate and time-point (allowing strains to be different within the batch)
* crispri_strains: z-scores for CRISPRi strains, subset by transcriptional profiles used in the CRISPRi reference set

## Reference Files
Reference files used for annotation or inputs to functions
* Compound metadata for reference compounds
* kabx_rocinfo...zscore...rds is an input file for running target prediction function for queries that are z-scores, for assigning target positive predictive values
* kabx_rocinfo...zscore...mechanism...rds is an input file for mechanism prediction for queries that are z-scores, for assigning mechanism positive predictive values
* kabx_rocinfo...nzscore_vs_zscore...rds is an input file for running target prediction function for queries that are nz-scores, for assigning target positive predictive values
* kabx_rocinfo...nzscore_vs_zscore...mechanism...rds is an input file for mechanism prediction for queries that are nz-scores, for assigning mechanism positive predictive values
