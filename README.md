# psa_rnaseq_manuscript_rproject
Scripts used for PsA RNAseq manuscript, included in an Rproject file.

# Scripts
Scripts include:
<ol>
  <li>01_qc_and_vst.Rmd: Script for QC reads and calcuating vst</li>
  <li>02a_calculate_mod_zscore.Rmd: Calculate moderated z-scores (if you have variety of drug treatments)</li>
  <li>02a_calculate_mod_nzscore.Rmd: Calculate moderated nz-scores (z-scores relative to negative control)</li>
  <li>03_target_predictions.R: Predict targets of query compounds</li>
</ol>

# Functions
Functions are R functions utilized by the scripts above.

# Data Tables
The data_tables folder includes the z-scores or nz-scores of the various data sets in the manuscript
* reference_set: z-score values of treatments using drugs with known mechanisim. Here, batch is defined as 384-well plate, strain, and time point
* internal_test_compounds: z-score or nz-score values of treatments of compounds with unknown mechanism. Tested internally.
* external_data: nz-score values of publicly available data GSE166602
* hypomorph_zscore_across_strains: z-score values of a subset of the reference set to assess hypomorph transcriptional signatures, where batch is re-defined as 384 well plate and time-point (allowing strains to be different within the batch)
* crispri_strains: nz-scores for CRISPRi strains

# Reference Files
Reference files used for annotation or calculating positive predictive value
* Compound metadata for reference compounds
* kabx_rocinfo...zscore.rds is an input file for running target prediction function for queries that are z-scores, for assigning positive predictive values
* kabx_rocinfo...nzscore_vs_zscore.rds is an input file for running target prediction function for queries that are nz-scores, for assigning positive predictive values

