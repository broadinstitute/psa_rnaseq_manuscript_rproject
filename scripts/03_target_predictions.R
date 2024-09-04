#Use moc1430_1779_0066_0110_2kd_modzscore_ref_set_n250x5679.gct reference
#This example uses external data as query

library(cmapR)
library(dplyr)
library(ggplot2)
library(tidyr)

source('./functions/Functions_psa_MOA.R')
source('./functions/Functions_general.R')

#####User input parameters#################
metric = "nzscore" #Choose zscore or nzscore based on the metric used for query data, only affects the positive predictive values listed in the final table.
outdir = '/broad/hptmp/jbagnall/PRJNA291292/'
screen_name = "PRJNA291292"
project_id0 = "PRJNA291292"
date0 = "231128"
query_data_path = './data_tables/external_data/PRJNA291292_mod_nzscore_n15x5077.gct'

###Get path to reference data######################
ref_data_path = './data_tables/antimicrobial_reference_set/moc1430_1779_0066_0110_2kd_modzscore_ref_set_n250x5679.gct'

save_combined_path = paste0(outdir, screen_name, '_ref_combined_modzscore_', date0)
save_cormat_path = paste0(outdir, screen_name, '_ref_combined_modzscore_cor_mat_', date0)
combine_and_make_cor_matrix(ref_data_path = ref_data_path, query_data_path = query_data_path, save_combined_path = save_combined_path, save_cor_path = save_cormat_path)

#get list of reference IDs
ref_data = parse_gctx(ref_data_path)
reference_ids = ref_data@cid
remove(ref_data)

#########Load annotation############
anno_path = './reference_files/Metadata_compound_240512.csv'
anno = read.csv(anno_path, stringsAsFactors = F)

########Make correlation matrix############
#Find exact path to correlation file
files_in_outdir = list.files(outdir)
cor_file_name = files_in_outdir[grepl(paste0(screen_name, '_ref_combined_modzscore_cor_mat_', date0), files_in_outdir)]
cormat_path = paste0(outdir,cor_file_name)
cor_mat = parse_gctx(cormat_path)
cdesc = cor_mat@cdesc
rdesc = cor_mat@rdesc
cor_mat = cor_mat@mat

#Extract query IDs from column metadata 
# query_df = filter(cdesc, tolower(project_id) == tolower(project_id0) & !(tolower(pert_id) %in% c("water", "dmso")))
query_data = parse_gctx(query_data_path)
query_df = query_data@cdesc
query_ids = query_df$id
query_pert_ids = unique(query_df$pert_id)

######Predict target#####################
colname1 = "pert_target" #This is the column we will try to predict on
colname2 = "pert_mechanism"

#Keep query compounds in reference even if they're the same as the query compounds
remove_query0 = FALSE #keep as false unless doing a leave one out cross validation
if(remove_query0){
  remove_query_label = "remove_query_comp"
}else{
  remove_query_label = "keep_query_comp"
}

filename1 = paste0(outdir, screen_name, '_refbased_moa_target_', remove_query_label, "_", date0, '.rds')
if(file.exists(filename1)){
  print("Warning: file exists, will overwrite")
  print(filename1)
}

#This function for each treatment (comp_conc), gives the maximally correlated ref compound per target category (colname1)
targetCorrelation(cormatrix_path = cormat_path, query_ids = query_ids, reference_ids = reference_ids, sample_meta = cdesc, compound_meta = anno, target_colname1 = colname1, target_colname2 = colname2, removequery = remove_query0, savefilename = filename1)

####Make reference-based MOA predictions
#Chooses which file to use to report positive predictive values (PPV) from leave-one-compound-out analysis, based on zscores or nzscores
if(metric == "nzscore"){
  roc_path1 = './reference_files/kabx_rocinfo_avg_cor_compound_include_singles_queryposcon_nzscore_vs_zscore_230919_pert_target_add_NA.rds'
  roc_path2 = './reference_files/kabx_rocinfo_avg_cor_queryposcon_nzscore_vs_zscore_compound_include_singles_240328_pert_target_mechanism_add_NA.rds'
}else if(metric == "zscore"){
  roc_path1 = './reference_files/kabx_rocinfo_avg_cor_by384well_bystrain_compound_include_singles_230918_zscore_vs_zscore_pert_target_add_NA.rds'
  roc_path2 = './reference_files/kabx_rocinfo_avg_cor_by384well_bystrain_compound_include_singles_240328_zscore_vs_zscore_pert_target_mechanism_add_NA.rds'
}else{
  print("Could not match metric values. PPV will be based on zscores.")
  roc_path1 = './reference_files/kabx_rocinfo_avg_cor_by384well_bystrain_compound_include_singles_230918_zscore_vs_zscore_pert_target_add_NA.rds'
  roc_path2 = './reference_files/kabx_rocinfo_avg_cor_by384well_bystrain_compound_include_singles_240328_zscore_vs_zscore_pert_target_mechanism_add_NA.rds'
}

avg_doses1 = T #True averages the maximal correlation of each target category across treatments (e.g. doses, time points, etc.)
if(avg_doses1){
  avg_label = "avg_doses"
}else{
  avg_label = "max_across_trts"
}

#Make a list of predictions
savefilepath = outdir
savefilename1 = paste0(savefilepath, paste(screen_name,'refbasedMOA', avg_label, remove_query_label, date0, colname1, sep = "_"))
if(file.exists(savefilename1)){
  print("Warning: file exists, will overwrite")
  print(savefilename1)
}

#Make prediction table
refBasedList_compound(target_list_path = filename1, rank_thresh = 1, kabx_annopath = anno, target_col = colname1, roc_path = roc_path1, roc_path2 = roc_path2, normalized_cor = F, savefilename = savefilename1, print_plot = F, avg_doses = avg_doses1)
print(paste0("Output prediction files start with: ", savefilename1))

##########Plots###########
cor_df_path = paste(savefilename1, "_unknown_all_target_cor.rds", sep = "")

for(num in 1:length(query_pert_ids)){
  pert_id0 = query_pert_ids[num]
  savefilepath = paste0(outdir, pert_id0, '_', metric, '_', date0, '.pdf') 
  temp = plot_refBasedList(pert_id0 = pert_id0, all_max_cor_df = cor_df_path, kabx_anno_path = anno_path, text_size = 12, save_file_path = savefilepath)
}
