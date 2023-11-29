#Use moc1430_1779_0066_0110_2kd_modzscore_ref_set_n250x5679.gct reference
#This example uses external data as query

library(cmapR)
library(dplyr)
library(ggplot2)
library(tidyr)

source('/home/unix/jbagnall/git/psa_rnaseq_manuscript/functions/Functions_psa_MOA.R')
source('/home/unix/jbagnall/git/psa_rnaseq_manuscript/functions/Functions_general.R')

outdir = '/broad/hptmp/jbagnall/GSE166602/'
screen_name = "GSE166602"
project_id0 = "PRJNA701395"
date0 = "231128"

###Get paths to data######################
ref_data_path = '/idi/cgtb/jbagnall/psa_RNAseq/moc1430_1779_0066_0110_2kd_combined/zscore_mocp0110_by_strain_pa14subset/moc1430_1779_0066_0110_2kd_mod_zscore_by384well_bystrain_n333x5679.gct'
query_data_path = '/idi/cgtb/jbagnall/psa_RNAseq/external_data/GSE139257/nzscore/GSE139257_modzscore_remove_low_count_samples_and_genes_rep_collapse_clipped500_n4x5893.gct'

save_combined_path = paste0(outdir, screen_name, '_moc1430_1779_0066_0110_2kd_combined_modzscore_', date0)
save_cormat_path = paste0(outdir, screen_name, '_moc1430_1779_0066_0110_2kd_combined_modzscore_cor_mat_', date0)
combine_and_make_cor_matrix(ref_data_path = ref_data_path, query_data_path = query_data_path, save_combined_path = save_combined_path, save_cor_path = save_cormat_path)


#check reference set from moc1430_1779_0066_01102kd, should be the same as removing the mocp0110 from the other reference set
ref1_path = '/idi/cgtb/jbagnall/psa_RNAseq/moc1430_1779_0066_0110_2kd_combined/zscore_mocp0110_by_strain_pa14subset/moc1430_1779_0066_0110_2kd_zscore_by384well_bystrain_poscon_reference_w_negcon_ids_remove_inactives_230918.txt'
ref1 = read.table(ref1_path)

#notes: 
#Remove "transcriptional inactives" based on max(-log10(padj)) < 20, where the threshold was taken between water vs dmso
#In addition to removing nonactive treatments from poscon categories, include negcon when predicting on unknowns


#########Load annotation############
anno_path = '/idi/cgtb/jbagnall/psa_RNAseq/annotation/Metadata_compound_JB_230626.csv'
anno_collapse = read.csv(anno_path, stringsAsFactors = F)

########Make correlation matrix############
cormat_path = '/idi/cgtb/jbagnall/psa_RNAseq/external_data/GSE139257/GSE1392257_moc1430_1779_0066_0110_2kd_combined_modzscore_cor_mat_231101_n337x337.gctx'
cor_mat = parse_gctx(cormat_path)
cdesc = cor_mat@cdesc
rdesc = cor_mat@rdesc
cor_mat = cor_mat@mat


#######using Gates refbasedmoa scripts###################################

######Calculate correlation values
colname1 = "pert_target"
colname2 = "pert_mechanism"
filename1 = paste0(outdir, screen_name, '_refbased_moa_target_keep_query_comp_', date0, '.rds')
if(file.exists(filename1)){
  print("Warning: file exists, will overwrite")
  print(filename1)
}

reference_ids = read.table(ref1_path)
reference_ids = reference_ids$V1

query_ids = filter(cdesc, project_id == project_id0 & pert_id != "dmso")
query_ids = query_ids$id

remove_query0 = FALSE

#This function for each treatment (comp_conc), gives the maximally correlated kabx compound per target category (colname1)
targetCorrelation(cormatrix_path = cormat_path, query_ids = query_ids, reference_ids = reference_ids, metadata = cdesc, target_colname1 = colname1, target_colname2 = colname2, removequery = remove_query0, savefilename = filename1)

temp = readRDS(filename1)

####Make reference-based MOA predictions
#This needs to change to nz-score vs zscore of the reference set alone
roc_path1 = '/idi/cgtb/jbagnall/psa_RNAseq/moc1430_1779_0066_0110_2kd_combined/zscore_mocp0110_by_strain_pa14subset/nzscore/kabx_rocinfo_avg_cor_compound_include_singles_queryposcon_nzscore_vs_zscore_230919_pert_target.rds'
tani_roci_path1 = ""
avg_doses1 = T #False to use Max instead of Average. Usually it's TRUE, we use average

#Make a list of predictions
savefilepath = outdir
savefilename1 = paste0(savefilepath, paste(screen_name,'refbasedMOA', 'avg_doses_keep_query_comp', date0, colname1, sep = "_"))
if(file.exists(savefilename1)){
  print("Warning: file exists, will overwrite")
  print(savefilename1)
}

# savefilename1 = paste0(savefilepath, paste(screen_name,'refbasedMOA', method0, 'avg_doses_oprl_hypo', date0, colname1, sep = "_"))
# if(file.exists(paste(savefilename1, "unknown.csv", sep = "_"))){
#   print("Warning: file exists, will overwrite")
#   print(savefilename1)
# }

#average actual correlation values, don't weight each treatment the same. Use actual correlation values makes more sense here so ill-behaving strains don't get equal weight
anno_collapse2 = read.csv(anno_path, stringsAsFactors = F)
refBasedList_compound(target_list_path = filename1, rank_thresh = 1, kabx_annopath = anno_collapse2, target_col = colname1, roc_path = roc_path1, tani_path = '', tani_roc_path = '', normalized_cor = F, savefilename = savefilename1, print_plot = T, avg_doses = avg_doses1)

temp1 = refBasedList_compound(target_list_path = filename1, rank_thresh = 1, kabx_annopath = anno_collapse2, target_col = colname1, roc_path = roc_path1, tani_path = '', tani_roc_path = '', normalized_cor = F, savefilename = NA, print_plot = T)




##########Plots###########
# #Figure out which ones are correct
# target_maxcor_table = readRDS(filename1)
# target_maxcor_table = mutate(target_maxcor_table, query_pert_id = ifelse(grepl("BRD", query_compound), substr(query_compound, 1, 13), query_compound))
# target_maxcor_table = left_join(target_maxcor_table, anno_collapse2, by = c("query_pert_id" = "pert_id"))
# targets_test= read.csv('/idi/cgtb/jbagnall/psa_RNAseq/poscon_correlation/poscon_target_prediction_maxcor_test_220606.csv', stringsAsFactors = F)
targets = read.csv(paste(savefilename1, "_known.csv", sep = ""), stringsAsFactors = F)
sum(targets$same_target)
dim(targets)
length(unique(targets$query_pert_id))
length(unique(targets$query_target))
length(unique(targets$predicted_pert_mechanism))
length(unique(targets$predicted_pert_target))
length(unique(targets$predicted_pharmaceutical_class))
colname_target = colname1

target_list = select(targets, query_target, query_pert_id) %>% distinct()
target_list2 = sort(table(targets$query_target))
singletons = target_list2[target_list2 == 1]


#bar plot for each category
target_summary = targets %>%
  group_by(query_target) %>%
  summarize(called_correctly = sum(same_target), called_correctly_cor_06 = sum(same_target & target_correlation >= 0.6), called_correctly_cor_05 = sum(same_target & target_correlation >= 0.5), called_incorrectly = sum(!same_target), total_n  = n()) %>%
  ungroup()
target_summary$check = target_summary$total_n - target_summary$called_correctly
all(target_summary$called_incorrectly == target_summary$check)

sum(target_summary$called_correctly_cor_05)
sum(target_summary$total_n)
temp = filter(target_summary, total_n > 1)
sum(temp$called_correctly_cor_05)
sum(temp$total_n)


ggplot(target_summary) +
  geom_bar(aes(x = reorder(query_target, -total_n), y = total_n), stat = "identity", fill = "gray")+
  geom_bar(aes(x = reorder(query_target, -total_n), y = called_correctly), stat = "identity", fill = "blue")+
  theme_bw(base_size = 14)+
  labs(x = colname_target, y = "# compounds", title = "KABX (blue = predicted correctly)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(target_summary) +
  geom_bar(aes(x = reorder(query_target, -total_n), y = total_n), stat = "identity", fill = "gray")+
  geom_bar(aes(x = reorder(query_target, -total_n), y = called_correctly_cor_06), stat = "identity", fill = "blue")+
  theme_bw(base_size = 14)+
  labs(x = colname_target, y = "# compounds", title = "KABX (blue = correlation >= 0.6 & predicted correctly)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(target_summary) +
  geom_bar(aes(x = reorder(query_target, -total_n), y = total_n), stat = "identity", fill = "gray")+
  geom_bar(aes(x = reorder(query_target, -total_n), y = called_correctly_cor_05), stat = "identity", fill = "blue")+
  theme_bw(base_size = 14)+
  labs(x = colname_target, y = "# compounds", title = "KABX (blue = correlation >= 0.5 & predicted correctly)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


###Correlation plots#####################
#######Test correlation report functions
#Add pert_time and pert_itime to the cdesc
cdesc = mutate(cdesc, pert_time = ifelse(project_id==project_id0, 90, pert_time))
cdesc = mutate(cdesc, pert_itime = ifelse(project_id==project_id0, "90min", pert_itime))
cdesc$dose_rank = as.numeric(cdesc$dose_rank_global)
cdesc = mutate(cdesc, pert_iname = ifelse(is.na(pert_iname), pert_id, pert_iname))

comp1 = "ciprofloxacin" #This only works becuase strain_id is still in the id
conc1 = "0.06ug_mL"
strain_id1 = "PA14"
pert_itime1 = "90min"
project_id1 = project_id0
pert_iname1 = "ciprofloxacin"
target1 = ""
description1 = ""
comp_anno_path = '/idi/cgtb/jbagnall/psa_RNAseq/annotation/Metadata_compound_JB_230626.csv'
plot_size_col = "dose_rank"
plot_shape_col = "strain_id"
plot_fill_col = "pert_time"
compcorthreshold1 = 0



#If you want to keep the CRISPRi strains in the reference set, change pert_id, otherwise they filter out with pert_id = ara
# colmeta_combined = mutate(colmeta_combined,  pert_id = ifelse(strain_id %in% c("PA14_c26", "PA14_murA_c119", "PA14_lpxC_c165"), strain_id, pert_id))
corlist1 = corList(pert_id1 = comp1, pert_idose1 = conc1, project_id1 = project_id1, strain_id1 = strain_id1, pert_itime1 = pert_itime1, 
                   cormatrix_path = cormat_path, reference_id_path = ref1_path, sample_metadata_path = cdesc, 
                   remove_query_id_from_ref = TRUE, remove_query_pert_id_from_ref = FALSE) #Doesn't include more because pa69180 was a test compound and not included in the reference set

page0 = corReport_with_negcon(comp1 = comp1, conc1 = conc1, project_id1 = project_id1, strain_id1 = strain_id1, pert_itime1 = pert_itime1, pert_iname1 = pert_iname1, 
                              target1 = target1, description1 = description1, cormatrix_path = cormat_path, compound_annotation_path = NA, reference_id_path = ref1_path, sample_metadata_path = cdesc,
                              remove_query_id_from_ref = TRUE, remove_query_pert_id_from_ref = FALSE, compoundnames = T, plotCorThreshold = compcorthreshold1,
                              target_column = "pert_target", target_column2 = "pert_mechanism", project_ids = unique(cdesc$project_id), 
                              strain_ids = unique(cdesc$strain_id), dose_list = sort(unique(cdesc$dose_rank)), time_list = sort(as.numeric(unique(cdesc$pert_time))), 
                              plot_shape_col = plot_shape_col, plot_size_col = plot_size_col, plot_fill_col = plot_fill_col)

grid.newpage()
grid.draw(page0)


##############Query external compounds#########################################

col_meta = cdesc
savefilefolder = paste0(outdir, 'zscore_correlation_analysis_nzscore_vs_ref_zscore_moc1430_1779_0066_0110_2kd_', date0, '/')
if(!dir.exists(savefilefolder)){
  dir.create(savefilefolder, showWarnings = T)
}else{
  print("directory already exists")
}

query_id_df = filter(cdesc, project_id == "prjna701395" & pert_id != "dmso")
pert_list = unique(query_id_df$pert_id)

#Since there's only a few compounds, I'll put them all in one file
filename1 = paste0(savefilefolder,'psa_rnaseq_correlation_', gsub('[[:space:]]', replacement = '_', gsub('[[:punct:]]', replacement = '_', screen_name)), ".pdf")
pdf(filename1, width = 12, height = 8)


for(pert_num in 1:length(pert_list)){
  pert_id1 = pert_list[pert_num]
  treatments = filter(col_meta, pert_id == pert_id1)
  treatments = treatments %>%
    dplyr::arrange(pert_time, pert_dose, strain_id)
  
  #For some reason the arrangement isn't working so do it manually:
  pert_time_list = sort(unique(treatments$pert_time))
  
  #loop through each time
  for(timenum in 1:length(pert_time_list)){
    pert_time1 = pert_time_list[timenum]
    treatments_filt = filter(treatments, pert_time == pert_time1)
    treatments_filt = arrange(treatments_filt, pert_dose, strain_id)
    dose_list = unique(treatments_filt$pert_idose)
    
      #loop through each dose
      for(dosenum in 1:length(dose_list)){
        conc1 = dose_list[dosenum]
        treatments_dose_filt = filter(treatments_filt, pert_idose == conc1) %>% arrange(strain_id)
        strain_list = unique(treatments_dose_filt$strain_id)
        
        #loop through each strain
        for(strain_num in 1:length(strain_list)){
          strain_id1 = strain_list[strain_num]
          treatments_strain_filt = filter(treatments_dose_filt, strain_id == strain_id1)
          pert_iname1 = treatments_strain_filt$pert_iname
          genetarget1 = ""
          targetdescription1 = ""
          
          
          # genedescription1 = paste(genetarget1, targetdescription1, sep = ":")
          page0 = corReport_with_negcon(comp1 = pert_id1, conc1 = conc1, project_id1 = project_id0, strain_id1 = strain_id1, pert_itime1 = pert_itime1, pert_iname1 = pert_iname1, 
                                        target1 = genetarget1, description1 = targetdescription1, cormatrix_path = cormat_path, compound_annotation_path = NA, reference_id_path = ref1_path, sample_metadata_path = col_meta,
                                        remove_query_id_from_ref = TRUE, remove_query_pert_id_from_ref = FALSE, compoundnames = T, plotCorThreshold = compcorthreshold1,
                                        target_column = "pert_target", target_column2 = "pert_mechanism", project_ids = unique(col_meta$project_id), 
                                        strain_ids = unique(col_meta$strain_id), dose_list = sort(unique(col_meta$dose_rank)), time_list = sort(as.numeric(unique(col_meta$pert_time))), 
                                        plot_shape_col = plot_shape_col, plot_size_col = plot_size_col, plot_fill_col = plot_fill_col)
          
          grid.draw(page0)
          #if(treatments_strain_filt$inactive == "TRUE"){
          #  grid.text("weak signal", x = 0.95, y = 0.95, draw = T)
          #}
          if((strain_num < length(strain_list)) | (dosenum < length(dose_list)) | (timenum < length(pert_time_list)) | (pert_num < length(pert_list))){
            grid.newpage()
          }
        }
      }
    }
  }
dev.off()



