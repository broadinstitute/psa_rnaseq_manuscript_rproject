#Functions for target/MoA prediction from RNAseq data
#jbagnall@broadinstitute.org
#Nov. 30, 2023

library(dplyr)
library(tidyr)
library(cmapR)
library(ggplot2)

targetCorrelation = function(cormatrix_path, query_ids = NA, reference_ids = NA, sample_meta, compound_meta, target_colname1 = "pert_target", target_colname2 = "pert_mechanism", savefilename = NA, removequery = F){
  #For reference-based MOA predictions
  #cormatrix_path is the correlation matrix between the query treatments (on columns) and the kabx2 knowns (on rows)
  #query_ids are the ids of the query treatments, if NA then takes all columns
  #reference_ids are the ids of the reference treatments (should match the ids in cormatrix), if NA then takes all rows
  #Go through all treatments of screen that satisfy certain conditions (repcorthreshold, nonWK)
  #Find maxmimally correlated treatment from each target class
  #Includes normalized and non-normalized values of the maximum correlation from each target class to the treatment
  #removeanyWK only feeds into calling active treatments, not into ranking of targets (there, well killers are always included)
  #sample_meta includes sample metadata (id, pert_id, pert_dose, pert_time, strain_id, project_id)
  #compound_meta includes compound metadata (target_colname1, target_colname2, etc.)
  #Return: data frame with all treatments and their maximum correlation values to each target class
  
  #correlation matrix, rows are kabx, and columns are queries
  col_meta = read_gctx_meta(cormatrix_path, dim = "col")
  row_meta = read_gctx_meta(cormatrix_path, dim = "row")
  
  
  if(all(is.na(query_ids))){
    query_ids = col_meta$id
  }
  if(all(is.na(reference_ids))){
    reference_ids = row_meta$id
  }
  
  #Load metadata
  if(is.character(sample_meta)){
    sample_meta = read.csv(sample_meta, stringsAsFactors = F)
  }else{
    #assumes the data frame or table has been passed in
    sample_meta = sample_meta
  }
  
  if(is.character(compound_meta)){
    compound_meta = read.csv(compound_meta, stringsAsFactors = F)
  }else{
    #assumes the data frame or table has been passed in
    compound_meta = compound_meta
  }
  
  cormat = parse_gctx(cormatrix_path, cid = which(col_meta$id %in% query_ids), rid = which(row_meta$id %in% reference_ids)) 
  cormat_df = as.data.frame(cormat@mat)
  cormat_df$kabx_id = rownames(cormat_df)
  cormat_lf = tidyr::gather(cormat_df, key = "query_id", value = "CompoundCorrelation", -kabx_id)

  
  colnames(compound_meta)[colnames(compound_meta)==target_colname1] = "target"
  colnames(compound_meta)[colnames(compound_meta)==target_colname2] = "target2"
  
  #add sample metadata for query
  metadata_select = dplyr::select(sample_meta, id, pert_id, pert_dose, pert_time, strain_id, project_id)
  query_sample_meta = metadata_select
  colnames(query_sample_meta) = paste0("query_", colnames(query_sample_meta))
  cormat_lf_anno = dplyr::left_join(cormat_lf, query_sample_meta, by = "query_id")
  rm(query_sample_meta)
  
  #add sample and compound metadata for poscons
  colnames(metadata_select) = paste0("kabx_", colnames(metadata_select))
  colnames(compound_meta) = paste0("kabx_", colnames(compound_meta))
  cormat_lf_anno = dplyr::left_join(cormat_lf_anno, metadata_select, by = "kabx_id")
  cormat_lf_anno = dplyr::left_join(cormat_lf_anno, compound_meta, by = "kabx_pert_id")
  
  #Double check to remove compounds with no annotation
  #Make sure there are no redundancies, reduce to only target information
  cormat_lf_anno = filter(cormat_lf_anno, kabx_target != "-666" & kabx_target != "NA" & kabx_target != "unknown" & kabx_target != "whole cell only")
  cormat_lf_anno = filter(cormat_lf_anno, !is.na(kabx_target))
  cormat_lf_anno = distinct(cormat_lf_anno)
  
  
  if(removequery){#remove query pert_id, only needed if testing with knowns for cross validation
    cormat_lf_anno = filter(cormat_lf_anno, kabx_pert_id != query_pert_id) 
  }
  
  #pick top correlation for each treatment from each target class
  cormat_lf_treatments = cormat_lf_anno %>%
    group_by(query_id, kabx_target) %>%
    dplyr::slice(which.max(CompoundCorrelation)) %>%
    ungroup() %>%
    group_by(query_id) %>%
    mutate(max_compcor = max(CompoundCorrelation, na.rm = T)) %>%
    ungroup() %>%
    mutate(compcor_norm = CompoundCorrelation/max_compcor)
  
  if(is.na(savefilename)){
    return(cormat_lf_treatments)
  }else{
    saveRDS(cormat_lf_treatments, savefilename)  
    print(paste0("File saved to: ", savefilename))
  }

}

refBasedList_compound = function(target_list_path, rank_thresh = 1, normalized_cor = F, kabx_annopath= "", target_col = "pert_target", roc_path = '', roc_path2 = '', savefilename, avg_doses = T, print_plot = F, save_all_target_cor = T){
  #6/7/22 List prediction for each compound(pert_id)
  #target_list is the output of targetCorrelation, where every query_pert_id has been correlated to a maximum from a target category from reference
  #it can either be a list or a data frame
  #normalized_cor == T, then uses ranks for the mean normalized correlation values (instead of the mean of the actual correlation values)
  #target_col should match the name of target column1 from targetCorrelation (used to define target)
  #annotation file requires pert_type = poscon, negcon or test. defines "unknowns" as "test" compounds
  #if avg_doses = T, then will use the mean of the correlations across doses (this is default). If F, then will use the max across all doses.
  #if save_all_target_cor = T, then will save RDS files of max correlations to all targets, useful for plotting purposes later
  #roc_path should give the positive predictive value of predicting the target_col (usually pert_target), given the correlation threshold. Will ignore if no input.
  #roc_path2 should give the positive predictive value of predicting the mechanism (based directly on the target prediction), given the correlation threshold. Will ignore if no input.
  
  #Outputs target predictions for each compound in a .csv file
  #subsets known (in reference set) vs unknowns (not in reference set)
  #saves 2 csv files, one with all knowns, one with unknowns
  #List all correlation values and the confidence score (positive predictive value, what fraction with that correlation value (or higher) were called correctly in the reference set)
  #Output column explanations in csv file:
  #query_pert_id: query compound
  #predicted_target: predicted target
  #target_correlation: maximal correlation between query compound and the treatments in the predicted target category
  #normalized_target_correlation: target_correlation divided by the maximal correlation between the query treatment and any target category
  #positive_predictive_value: fraction of compounds whose targets would be predicted correctly in the reference set (Leave-one-compound-out-cross_validation), with the listed correlation value or higher
  #neighbor_pert_id: the compound from the reference set with the highest correlation to the query
  #predicted_pert_target: the target of the neighbor_pert_id
  #predicted_pharmaceutical_class: the pharmaceutical class of the neighbor_pert_id
  
  
  if(is.character(kabx_annopath)){
    kabx_annotation = read.csv(kabx_annopath, stringsAsFactors = F)
  }else{
    #assumes the data frame or table has been passed in
    kabx_annotation = kabx_annopath
  }
  
  kabx_anno_collapse = kabx_annotation %>%
    group_by(pert_id) %>%
    summarise(pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism ), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|"), pert_type = paste(unique(pert_type), collapse = "|")) %>%
    ungroup() %>% 
    distinct()
  
  target_df = readRDS(target_list_path)
  if(any(class(target_df) == "list")){
    target_df = do.call(rbind, target_df)
  }
  
  if(any(colnames(target_df) == "pert_iname")){
    colnames(target_df)[colnames(target_df)=="pert_iname"] = "kabx_pert_iname"
  }
  
  target_df = left_join(target_df, kabx_anno_collapse, by = c("query_pert_id" = "pert_id"))
  colnames(target_df)[colnames(target_df) == target_col] = "target_col"
  
  #separate the knowns from the unknowns
  target_knowns = filter(target_df, ((!target_col %in% c("-666", "NA", "unknown", "whole cell only")) & !is.na(target_col) & (pert_type != "test")))
  target_unknowns = filter(target_df, target_col %in% c("-666", "NA", "unknown", "whole cell only") | is.na(target_col) | pert_type == "test")  
  
  #change column names back to originals
  colnames(target_knowns)[colnames(target_knowns) == "target_col"] = target_col
  colnames(target_unknowns)[colnames(target_unknowns) == "target_col"] = target_col
  
  #Find negative control treatments to remove? These turn out to be different than the ones we removed from poscon reference because those we allowed to 
  #correlate to other doses of the same drug, but here we remove the query drug entirely...
  negcon = target_df %>% 
    group_by(query_id, query_project_id, query_pert_id, query_pert_dose, query_pert_time, query_strain_id) %>%
    dplyr::slice(which.max(CompoundCorrelation)) %>%
    ungroup() %>%
    filter(kabx_target %in% c("Negative control", "negative control", "negcon", "Negative Control")) #%>%
  # mutate(query_comp_conc = paste(project_id, query_pert_id, query_pert_dose, query_pert_time, query_strain_id, sep = ":"))
  
  if (dim(target_unknowns)[1]>0){
    
    if(avg_doses){
      unknown_df = target_unknowns %>%
        dplyr::filter(!(query_id %in% unique(negcon$query_id))) %>% #remove treatments that are most correlated to negcons
        group_by(query_pert_id, kabx_target) %>% #each project, pert_id, strain (average over dose and time)
        summarise(mean_compcor = mean(CompoundCorrelation, na.rm = T), mean_compcor_norm = mean(compcor_norm, na.rm = T), pert_id = paste(unique(kabx_pert_id), collapse = "|")) %>%
        ungroup() %>%
        group_by(query_pert_id) %>% #collapse this to pert_id only once confident that project & strain can be disregarded
        arrange(desc(mean_compcor_norm)) %>%
        mutate(target_rank_norm = 1:n()) %>%
        arrange(desc(mean_compcor)) %>%
        mutate(target_rank_notnorm = 1:n()) %>%
        # filter(mean_compcor_norm >= 0.85) %>%
        ungroup() 
    }else{
      unknown_df = target_unknowns %>%
        dplyr::filter(!(query_id %in% unique(negcon$query_id))) %>% #remove treatments that are most correlated to negcons
        group_by(query_pert_id, kabx_target) %>% #each project, pert_id, strain (average over dose and time)
        summarise(mean_compcor = max(CompoundCorrelation, na.rm = T), mean_compcor_norm = max(compcor_norm, na.rm = T), pert_id = paste(unique(kabx_pert_id), collapse = "|")) %>%
        ungroup() %>%
        group_by(query_pert_id) %>% #collapse this to pert_id only once confident that project & strain can be disregarded
        arrange(desc(mean_compcor_norm)) %>%
        mutate(target_rank_norm = 1:n()) %>%
        arrange(desc(mean_compcor)) %>%
        mutate(target_rank_notnorm = 1:n()) %>%
        # filter(mean_compcor_norm >= 0.85) %>%
        ungroup() 
    }
    
    if(save_all_target_cor){
      if(is.na(savefilename) | savefilename == ""){
        print("Warning: Not saving file, filepath is NA")
      }else{
        unknown_df_save = unknown_df
        unknown_df_save = unknown_df_save %>%
          unnest(pert_id = strsplit(pert_id, "|", fixed = T)) %>%
          mutate(neighbor_pert_id = ifelse(grepl("BRD", pert_id), substr(pert_id, 1, 13), pert_id))
        unknown_df_save = left_join(unknown_df_save, kabx_anno_collapse, by = c("neighbor_pert_id" = "pert_id"))
        
        unknown_df_save = unknown_df_save %>%
          group_by(query_pert_id, kabx_target, mean_compcor, mean_compcor_norm) %>%
          summarise(neighbor_pert_id = paste(unique(neighbor_pert_id), collapse = "|"), pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|")) %>%
          ungroup()
        colnames(unknown_df_save) = c("query_pert_id", "predicted_target", "target_correlation", "normalized_target_correlation", "neighbor_pert_id", "neighbor_pert_iname", "predicted_pert_mechanism", "predicted_pert_target", "predicted_pharmaceutical class")
        
        saveRDS(unknown_df_save, paste(savefilename, "_unknown_all_target_cor.rds", sep = ''))
      }
    }

    
    if(normalized_cor){
      unknown_df_output = filter(unknown_df, target_rank_norm <= rank_thresh)
      
      #Remove final negcons (cases where on average, negcons was most correlated, but not for individual treatments)
      negcon1 = unknown_df_output %>%
        filter((target_rank_norm == 1) & (kabx_target %in% c("Negative control", "negative control", "negcon", "Negative Control")))
      
    }else{
      unknown_df_output = filter(unknown_df, target_rank_notnorm <= rank_thresh)
      negcon1 = unknown_df_output %>%
        filter((target_rank_notnorm == 1) & (kabx_target %in% c("Negative control", "negative control", "negcon", "Negative Control")))
    }
    
    #Remove final negcons (cases where on average, negcons was most correlated, but not for individual treatments)
    unknown_df_output = filter(unknown_df_output, !(query_pert_id %in% negcon1$query_pert_id))
    
    
    if(print_plot){
      p0 = ggplot(unknown_df_output, aes(x = mean_compcor_norm, y = mean_compcor))+
        geom_point()+
        labs(x = "Mean Normalized Correlation", y = "Mean Correlation", title = paste("Rank <=", rank_thresh, sep = " "))+
        theme_bw(base_size = 16)
      print(p0)
    }
    
    #Add FDR from reference set
    if(is.na(roc_path)|roc_path == ""){
      unknown_df_output$threshold_idx = NA
      unknown_df_output$fdr = NA
      unknown_df_output$precision = NA
      
    }else{
      kabx2_rocinfo = readRDS(roc_path)
      
      # test = unknown_df_output
      unknown_df_output$threshold_idx = sapply(unknown_df_output$mean_compcor, function(x){max(which(x > kabx2_rocinfo$threshold))})
      unknown_df_output$fdr = kabx2_rocinfo$fdr[unknown_df_output$threshold_idx]
      unknown_df_output$precision = kabx2_rocinfo$precision[unknown_df_output$threshold_idx]
      
    }
    
    #Add more annotations based on neighbor (target) compounds
    #separate out each neighbor compound
    unknown_df_output_anno = unknown_df_output %>%
      unnest(pert_id = strsplit(pert_id, "|", fixed = T)) %>%
      mutate(neighbor_pert_id = ifelse(grepl("BRD", pert_id), substr(pert_id, 1, 13), pert_id))
    unknown_df_output_anno = left_join(unknown_df_output_anno, kabx_anno_collapse, by = c("neighbor_pert_id" = "pert_id"))
    
    
    #collapse neighbors again
    unknown_df_output_anno = select(unknown_df_output_anno, -target_rank_norm, -target_rank_notnorm, -threshold_idx, -fdr)
    
    unknown_df_output_anno = unknown_df_output_anno %>%
      group_by(query_pert_id, kabx_target, mean_compcor, mean_compcor_norm, precision) %>%
      summarise(neighbor_pert_id = paste(unique(neighbor_pert_id), collapse = "|"), pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|")) %>%
      ungroup()
    colnames(unknown_df_output_anno) = c("query_pert_id", "predicted_target", "target_correlation", "normalized_target_correlation", paste("positive_predictive_value", target_col, sep = "_"), "neighbor_pert_id", "neighbor_pert_iname", "predicted_pert_mechanism", "predicted_pert_target", "predicted_pharmaceutical class")
    
    #Add additional PPV if desired (e.g. for mechanism)
    if(!is.na(roc_path2) & roc_path2 != ""){
      kabx2_rocinfo2 = readRDS(roc_path2)
      unknown_df_output_anno$threshold_idx2 = sapply(unknown_df_output_anno$target_correlation, function(x){max(which(x > kabx2_rocinfo2$threshold))})
      unknown_df_output_anno$positive_predictive_value2 = kabx2_rocinfo2$precision[unknown_df_output_anno$threshold_idx2]
      unknown_df_output_anno = dplyr::select(unknown_df_output_anno, -threshold_idx2)
    }
    
    unknown_df_output_anno = arrange(unknown_df_output_anno, desc(target_correlation))
    
    write.csv(unknown_df_output_anno, paste(savefilename, "_unknown.csv", sep = ''), row.names = F)

    #Distribution of predicted targets
    if(print_plot){
      target_summary = unknown_df_output_anno %>%
        group_by(predicted_pert_mechanism) %>%
        summarise(numcomps = n()) %>%
        ungroup()
      
      p1 = ggplot(target_summary, aes(x = reorder(predicted_pert_mechanism, numcomps), y = numcomps))+
        geom_bar(stat = "identity")+
        labs(x = "Predicted mechanism", y = "# Compounds") +
        theme_bw(base_size = 16)+
        theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
      print(p1)
      
      # if(tani_path != ''){
      #   p2 = ggplot(unknown_df_output_anno, aes(x = target_correlation, y = max_tani))+
      #     geom_point()+
      #     labs(x = "Average correlation", y = "Max Tanimoto", title = "Unknowns")+
      #     xlim(0,1)+
      #     ylim(0,1)+
      #     theme_bw(base_size = 16)
      #   print(p2)
      # }
    }
  }
  
  ####annotate and predict on knowns###
  #predict MOA for unknowns
  #Note that the knowns were not fairly predicted because the same pert_id may not have been allowed to be included in the correlation search?
  if (dim(target_knowns)[1]>0){
    if(avg_doses){
      known_df = target_knowns %>%
        filter(!(query_id %in% unique(negcon$query_id))) %>%
        group_by(query_pert_id, kabx_target) %>%
        summarise(mean_compcor = mean(CompoundCorrelation, na.rm = T), mean_compcor_norm = mean(compcor_norm, na.rm = T), pert_id = paste(unique(kabx_pert_id), collapse = "|")) %>%
        ungroup() %>%
        group_by(query_pert_id) %>%
        arrange(desc(mean_compcor_norm)) %>%
        mutate(target_rank_norm = 1:n()) %>%
        arrange(desc(mean_compcor)) %>%
        mutate(target_rank_notnorm = 1:n()) %>%
        # filter(mean_compcor_norm >= 0.85) %>%
        ungroup() 
    }else{
      #use maximum across doses
      known_df = target_knowns %>%
        filter(!(query_id %in% unique(negcon$query_id))) %>%
        group_by(query_pert_id, kabx_target) %>%
        summarise(mean_compcor = max(CompoundCorrelation, na.rm = T), mean_compcor_norm = max(compcor_norm, na.rm = T), pert_id = paste(unique(kabx_pert_id), collapse = "|")) %>%
        ungroup() %>%
        group_by(query_pert_id) %>%
        arrange(desc(mean_compcor_norm)) %>%
        mutate(target_rank_norm = 1:n()) %>%
        arrange(desc(mean_compcor)) %>%
        mutate(target_rank_notnorm = 1:n()) %>%
        # filter(mean_compcor_norm >= 0.85) %>%
        ungroup() 
    }

    if(save_all_target_cor){
      if(is.na(savefilename) | savefilename == ""){
        print("Warning: Not saving file, filepath is NA")
      }else{
        known_df_save = known_df
        known_df_save = known_df_save %>%
          unnest(pert_id = strsplit(pert_id, "|", fixed = T)) %>%
          mutate(neighbor_pert_id = ifelse(grepl("BRD", pert_id), substr(pert_id, 1, 13), pert_id))
        known_df_save = left_join(known_df_save, kabx_anno_collapse, by = c("neighbor_pert_id" = "pert_id"))
        
        known_df_save = known_df_save %>%
          group_by(query_pert_id, kabx_target, mean_compcor, mean_compcor_norm) %>%
          summarise(neighbor_pert_id = paste(unique(neighbor_pert_id), collapse = "|"), pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|")) %>%
          ungroup()
        colnames(known_df_save) = c("query_pert_id", "predicted_target", "target_correlation", "normalized_target_correlation", "neighbor_pert_id", "neighbor_pert_iname", "predicted_pert_mechanism", "predicted_pert_target", "predicted_pharmaceutical class")
        
        saveRDS(known_df_save, paste(savefilename, "_known_all_target_cor.rds", sep = ''))
      }
    }
    
    
    if(normalized_cor){
      known_df_output = filter(known_df, target_rank_norm <= rank_thresh)
      #Remove final negcon, those whose average top rank was negcon
      negcon2 = known_df_output %>%
        filter((target_rank_norm == 1) & (kabx_target %in% c("Negative control", "negative control", "negcon", "Negative Control")))
      
    }else{
      known_df_output = filter(known_df, target_rank_notnorm <= rank_thresh)
      negcon2 = known_df_output %>%
        filter((target_rank_notnorm == 1) & (kabx_target %in% c("Negative control", "negative control", "negcon", "Negative Control")))
    }
    
    #Remove final negcon, those whose average top rank was negcon
    known_df_output = filter(known_df_output, !(query_pert_id %in% negcon2$query_pert_id))
    
    if(is.na(roc_path)| roc_path == ""){
      known_df_output$threshold_idx = NA
      known_df_output$fdr = NA
      known_df_output$precision = NA
      
    }else{
      kabx2_rocinfo = readRDS(roc_path)
      known_df_output$threshold_idx = sapply(known_df_output$mean_compcor, function(x){max(which(x > kabx2_rocinfo$threshold))})
      known_df_output$fdr = kabx2_rocinfo$fdr[known_df_output$threshold_idx]
      known_df_output$precision = kabx2_rocinfo$precision[known_df_output$threshold_idx]
      
    }
    
    #Add annotations to the knowns
    kabx_anno_collapse_knowns = kabx_anno_collapse
    colnames(kabx_anno_collapse_knowns) = paste0("query_", colnames(kabx_anno_collapse_knowns))
    known_df_output_anno = left_join(known_df_output, kabx_anno_collapse_knowns, by = "query_pert_id")
    
    
    #Annotate predictions
    #separate out each neighbor compound
    known_df_output_anno = known_df_output_anno %>%
      unnest(pert_id = strsplit(pert_id, "|", fixed = T)) %>%
      mutate(neighbor_pert_id = ifelse(grepl("BRD", pert_id), substr(pert_id, 1, 13), pert_id))
    known_df_output_anno = left_join(known_df_output_anno, kabx_anno_collapse, by = c("neighbor_pert_id" = "pert_id"))
    
    #Add tanimoto information if given
    # if(tani_path != ''){
    #   row_meta = read.gctx.meta(tani_path, dimension = "row")
    #   col_meta = read.gctx.meta(tani_path, dimension = "col")
    #   
    #   known_df_output_anno = known_df_output_anno %>%
    #     group_by(query_pert_id, neighbor_pert_id) %>%
    #     mutate(tanimoto = ifelse((query_pert_id %in% row_meta$id) & (neighbor_pert_id %in% col_meta$id), as.numeric(parse.gctx(tani_path, rid = which(row_meta$id == query_pert_id), cid = which(col_meta$id == neighbor_pert_id))@mat), NaN)) %>%
    #     ungroup()
    # }else{
    #   known_df_output_anno$tanimoto = NaN
    # }
    
    #collapse neighbors again
    known_df_output_anno = select(known_df_output_anno, -target_rank_norm, -target_rank_notnorm, -threshold_idx, -fdr)
    
    known_df_output_anno = known_df_output_anno %>%
      group_by(query_pert_id, kabx_target, mean_compcor, mean_compcor_norm, precision, query_pert_iname, query_pert_mechanism, query_pert_target, query_pharmaceutical_class) %>%
      summarise(neighbor_pert_id = paste(unique(neighbor_pert_id), collapse = "|"), pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|")) %>%
      ungroup()
    
    colnames(known_df_output_anno) = c("query_pert_id", "predicted_target", "target_correlation", "normalized_target_correlation", paste("positive_predictive_value", target_col, sep = "_"), "query_pert_iname", "query_pert_mechanism", "query_pert_target", "query_pharmaceutical_class", "neighbor_pert_id", "neighbor_pert_iname", "predicted_pert_mechanism", "predicted_pert_target", "predicted_pharmaceutical_class")
    
    #Add additional PPV if desired (e.g. for mechanism)
    if(!is.na(roc_path2) & roc_path2 != ""){
      kabx2_rocinfo2 = readRDS(roc_path2)
      known_df_output_anno$threshold_idx2 = sapply(known_df_output_anno$target_correlation, function(x){max(which(x > kabx2_rocinfo2$threshold))})
      known_df_output_anno$positive_predictive_value2 = kabx2_rocinfo2$precision[known_df_output_anno$threshold_idx2]
      known_df_output_anno = dplyr::select(known_df_output_anno, -threshold_idx2)
    }
    
    known_df_output_anno = arrange(known_df_output_anno, desc(target_correlation))
    
    colnames(known_df_output_anno)[colnames(known_df_output_anno) == paste0("query_", target_col)] = "query_target"
    known_df_output_anno = known_df_output_anno %>%
      group_by(query_pert_id, predicted_target, neighbor_pert_id) %>%
      mutate(same_target = grepl(predicted_target, query_target, fixed = T)) %>%
      ungroup
    colnames(known_df_output_anno[colnames(known_df_output_anno) == "query_target"]) = paste0("query_", target_col)
    write.csv(known_df_output_anno, paste(savefilename, "_known.csv", sep = ""), row.names = F)
  }
}

plot_refBasedList = function(pert_id0, all_max_cor_df, kabx_anno_path, text_size = 14, save_file_path = NA){
  #Inputs results from MOA prediction from refBasedList_compound
  #pert_id is the query compound 
  #all_max_cor_df is a data frame output when running refBasedList_compound, which has the maximal correlation from each target category to the query (after collapsing doses etc.)
  #kabx_anno_path has pert_id and pert_iname descriptions for the compounds
  #Outputs a plot of all correlations to the query from high to low
  #if save_file_path is not NA, then saves first plot to pdf
  
  plot_list = list()
  
  max_cor_df = readRDS(all_max_cor_df)
  max_cor_df_filt = dplyr::filter(max_cor_df, query_pert_id == pert_id0)

  max_cor_df_filt_max = max_cor_df_filt %>% dplyr::slice(which.max(target_correlation))
  closest_compound = unique(max_cor_df_filt_max$neighbor_pert_iname)
  
  if(is.character(kabx_anno_path)){
    kabx_annotation = read.csv(kabx_anno_path, stringsAsFactors = F)
  }else{
    #assumes the data frame or table has been passed in
    kabx_annotation = kabx_anno_path
  }
  
  kabx_anno_collapse = kabx_annotation %>%
    group_by(pert_id) %>%
    summarise(pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism ), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|"), pert_type = paste(unique(pert_type), collapse = "|")) %>%
    ungroup() %>% 
    distinct()
  
  if (!("query_pert_iname" %in% colnames(max_cor_df))){
    kabx_query = kabx_anno_collapse
    colnames(kabx_query) = paste0("query_", colnames(kabx_query))
    max_cor_df_filt = dplyr::left_join(max_cor_df_filt, kabx_query, by = c("query_pert_id"))
  }
  
  #changes pert_iname to pert_id if it's NA
  max_cor_df_filt = dplyr::mutate(max_cor_df_filt, query_pert_iname = ifelse(is.na(query_pert_iname), query_pert_id, query_pert_iname))
  
  pert_iname0 = unique(max_cor_df_filt$query_pert_iname)
  max_cor_df_filt = arrange(max_cor_df_filt, desc(target_correlation))

  #color code the x labels by mechanism
  mech_color_df = data.frame(pert_mechanism = c("Negative control", "DNA synthesis", "Membrane Integrity", "Peptidoglycan synthesis", "Protein synthesis"), mech_color = c("black", "red", "blue","green4","magenta3"))
  target_mech_list = select(max_cor_df_filt, predicted_target, predicted_pert_mechanism) %>% distinct()
  target_mech_list = left_join(target_mech_list, mech_color_df, by = c("predicted_pert_mechanism" = "pert_mechanism"))
  color_vec = target_mech_list$mech_color
  
  #make a factor for the mechanism
  rank_mech = max_cor_df_filt %>%
    group_by(predicted_pert_mechanism) %>%
    dplyr::slice(which.max(target_correlation)) %>%
    ungroup %>%
    arrange(desc(target_correlation)) %>%
    mutate(rank = 1:n()) %>%
    ungroup()
  
  rank_vec = rank_mech$predicted_pert_mechanism
  names(rank_vec) = rank_mech$rank
  
  max_cor_df_filt$predicted_pert_mechanism_factor = factor(max_cor_df_filt$predicted_pert_mechanism, levels = rank_vec)
  
  # p0 = ggplot(max_cor_df_filt, aes(x = reorder(predicted_target, -target_correlation), y = target_correlation))+
  #   geom_point(size = 2)+
  #   geom_text(label = paste0("Closest compound: \n", closest_compound), x = 18, y = 0.95, size = 4.5, font_face = "plain")+
  #   geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray")+
  #   labs(x = "Target", y = "Pearson Correlation", title = pert_iname0)+
  #   ylim(-0.2,1)+
  #   theme_bw(base_size = 14)+
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = color_vec))+
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))
  

  
  #Facet by mechanism
  p0 = ggplot(max_cor_df_filt, aes(x = reorder(predicted_target, -target_correlation), y = target_correlation, color = predicted_pert_mechanism))+
    geom_point(size = 2)+
    #geom_text(label = paste0("Closest compound: \n", closest_compound), x = 18, y = 0.95, size = 4.5, font_face = "plain")+
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray")+
    labs(x = "Target", y = "Similarity Score", title = paste0(pert_iname0, " (Closest drug: ", closest_compound, ")"), color = "Mechanism")+
    facet_grid(~predicted_pert_mechanism_factor, space="free_x", scales = "free_x")+
    ylim(-0.2,1)+
    scale_color_manual(values = c("Negative control" = "black", "DNA synthesis" = "red", "Membrane Integrity" = "blue", "Peptidoglycan synthesis" = "green4", "Protein synthesis" = "magenta3"))+
    theme_bw(base_size = text_size)+
    theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = element_blank())
  
  # print(p0)
  
  
  #Color the points by Mechanism
  p1 = ggplot(max_cor_df_filt, aes(x = reorder(predicted_target, -target_correlation), y = target_correlation, color = predicted_pert_mechanism))+
    geom_point(size = 2)+
    #geom_text(label = paste0("Closest compound: \n", closest_compound), x = 18, y = 0.95, size = 4.5, font_face = "plain")+
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray")+
    labs(x = "Target", y = "Similarity Score", title = pert_iname0, color = "Mechanism")+
    ylim(-0.2,1)+
    scale_color_manual(values = c("Negative control" = "black", "DNA synthesis" = "red", "Membrane Integrity" = "blue", "Peptidoglycan synthesis" = "green4", "Protein synthesis" = "magenta3"))+
    theme_bw(base_size = text_size)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = color_vec))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # print(p1)
  
  plot_list = list(p0, p1)
  
  if(!is.na(save_file_path)){
    pdf(save_file_path, width = 7, height = 5)
    print(p0)
    dev.off()
  }else{
    return(plot_list)
  }
  
}
