#General functions for creating directories, saving gct(x) files and combining them.
#jbagnall@broadinstitute.org
#Nov. 30, 2023

library(cmapR)
library(dplyr)

create_dir = function(path_directory){
  #Checks if the directory exists, if not creates it.
  #It also makes sure there is a '/' at the end of the path string
  
  end_char = substr(path_directory, nchar(path_directory), nchar(path_directory))
  if(end_char != "/"){
    path_directory = paste0(path_directory, '/')
  }
  
  if(!dir.exists(path_directory)){
    print(paste("Creating directory:", path_directory))
    dir.create(path_directory)
  }else{
    print("Directory already exists")
  }
}

save_to_gctx = function(data_mat, data_cdesc = NA, data_rdesc = NA, savefilepath = NA){
  #Write gctx file
  
  if(is.na(savefilepath)){
    stop("Error: No file path given to save the file")
  }
  
  data_gct = new("GCT")
  data_gct@rid = rownames(data_mat)
  data_gct@cid = colnames(data_mat)
  data_gct@mat = data_mat
  
  if(!all(is.na(data_rdesc))){
    data_gct@rdesc = data_rdesc
  }
  
  if(!all(is.na(data_cdesc))){
    data_gct@cdesc = data_cdesc
  }
  
  write_gctx(ds = data_gct, ofile = savefilepath)
}

save_to_gct = function(data_mat, data_cdesc = NA, data_rdesc = NA, savefilepath = NA){
  #Write gct file (useful for visualizing in morpheus)
  
  if(is.na(savefilepath)){
    stop("Error: No file path given to save the file")
  }
  
  data_gct <- new("GCT", mat = data_mat)
  # data_gct@rid = rownames(data_mat)
  # data_gct@cid = colnames(data_mat)
  # data_gct@mat = data_mat
  
  if(!all(is.na(data_rdesc))){
    data_gct@rdesc = data_rdesc
  }
  
  if(!all(is.na(data_cdesc))){
    data_gct@cdesc = data_cdesc
  }
  
  write_gct(ds = data_gct, ofile = savefilepath)
}

combine_gct = function(ref_data_path, query_data_path, save_file_path, save_format_gct = TRUE){
  #Merges 2 gct files, e.g. query data + reference data (in order to do target prediction)
  #Input path to reference data (gct), usually moderated zscores
  #Input path to query data (gct), usually moderated zscores or nzscores
  #Input save file path and stem of name
  #If save_format_gct is TRUE, saves as gct, otherwise gctx
  #Outputs a gctx correlation matrix (all by all, symmetric) saved to save_file_path (gctx)
  
  #Load reference data
  ref_data = parse_gctx(ref_data_path)
  ref_data_mat = ref_data@mat
  ref_data_cdesc = ref_data@cdesc
  ref_data_rdesc = ref_data@rdesc
  
  #Load query data
  query_data = parse_gctx(query_data_path)
  query_data_mat = query_data@mat
  query_data_cdesc = query_data@cdesc
  query_data_rdesc = query_data@rdesc
  
  #Intersect the rows (genes)
  genes_intersect = intersect(rownames(ref_data_mat), rownames(query_data_mat))
  
  #Combine data sets
  ref_data_mat_sub = ref_data_mat[genes_intersect,]
  
  query_data_mat_sub = query_data_mat[genes_intersect,]
  
  if(!all(rownames(ref_data_mat_sub) == rownames(query_data_mat_sub))){
    stop("Rownames do not match. Cannot concatenate matrices")
  }
  
  if(all(ref_data_mat_sub == query_data_mat_sub)){
    print("Warning: query and reference are the same. Returning reference matrix")
    combined_mat = ref_data_mat_sub
  }else if(any(colnames(query_data_mat_sub) %in% colnames(ref_data_mat_sub))){
    print("Warning: some column names in query were the same as reference. Resulting column names may be modified.")
    combined_mat = base::cbind(ref_data_mat_sub, query_data_mat_sub)
  }else{
    combined_mat = base::cbind(ref_data_mat_sub, query_data_mat_sub)
  }
  num_row = dim(combined_mat)[1]
  num_col = dim(combined_mat)[2]
  
  #Combine column metadata
  if(!all(colnames(ref_data_cdesc) == colnames(query_data_cdesc))){
    print("Warning: cdesc column names do not match, but will combine using plyr::rbind.fill. Expect NAs.")
  }
  combined_cdesc = plyr::rbind.fill(ref_data_cdesc, query_data_cdesc)
  
  #Add metadata for global dose rank
  if("dose_rank" %in% colnames(combined_cdesc)){
    combined_cdesc = combined_cdesc %>%
      group_by(pert_id) %>%
      arrange(pert_dose) %>%
      mutate(dose_rank_global = dense_rank(pert_dose)) %>%
      ungroup()
    "Calculated dose_rank_global and added to column metadata."
  }else{
    print("Dose_rank is not in cdesc. No global dose rank calculated.")
  }
  
  combined_cdesc = as.data.frame(combined_cdesc)
  
  #Combine row metadata
  if(!all(colnames(ref_data_rdesc) == colnames(query_data_rdesc))){
    print("Warning: rdesc column names do not match, only keeping rdesc from ref_data")
  }

  combined_rdesc = data.frame(gene_id = genes_intersect)
  combined_rdesc = dplyr::left_join(combined_rdesc, ref_data_rdesc, by = "gene_id")
  
  #Reorder to be same as column names in matrix
  match_idx = match(colnames(combined_mat), combined_cdesc$id)
  combined_cdesc = combined_cdesc[match_idx,]
  
  if(save_format_gct){
    save_to_gct(data_mat = combined_mat, data_cdesc = combined_cdesc, data_rdesc = combined_rdesc, savefilepath = save_file_path)
    return(paste0(save_file_path, "_n",num_col, 'x', num_row, '.gct'))
  }else{
    save_to_gctx(data_mat = combined_mat, data_cdesc = combined_cdesc, data_rdesc = combined_rdesc, savefilepath = save_file_path)
    return(paste0(save_file_path, "_n", num_col, 'x', num_row, '.gctx'))
  }

  
}

save_cor_matrix = function(input_data_path, save_file_path, save_format_gct = FALSE, cor_method = "pearson"){
  #Calculate correlation matrix and save
  #Input GCT file, usually moderated zscores
  #Save_file_path is where to save the file
  #save_format_gct defaults to FALSE to save as gctx.
  #cor_method is the correlation method
  #Output correlation matrix, all by all, symmetric, to save_file_path
  
  data = parse_gctx(input_data_path)
  data_mat = data@mat
  data_cdesc = data@cdesc
  data_rdesc = data@rdesc
  
  cor_mat = cor(data_mat, method = cor_method)
  num_row = dim(cor_mat)[1]
  num_col = dim(cor_mat)[2]
  
  if(save_format_gct){
    save_to_gct(data_mat = cor_mat, data_cdesc = data_cdesc, savefilepath = save_file_path)
    return(paste0(save_file_path, "_n", num_col, 'x', num_row, '.gct'))
  }else{
    save_to_gctx(data_mat = cor_mat, data_cdesc = data_cdesc, savefilepath = save_file_path)
    return(paste0(save_file_path, "_n", num_col, 'x', num_row, '.gctx'))
  }

  
}

combine_and_make_cor_matrix = function(ref_data_path, query_data_path, save_combined_path, save_cor_path, combined_mat_format_gct = TRUE, cor_mat_format_gct = FALSE, cor_method = "pearson"){
  #Wrapper for combining & saving gct files and then saving the correlation matrix
  
  combined_path = combine_gct(ref_data_path = ref_data_path, query_data_path = query_data_path, save_file_path = save_combined_path, save_format_gct = combined_mat_format_gct)
  cormat_path = save_cor_matrix(input_data_path = combined_path, save_file_path = save_cor_path, save_format_gct = cor_mat_format_gct, cor_method = cor_method)
  
}
