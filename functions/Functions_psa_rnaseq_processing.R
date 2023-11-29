#Functions for RNAseq data
#11/92/21

library(dplyr)
library(ggplot2)

##########Raw data from MOC###########
extract_cds = function(countpath, savefilepath){
  #countpath is the path to the counts file
  #savefilepath, don't add suffix, saves an RDS of the coding genes (CDS) only & its associated metadata
  #metadata releveled for DESeq2
  if(is.character(countpath)){
    if(grepl(".tsv", tolower(countpath))){
      counts = read.table(countpath, sep = '\t', header = T, check.names = F) #check names being false maintains the same characters as input originally
    }else if(grepl(".csv", tolower(countpath))){
      counts = read.csv(countpath, header = T, check.names = F, stringsAsFactors = F)
    }else if(grepl(".rds", tolower(countpath))){
      counts = readRDS(countpath)
    }else{
      stop("Error: input file type not recognized.")
    }
    
  }else{
    #data already loaded
    counts = countpath
  }
  
  # Find coding sequences
  colnames(counts) = tolower(colnames(counts))
  if(any(colnames(counts) %in% c("geneid", "gene_id"))){
    idx = which(colnames(counts) %in% c("geneid", "gene_id"))
    colnames(counts)[idx] = "geneid"
    counts_cds = filter(counts, startsWith(geneid, "CDS"))
    rownames(counts_cds) = counts_cds$geneid
    counts_cds = counts_cds[,-1]
    saveRDS(counts_cds, paste0(savefilepath, '_counts_cds.rds'))
  }else{
    stop("Error: no gene ID column name found")
  }

}

extract_cds_and_metadata = function(countpath, savefilepath, reference_level = "dmsoA_30"){
  #countpath is the path to the counts file
  #savefilepath, don't add suffix, saves an RDS of the coding genes (CDS) only & its associated metadata
  #metadata releveled for DESeq2
  
  counts = read.table(countpath, sep = '\t', header = T)
  # Find coding sequences
  counts_cds = filter(counts, startsWith(Geneid, "CDS"))
  rownames(counts_cds) = counts_cds$Geneid
  counts_cds = counts_cds[,-1]
  saveRDS(counts_cds, paste0(savefilepath, '_counts_cds.rds'))
  
  # head(counts_cds)
  #This is specific to how Keith/Thulasi set up their experiment labels
  treatment_time = colnames(counts_cds)
  metadata = data.frame(sample_id = treatment_time)
  metadata = separate(metadata, col = sample_id, sep = "_", remove = F, into = c("time", "compound", "rep"))
  metadata$time = sub("T", replacement = "", x= metadata$time)
  metadata$time = as.numeric(metadata$time)
  metadata = mutate(metadata, compound_time = paste(compound, time, sep = '_'))
  rownames(metadata) = metadata$sample_id
  metadata$compound_time = as.factor(metadata$compound_time)
  metadata$compound_time = relevel(metadata$compound_time, ref = reference_level)
  saveRDS(metadata, paste0(savefilepath, '_column_metadata.rds'))
  
}

gene_stats = function(count_mat, col_meta, gene_name){
  #calculate IQR and sd for centered values (mean centered to 0 for each condition), and uncentered (across all conditions)
  #count_mat is the count matrix
  #col_meta is the column metadata that includes strain_pert as the condition, and replicate numbers
  #gene_name is any subset of the gene name (row_id)
  #return data frame
  
  gene_idx = which(grepl(gene_name, rownames(count_mat), fixed = T))
  counts_filt = data.frame(count = count_mat[gene_idx,])
  counts_filt$sample_id = rownames(counts_filt)
  counts_filt = left_join(counts_filt, col_meta, by = "sample_id")
  
  #calculate variance, std dev, and IQR
  counts_filt = counts_filt %>%
    group_by(strain_pert) %>%
    mutate(mean_count = mean(count, na.rm = T), median_count = median(count, na.rm = T)) %>%
    ungroup() %>%
    mutate(count_centered = count - mean_count)
  
  df = data.frame(row_id = rownames(count_mat)[gene_idx], max_count = max(counts_filt$count, na.rm = T), mean_count = mean(counts_filt$count, na.rm = T), median_count = median(counts_filt$count, na.rm = T), sd_count_centered = sd(counts_filt$count_centered, na.rm = T), sd_count = sd(counts_filt$count, na.rm = T), iqr_count_centered = IQR(counts_filt$count_centered, na.rm = T), iqr_count = IQR(counts_filt$count, na.rm = T), range_count = max(counts_filt$count, na.rm = T) - min(counts_filt$count, na.rm = T)) 
  df = mutate(df, ratio_sd_count = sd_count/sd_count_centered, ratio_iqr_count = iqr_count/iqr_count_centered)
  return(df)
}

convert_gene_id = function(rds_file_path, gene_meta_path){
  #For data that has old gene IDs (ORF IDS: e.g. PA14_RS09865, convert to PA14_number)
  #Input: rds_file_path is path to file with "gene_id" that are actually orf_ids
  #gene_metadata_path is path to file with gene_ids and orf_ids
  #Output: dataframe with new gene_ids
  #2/21/23
  
  df = readRDS(rds_file_path)
  gene_anno = readRDS(gene_meta_path)
  gene_anno_sub = select(gene_anno, gene_id, orf_id) %>% distinct()
  
  if(any(df$gene_id %in% gene_anno_sub$orf_id)){
    colnames(df)[colnames(df) == "gene_id"] = "orf_id"
    df = left_join(df, gene_anno_sub, by = "orf_id")
  }else if(any(df$gene_id %in% gene_anno_sub$gene_id)){
    print("Gene IDs already converted.")
  }else{
    stop("Error: cannot convert gene names, unidentified format.")
  }
  
  return(df)
  
}
