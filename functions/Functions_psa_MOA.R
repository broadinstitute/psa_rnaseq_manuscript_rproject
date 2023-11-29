library(dplyr)
library(tidyr)
library(cmapR)
library(ggplot2)
library(grid)
library(gridExtra)

compoundCorrelationdf <- function(comp1, conc1, pert_itime1 = "0min", project_id1, strain_id1, cormatrix_path, cormatrix_path2 = NA, comptable_path, comptable_path2 = NA, tanimoto_path = NA, repCorThreshold = 0.4){
  #Input compound and concentration
  #Input either the compoundtable path or the data frame
  #Input correlation matrix path gctx (has query compounds on columns and kabx2 columns on rows)
  #returns correlation between query comp_conc and other comp_conc in the screen (comptable)
  #comp1 = pert_id, conc1 = pert_idose
  #output data frame  of compound correlations to the query comp_conc, listed in order from high to low
  #For reference based predictions, correlation matrix: unknowns on columns & knowns on rows
  #for refernce-free predictions, correlation matrix should be square: unknowns on columns & on rows
  #comptablepath is the compound table for the unknown screen
  #1/29/20 add comptable_path2 is optional, and is the compound table for the kabx screen, if it is different from comptable_path are different.
  #tanimoto only works if both the query and the kabx are in the same tanimoto matrices
  #1/30/20 add cormatrix_path2 which should have the query compound correlated to its other concentrations in the same screen. It is optional, since this may be included in the first cormatrixpath
  #Warning: If cormatrix_path2 is not present, it will try assign compoundcorrelation = 1 for the query. If query is already present in cormatrixpath1, it will assume it is the correct one (from the same screen)
  #Need to be careful if a compound from kabx2 was screened again in another screen--will extract information from cormatrix_path1, if cormatrix_path2 is NA, which if the treatment does exist in 2 screens, would be erroneous
  #timepoint1 = time point (pert_itime), now comp_conc is comp_conc_time
  
  #doesn't include the query compound's correlation to it's other doses if it's not in the columns...do we want to add this?
  
  query_compconc = paste(project_id1, comp1, conc1, pert_itime1, strain_id1, sep = ":") #comp1 = pert_id, conc1 = pert_idose
  
  cormatrix_rowmeta = read_gctx_meta(cormatrix_path, dim = "row")
  cormatrix_colmeta = read_gctx_meta(cormatrix_path, dim = "col")
  
  # #Refine compound table based on thresholds
  if(any(class(comptable_path) =="character")){
    comptable = readRDS(comptable_path)
  } else{
    #else assume that it's a table or data frame
    comptable = comptable_path
  }
  
   #check if column comp_conc exists, added 5/11/22
  comptable = mutate(comptable, comp_conc = paste(project_id, pert_id, pert_idose, pert_itime, strain_id, sep = ":"))
 
  #find query id--hopefully there's only 1
  query_id = unique(filter(comptable, comp_conc == query_compconc)$id)
  if(length(query_id) > 1){
    # query_id = query_id[1]
    print("There were multiple id's associated with this condition.")
  }
  
  keepcomps = comptable


  if(!(any(grepl(comp1, keepcomps$pert_id, fixed = T)))){
    print(paste("Query compound not found", sep = ""))
    # temp = filter(comptable, grepl(comp1, compound, fixed = T))
    # temp <- temp %>%
    #   slice(which.max(correlation))
    # keepcomps = rbind(keepcomps, temp)
    # print("Chose concentration with maximum replicate correlation")
  }
  
  #Look at all the concentrations for the given compound
  keepcomps_filtcomp = filter(keepcomps, grepl(comp1, pert_id, fixed = T))
  keepcomps_filtcomp = mutate(keepcomps_filtcomp, comp_conc = paste(project_id, pert_id, pert_idose, pert_itime, strain_id, sep =':'))
  # print(paste(keepcomps_filtcomp$pert_id, keepcomps_filtcomp$pert_idose, keepcomps_filtcomp$project_id, collapse = ", "))
  print(paste(keepcomps_filtcomp$id, collapse = ", "))
  
  #filter by compound and concentration
  if(!(conc1 %in% keepcomps_filtcomp$pert_idose)){
    print(paste("Error: concentration", conc1, "not found"))
    stop()
  }
  if(!(pert_itime1 %in% keepcomps_filtcomp$pert_itime)){
    print(paste("Error: timepoint", pert_itime1, "not found"))
    stop()
  }
  if(!(project_id1 %in% keepcomps_filtcomp$project_id)){
    print(paste("Error: project_id", project_id1, "not found"))
    stop()
  }
  
  #find exact query (compound, concentration, timepoint, strain_id)
  keepcomps_filt = filter(keepcomps, grepl(comp1, pert_id, fixed = T) & pert_idose == conc1 & pert_itime == pert_itime1 & project_id == project_id1 & strain_id == strain_id1)
  
  if(!is.na(comptable_path2)){
    remove(keepcomps)
    if(any(class(comptable_path2) =="character")){
      comptable2 = readRDS(comptable_path2)
    } else{
      #else assume that it's a table or data frame
      comptable2 = comptable_path2
    }

    
    comptable2 <- mutate(comptable2, comp_conc = paste(project_id, pert_id, pert_idose, pert_itime, strain_id, sep = ":"))
    comptable2 <- mutate(comptable2, pert_id = ifelse(grepl("BRD", pert_id), substr(pert_id, 1, 13), pert_id))
    
    comptable2$pert_dose = as.numeric(comptable2$pert_dose)
    keepcomps = comptable2
  }
  
  
  colid = which(cormatrix_colmeta$id %in% keepcomps_filt$id)
  rowid = which(cormatrix_rowmeta$id %in% keepcomps$id) #comes from second comptable_path2 if given
  cor_matrix = parse_gctx(cormatrix_path, rid = rowid, cid = colid)
  cor_matrix = cor_matrix@mat
  
  #return list of most compounds ranked by correlation
  comp_conc1 = colnames(cor_matrix)[grepl(comp1, colnames(cor_matrix), fixed = T)] #only gives the comp_conc because should only have one column from keepcomps_filt
  temp <- as.data.frame(cor_matrix[,comp_conc1])
  colnames(temp)[1]= "CompoundCorrelation" 
  temp$id =  rownames(temp)
  # if(!(query_compconc %in% temp$comp_conc) & !is.na(cormatrix_path2)){ 
  
  #Remove the correlation of 1?
  if(!is.na(cormatrix_path2)){
    #Assumes if the query was not in the first cormatrix, it is in the second, if it is supplied
    row_meta2 = read_gctx_meta(cormatrix_path2, dim = "row")
    col_meta2 = read_gctx_meta(cormatrix_path2, dim = "col")
    #Choose all concentrations of query compound along rows
    cormat2 = parse_gctx(cormatrix_path2, rid = which(grepl(comp1, row_meta2$id, fixed = T)), cid = which(grepl(tolower(query_compconc), tolower(col_meta2$id))))
    cormat2 = as.data.frame(cormat2@mat)
    colnames(cormat2)[1] = "CompoundCorrelation"
    cormat2$id = rownames(cormat2)
    temp = rbind(cormat2, temp) %>% distinct()
  }else if(!(any(grepl(paste(query_id, collapse = "|"),temp$id)))){ #if no second correlation matrix is supplied
    query_df = data.frame(CompoundCorrelation = 1, comp_conc = query_compconc)
    query_df = left_join(query_df, select(comptable, comp_conc, id), by= "comp_conc") #get the id
    temp = rbind(select(query_df, CompoundCorrelation, id), temp)
  }
  
  # else{
  #   #artificially change query correlation to 1 (if it isn't 1, that means the same treatment was screend in two different screens)
  #   temp$CompoundCorrelation[grepl(query_compconc, temp$id)] = 1
  # }
  
  # temp <- separate(temp, id, into = c("project_id", "pert_id", "pert_idose", "strain_id", "project_id2", "plate_id"), sep = ":", remove = F)
  # temp = mutate(temp, plate_id = paste(project_id2, plate_id, sep = ":")) %>% select(-project_id2)
  # temp$concentration <- sub("\\.[[:alpha:]]", "", temp$concentration)
  # temp$pert_dose <- as.numeric(temp$pert_dose)
  # temp <- mutate(temp, pert_id = ifelse(grepl("BRD", pert_id), substr(pert_id, 1, 13), pert_id))
  colnames(temp)[1]= "CompoundCorrelation" 
  finaltable <- arrange(temp, desc(CompoundCorrelation))
  finaltable <- mutate(finaltable, order= 1:n(), CompCorPercentile = 100*order/dim(finaltable)[1]) %>% select(-order)
  
  #Add tanimoto coefficient if available
  if(!is.na(tanimoto_path)){
    row_meta_tani <- read.gctx.meta(tanimoto_path, dimension = "row")
    col_meta_tani <- read.gctx.meta(tanimoto_path, dimension = "col")
    pert_id1 = ifelse(grepl("BRD", comp1), substr(comp1, 1, 13), comp1)
    colid_tani = match(pert_id1, col_meta_tani$id) #may need to add fixed = T
    rowid_tani = match(finaltable$pert_id, row_meta_tani$id) #may need to add fixed = T
    # row_meta_tani <- mutate(row_meta_tani, rownum = 1:n())
    if(!is.na(colid_tani)){
      tani_vector <- parse.gctx(tanimoto_path, cid = colid_tani)
      tani_vector <- tani_vector@mat
      tani_vector_sub <- tani_vector[rowid_tani]
      finaltable$tanimoto = tani_vector_sub
    }else{
      print("Tanimoto coefficient not found")
    }
    
    # tani_vector_df <- data.frame(tani_vector)
    # tani_vector_df$pert_id <- rownames(tani_vector_df)
  }
  
  #Add compound table information (correlation, well count)
  finaltable2 <- left_join(finaltable, keepcomps, by = "id") #requires full id to be the same. May want to relax to pert_id
  finaltable2 = arrange(finaltable2, desc(CompoundCorrelation))
  
  return(finaltable2)
}

compCorFilter <- function(comp1, conc1, pert_itime1 = "0min", project_id1, strain_id1, compcortable, cormatrix_path, cormatrix_path2 = NA, metric = "median", compCorThreshold = 0.5, metricThreshold = 0.5, removequery = F, printPlots = F){
  #3/23/23 updated to accomodate new doses with 5 digits after the decimal, and change of project_id to have underscore instead of hyphen
  #Find Compound Community based on compound correlations
  #From a list of correlations (compcortable) between compounds and a query compound (output of compoundCorrelationdf),
  #filter out the ones of interest that are most similar to the rest of the group.
  #The column name in compcortable should be "CompoundCorrelation"
  #metric can be median, mean, min
  #compCorThreshold is the threshold for the compound correlations to the query compound
  #metricThreshold is the threshold for the metric, below which comp_conc are filtered out as not part of the group
  #cormatrix_path is the path to the gctx file for the similarity matrix containing the query compound vs kabx
  #cormatrix_path2 is optional to contain the similarity matrix for kabx against itself (if this is NA, then it assumes cormatrix_path includes all pairwise comparisons)
  #removequery = F. if it is TRUE, then removes the query compound (and all of its concentrations) prior to forming the community
  #assumes replicates are collapsed. Otherwise may need to delineate by project id.
  #assumes condition_id is now part of column metadata of cormatrix_path
  #id_col is the column name that the id should be matched to, usually condition_id or comb_id
  
  #Output: Returns data frame with median, average, min compound correlations between each compound and the rest of the group, 
  #CompoundCorrelation column corresponds to correlation with original query compound
  #1/10/20 edit. Second item in return list is the compound correlations within the entire community
  #1/29/20 edit. add comp1, conc1, which are the query compound and concentration, need to obtain from other correlation matrix
  #2/13/20 edit to accomodate cases where a known is correlated with a known from a different screen
  
  #*Note: As written, this function is not suitable for cases where want to keep correlation for different doses of the same query compound if it is not in the kabx2 screen. 
  #Should probably only be used for reference-based MOA
  
  #query
  # query_compconc = paste(project_id1, comp1, conc1, pert_itime1, strain_id1, sep = ":")
  query_compconc = tolower(paste(project_id1, comp1, conc1, strain_id1, pert_itime1, sep = ":")) #3/23/23
  
  cormatrix_rowmeta = read_gctx_meta(cormatrix_path, dim = "row")
  cormatrix_colmeta = read_gctx_meta(cormatrix_path, dim = "col")

  finaltable_filt = filter(compcortable, CompoundCorrelation > compCorThreshold)
  finaltable_filt = mutate(finaltable_filt, comp_conc = paste(project_id, pert_id, pert_idose, pert_itime, strain_id, sep =":"))
  
  if(removequery){
    querycomp_df = filter(finaltable_filt, CompoundCorrelation == 1)
    if(nrow(querycomp_df) > 0){
      finaltable_filt = filter(finaltable_filt, pert_id != querycomp_df$pert_id) #removes all concentrations of the query pert_id 
      finaltable_filt = rbind(querycomp_df, finaltable_filt) #adds back only the query itself
    } #assumes if the correlation is not 1, and it has the same compound and concentration, it was screened at a different time, not the same and should be compared
    #     else{
    #       finaltable_filt = filter(finaltable_filt, comp_conc != compconc1)
    #     }
  }
  
  # if((query_compconc %in% cormatrix_rowmeta$id) & (query_compconc %in% cormatrix_colmeta$id) & (length(cormatrix_rowmeta$id)==length(cormatrix_colmeta$id))){
  #Assume a symmetric correlation matrix with unknown and kabx compounds
  if(all(cormatrix_colmeta$id %in% cormatrix_rowmeta$id)){
    #all columns are in the rows
    cormatrix_filt = parse_gctx(cormatrix_path, rid = match(finaltable_filt$id, cormatrix_rowmeta$id), cid =  match(finaltable_filt$id, cormatrix_colmeta$id))
    cormatrix_filt = cormatrix_filt@mat
    
    cormatrix_filt_upper = cormatrix_filt
    cormatrix_filt_upper[lower.tri(cormatrix_filt_upper, diag = T)] = NA
  }else{
    #must reconstruct the matrix from two correlation matrices, one with the query, and the other with kabx
    #retrieve query from first correlation matrix column
    finaltable_filt_noquery = dplyr::filter(finaltable_filt, pert_id != comp1) 
    
    if(nrow(finaltable_filt_noquery) > 0){ #the query has a community outside of itself, same pert_id is allowed...
      #allowed to correlate to same treatment condition if it passed the removequery stage
      # query_col = parse_gctx(cormatrix_path, rid = match(finaltable_filt$id, cormatrix_rowmeta$id, nomatch = 0), cid = which(grepl(query_compconc, cormatrix_colmeta$id)))
      query_col = parse_gctx(cormatrix_path, rid = match(finaltable_filt$id, cormatrix_rowmeta$id, nomatch = 0), cid = which(grepl(query_compconc, cormatrix_colmeta$condition_id))) #3/23/23
      
      query_col = query_col@mat
      
      kabx_rowmeta = read_gctx_meta(cormatrix_path2, dim = "row")
      kabx_colmeta = read_gctx_meta(cormatrix_path2, dim = "col")
      kabx_mat = parse_gctx(cormatrix_path2, rid = match(finaltable_filt$id, kabx_rowmeta$id, nomatch = 0), cid = match(finaltable_filt$id, kabx_colmeta$id, nomatch = 0))
      kabx_mat = kabx_mat@mat
      
      # length(query_col) == dim(kabx_mat)[1]
      # all(rownames(query_col) == rownames(kabx_mat))
      cormatrix_filt = cbind(query_col, kabx_mat)
      temp = t(rbind(1, query_col))
      cormatrix_filt = rbind(temp, cormatrix_filt)
      colnames(cormatrix_filt)[1] = query_compconc
      rownames(cormatrix_filt)[1] = query_compconc
      
      cormatrix_filt_upper = cormatrix_filt
      cormatrix_filt_upper[lower.tri(cormatrix_filt_upper, diag = T)] = NA
      cormatrix_filt_upper[1,1]= 1
    }else{
      #Not suitable for cases where want to keep correlation for different doses of the same compound. 
      #Would need to add another correlation matrix for the screen against itself.
      
      print("no neighbors")
      cormatrix_filt_upper= 1
      # cormatrix_filt = matrix(c(finaltable_filt$CompoundCorrelation), nrow = 1, ncol = 1)
      # rownames(cormatrix_filt)[1] = query_compconc #should use id
      # colnames(cormatrix_filt)[1] = query_compconc
      # 
      # cormatrix_filt_upper = matrix(c(1), nrow = 1, ncol = 1)
      # rownames(cormatrix_filt_upper)[1] = query_compconc
      # colnames(cormatrix_filt_upper)[1] = query_compconc
    }
  }
  
  #Replace the 1 for the query compound
  
  colmedian = apply(cormatrix_filt_upper, 2, FUN = substitute(metric), na.rm = T)
  test = apply(cormatrix_filt_upper, 2, FUN = median, na.rm = T)
  idxend_all = which(colmedian <= metricThreshold)  #Choose up to the first compound that has median low compCor with rest of compounds in group
  if(length(idxend_all) == 0){
    #Even the last element fulfills the compound correlation requirement
    # & !unname((colmedian[length(colmedian)] <= compCorThreshold))
    idxend = length(colmedian)
  }else if(any(is.na(unname(idxend_all)))){
    #Only the first element (query compound) remains. Nothing else is similar enough
    idxend = 1
  }else{
    idxend = idxend_all[1] - 1 #Decrease one index since you want one with colmedian > compCorThreshold
  }
  
  cormatrix_filt2 = as.matrix(cormatrix_filt_upper[1:idxend,1:idxend])
  rownames(cormatrix_filt2) = rownames(cormatrix_filt_upper)[1:idxend]
  colnames(cormatrix_filt2) = colnames(cormatrix_filt_upper)[1:idxend]
  
  if(printPlots){
    hist(cormatrix_filt_upper, breaks = 50,xlim = c(0,1), xlab = "Compound Correlation", ylab = "# Comp_Conc")
    hist(cormatrix_filt2, xlim = c(0,1), breaks = 50, xlab = "Compound Correlation of Subset", ylab = "# Comp_Conc")
    
  }
  
  #return correlations of the community for later plotting purposes. 
  if(idxend == 1){
    #No compounds in the community
    #return empty data frame
    cormatrixfilt2_df = data.frame(comp_conc = character(), CommunityCorrelation = numeric())
  }else{
    cormatrixfilt2_df = as.data.frame(cormatrix_filt2[,2:idxend]) #exclude query
    cormatrixfilt2_df = gather(cormatrixfilt2_df, key = "id", value = "CommunityCorrelation")
    cormatrixfilt2_df = filter(cormatrixfilt2_df, !is.na(CommunityCorrelation))
  }
  
  # pdf('/broad/hptmp/jbagnall/cormatrix_K71125014.pdf', width = 70, height = 50)
  # pheatmap(cormatrix_filt)
  # dev.off()
  
  #This gives correlations to all compounds within the community
  cormatrixfilt_df = as.data.frame(cormatrix_filt[1:idxend, 1:idxend])
  rownames(cormatrixfilt_df) = rownames(cormatrix_filt)[1:idxend]
  colnames(cormatrixfilt_df) = colnames(cormatrix_filt)[1:idxend]
  cormatrixfilt_df = gather(cormatrixfilt_df, key = "id", value = "CompoundCorrelation")
  # cormatrixfilt_df$comp_conc =factor(cormatrixfilt_df$comp_conc, levels = unique(cormatrixfilt_df$comp_conc))
  #   ggplot(cormatrixfilt_df, aes(x = comp_conc, y = CompoundCorrelation))+
  #     geom_point()+
  #     theme_bw()+
  #     theme(axis.text.x = element_text(angle = 90))
  
  numcompconc = length(unique(cormatrixfilt_df$id))
  cormatrixfilt_df = cormatrixfilt_df %>%
    filter(CompoundCorrelation != 1) %>%
    group_by(id) %>%
    summarise(medianCompCor = median(CompoundCorrelation), avgCompCor = mean(CompoundCorrelation, na.rm = T), minCompCor = min(CompoundCorrelation, na.rm = T), maxCompCor = max(CompoundCorrelation, na.rm = T)) %>%
    ungroup()
  #     summarise(medianCompCor = median(CompoundCorrelation), avgCompCor = (sum(CompoundCorrelation)-1)/(numcompconc -1), minCompCor = min(CompoundCorrelation, na.rm = T)) %>%
  #     ungroup()
  
  
  #Add cormatrix_filt_upper?
  #This gives the correlation to the compounds up to that point (not to all compounds in community)
  cormatrixfilt_df_sub = as.data.frame(cormatrix_filt_upper[1:idxend, 1:idxend])
  rownames(cormatrixfilt_df_sub) = rownames(cormatrix_filt)[1:idxend]
  colnames(cormatrixfilt_df_sub) = colnames(cormatrix_filt)[1:idxend]
  cormatrixfilt_df_sub = gather(cormatrixfilt_df_sub, key = "id", value = "CompoundCorrelation")
  numcompconc = length(unique(cormatrixfilt_df_sub$id))
  cormatrixfilt_df_sub = cormatrixfilt_df_sub %>%
    filter(!is.na(CompoundCorrelation)) %>%
    group_by(id) %>%
    summarise(medianCompCor_sub = median(CompoundCorrelation, na.rm = T), avgCompCor_sub = mean(CompoundCorrelation, na.rm = T), minCompCor_sub = min(CompoundCorrelation, na.rm = T), maxCompCor_sub = max(CompoundCorrelation, na.rm = T)) %>%
    ungroup()
  
  # comb_df = left_join(cormatrixfilt_df, select(finaltable_filt, CompoundCorrelation, comp_conc), by  = "comp_conc")
  comb_df = full_join(finaltable_filt, cormatrixfilt_df, by = "id") #changed from right join to full_join 2/28/20
  comb_df = left_join(comb_df, cormatrixfilt_df_sub, by = "id")
  
  # comb_df = left_join(cormatrixfilt_df, finaltable_filt, by  = "comp_conc")
  # comb_df$comp_conc =factor(comb_df$comp_conc, levels = unique(comb_df$comp_conc))
  #   ggplot(comb_df, aes(x = comp_conc, y = avgCompCor, color = as.factor(medianCompCor > 0.5)))+
  #     geom_point()+
  #     theme_bw()+
  #     theme(axis.text.x = element_text(angle = 90))
  
  comb_df = arrange(comb_df, desc(CompoundCorrelation))
  
  return_list = list()
  return_list[[1]] = comb_df
  return_list[[2]] = cormatrixfilt2_df
  return(return_list)
  
}

targetCorrelation = function(cormatrix_path, query_ids = NA, reference_ids = NA, metadata, target_colname1 = "pert_mechanism", target_colname2 = "pert_target", savefilename = NA, removequery = F){
  #6/3/22 Adapted from screening data function. For reference-based MOA predictions
  #cormatrix_path is the correlation matrix between the query treatments (on columns) and the kabx2 knowns (on rows)
  #query_ids are the ids of the query treatments, if NA then takes all columns
  #reference_ids are the ids of the reference treatments (should match the ids in cormatrix), if NA then takes all rows
  #Go through all treatments of screen that satisfy certain conditions (repcorthreshold, nonWK)
  #Find maxmimally correlated treatment from each target class
  #Includes normalized and non-normalized values of the maximum correlation from each target class to the treatment
  #removeanyWK only feeds into calling active treatments, not into ranking of targets (there, well killers are always included)
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
  
  cormat = parse_gctx(cormatrix_path, cid = which(col_meta$id %in% query_ids), rid = which(row_meta$id %in% reference_ids)) 
  cormat_df = as.data.frame(cormat@mat)
  cormat_df$kabx_id = rownames(cormat_df)
  cormat_lf = tidyr::gather(cormat_df, key = "query_id", value = "CompoundCorrelation", -kabx_id)

  
  colnames(metadata)[colnames(metadata)==target_colname1] = "target"
  colnames(metadata)[colnames(metadata)==target_colname2] = "target2"
  
  #add compound metadata for query
  metadata_select = select(metadata, id, pert_id, pert_dose, pert_time, strain_id, project_id)
  colnames(metadata_select) = paste0("query_", colnames(metadata_select))
  cormat_lf_anno = left_join(cormat_lf, metadata_select, by = "query_id")
  
  #metadata for poscons
  colnames(metadata) = paste0("kabx_", colnames(metadata))
  cormat_lf_anno = left_join(cormat_lf_anno, metadata, by = "kabx_id")
  
  #Double check to remove compounds with no annotation
  #Make sure there are no redundancies, reduce to only target information
  cormat_lf_anno = filter(cormat_lf_anno, kabx_target != "-666" & kabx_target != "NA" & kabx_target != "unknown" & kabx_target != "whole cell only")
  cormat_lf_anno = filter(cormat_lf_anno, !is.na(kabx_target))
  cormat_lf_anno = distinct(cormat_lf_anno)
  
  
  if(removequery){#remove query pert_id, only needed if testing with knowns
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
  }

}

refBasedList_compound = function(target_list_path, rank_thresh = 1, normalized_cor = F, kabx_annopath= "", target_col = "pert_mechanism", roc_path = '/idi/cgtb/jbagnall/psa_RNAseq/poscon_correlation/R_datafiles/kabx_rocinfo_compound_include_singles_220607_pert_mechanism.rds', tani_path = '', tani_roc_path = '', savefilename, avg_doses = T, print_plot = F, save_all_target_cor = T){
  #6/7/22 List prediction for each compound(pert_id)
  ###2/13/20 Need to adjust to new column names from targetCorrelation
  #target_list is the output of targetCorrelation, where every query_pert_id has been correlated to a maximum from a target category from kabx2
  #it can either be a list or a data frame
  #normalized_cor == T, then uses ranks for the mean normalized correlation values (instead of the mean of the actual correlation values)
  #Should I make these optional inputs to define a high confidence region? norm_compcor_thresh = 0.95, compcor_thresh = 0.8
  #if anno_name_change == T, then changes annotation columns to new names (1/31/22)
  #target_col should match the name of target column1 from targetCorrelation (used to define target)
  #annotatoin file requires pert_type = poscon, negcon or test. defines "unknowns" as "test" compounds
  #if avg_doses = T, then wil use the mean of the correlations across doses (this is default). If F, then will use the max across all doses.
  #if save_all_target_cor = T, then will save RDS files of max correlations to all targets, useful for plotting purposes later
  
  #outputs top ranked target predictions
  #subsets known vs unknowns
  #subsets high confidence vs everything else
  #saves 2 csv files, one with all knowns, one with unknowns
  #List all correlation values and our confidence score (what fraction with that correlation value that were called correctly in kabx2)
  
  if(is.character(kabx_annopath)){
    kabx_annotation = read.csv(kabx_annopath, stringsAsFactors = F)
  }else{
    #assumes the data frame or table has been passed in
    kabx_annotation = kabx_annopath
  }
  
  
  #   kabx_categories = kabx_annotation 
  #   # kabx_categories = replace(kabx_categories, kabx_categories == "-666", "unknown")
  #   kabx_categories = kabx_categories %>%  
  #     select(target_gene, target_details, target_protein_enzyme_or_complex, target_process) %>%
  #     distinct() %>%
  #     group_by(target_details) %>%
  #     summarise(target_gene = paste(unique(target_gene), collapse = "|"), target_protein_enzyme_or_complex = paste(unique(target_protein_enzyme_or_complex), collapse = "|"), target_process = paste(unique(target_process), collapse = "|")) %>%
  #     ungroup()
  
  
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
    
    #Add FDR from kabx2 set
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
    
    #Add tanimoto information if given
    if(tani_path != ''){
      row_meta = read.gctx.meta(tani_path, dimension = "row")
      col_meta = read.gctx.meta(tani_path, dimension = "col")
      
      
      unknown_df_output_anno = unknown_df_output_anno %>%
        group_by(query_pert_id, neighbor_pert_id) %>%
        mutate(tanimoto = ifelse((query_pert_id %in% row_meta$id) & (neighbor_pert_id %in% col_meta$id), as.numeric(parse.gctx(tani_path, rid = which(row_meta$id == query_pert_id), cid = which(col_meta$id == neighbor_pert_id))@mat), NaN)) %>%
        ungroup()
    }else{
      unknown_df_output_anno$tanimoto = NaN
    }
    
    #collapse neighbors again
    unknown_df_output_anno = select(unknown_df_output_anno, -target_rank_norm, -target_rank_notnorm, -threshold_idx, -fdr)
    
    unknown_df_output_anno = unknown_df_output_anno %>%
      group_by(query_pert_id, kabx_target, mean_compcor, mean_compcor_norm, precision) %>%
      summarise(neighbor_pert_id = paste(unique(neighbor_pert_id), collapse = "|"), pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|"), max_tani = max(tanimoto, na.rm = T), mean_tani = mean(tanimoto, na.rm = T), median_tani = median(tanimoto, na.rm = T)) %>%
      ungroup()
    colnames(unknown_df_output_anno) = c("query_pert_id", "predicted_target", "target_correlation", "normalized_target_correlation", "fraction_correct_in_kabx2", "neighbor_pert_id", "neighbor_pert_iname", "predicted_pert_mechanism", "predicted_pert_target", "predicted_pharmaceutical class", "max_tani", "mean_tani", "median_tani")
    
    
    
    if(tani_path != ''){
      #Add max_tanimoto fdr/precision information from kabx2 set
      tani_rocinfo = readRDS(tani_roc_path)
      
      # test = unknown_df_output
      unknown_df_output_anno$threshold_idx = sapply(unknown_df_output_anno$max_tani, function(x){max(which(x > tani_rocinfo$threshold))})
      unknown_df_output_anno$fraction_correct_kabx2_tani = tani_rocinfo$precision[unknown_df_output_anno$threshold_idx]
      unknown_df_output_anno = select(unknown_df_output_anno, -threshold_idx)
    }else{
      unknown_df_output_anno$fraction_correct_kabx2_tani = NA
    }
    unknown_df_output_anno = arrange(unknown_df_output_anno, desc(target_correlation))
    
    write.csv(unknown_df_output_anno, paste(savefilename, "_unknown.csv", sep = ''), row.names = F)
    
    
    #   unknown_df_output_anno = left_join(unknown_df_output, kabx_categories, by = c("target" = "target_details"))
    #   unknown_df_output_anno = select(unknown_df_output_anno, -target_rank_norm, -target_rank_notnorm, -threshold_idx, -fdr)
    #   colnames(unknown_df_output_anno) = c("query_pert_id", "predicted_target", "target_correlation", "normalized_target_correlation", "fraction_correct_in_kabx2", "predicted_target_gene", "predicted_target_protein_enzyme_or_complex", "predicted_target_process")
    #   unknown_df_output_anno = arrange(unknown_df_output_anno, desc(target_correlation))
    #   write.csv(unknown_df_output_anno, paste(savefilename, "_unknown.csv", sep = ''), row.names = F)
    
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
      
      if(tani_path != ''){
        p2 = ggplot(unknown_df_output_anno, aes(x = target_correlation, y = max_tani))+
          geom_point()+
          labs(x = "Average correlation", y = "Max Tanimoto", title = "Unknowns")+
          xlim(0,1)+
          ylim(0,1)+
          theme_bw(base_size = 16)
        print(p2)
      }
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
    if(tani_path != ''){
      row_meta = read.gctx.meta(tani_path, dimension = "row")
      col_meta = read.gctx.meta(tani_path, dimension = "col")
      
      known_df_output_anno = known_df_output_anno %>%
        group_by(query_pert_id, neighbor_pert_id) %>%
        mutate(tanimoto = ifelse((query_pert_id %in% row_meta$id) & (neighbor_pert_id %in% col_meta$id), as.numeric(parse.gctx(tani_path, rid = which(row_meta$id == query_pert_id), cid = which(col_meta$id == neighbor_pert_id))@mat), NaN)) %>%
        ungroup()
    }else{
      known_df_output_anno$tanimoto = NaN
    }
    
    #collapse neighbors again
    known_df_output_anno = select(known_df_output_anno, -target_rank_norm, -target_rank_notnorm, -threshold_idx, -fdr)
    
    known_df_output_anno = known_df_output_anno %>%
      group_by(query_pert_id, kabx_target, mean_compcor, mean_compcor_norm, precision, query_pert_iname, query_pert_mechanism, query_pert_target, query_pharmaceutical_class) %>%
      summarise(neighbor_pert_id = paste(unique(neighbor_pert_id), collapse = "|"), pert_iname = paste(unique(pert_iname), collapse = "|"), pert_mechanism = paste(unique(pert_mechanism), collapse = "|"), pert_target = paste(unique(pert_target), collapse = "|"), pharmaceutical_class = paste(unique(pharmaceutical_class), collapse = "|"), max_tani = max(tanimoto, na.rm = T), mean_tani = mean(tanimoto, na.rm = T), median_tani = median(tanimoto, na.rm = T)) %>%
      ungroup()
    
    colnames(known_df_output_anno) = c("query_pert_id", "predicted_target", "target_correlation", "normalized_target_correlation", "fraction_correct_in_kabx2", "query_pert_iname", "query_pert_mechanism", "query_pert_target", "query_pharmaceutical_class", "neighbor_pert_id", "neighbor_pert_iname", "predicted_pert_mechanism", "predicted_pert_target", "predicted_pharmaceutical_class", "max_tani", "mean_tani", "median_tani")
    
    known_df_output_anno = arrange(known_df_output_anno, desc(target_correlation))
    
    colnames(known_df_output_anno)[colnames(known_df_output_anno) == paste0("query_", target_col)] = "query_target"
    known_df_output_anno = known_df_output_anno %>%
      group_by(query_pert_id, predicted_target, neighbor_pert_id) %>%
      mutate(same_target = grepl(predicted_target, query_target, fixed = T)) %>%
      ungroup
    colnames(known_df_output_anno[colnames(known_df_output_anno) == "query_target"]) = paste0("query_", target_col)
    #   known_df_output_anno = left_join(known_df_output_anno, kabx_categories, by = c("target"="target_details"))
    #   known_df_output_anno = mutate(known_df_output_anno, same_target = target == query_target_details, same_process = query_target_process == target_process)
    #   colnames(known_df_output_anno) = c("query_pert_id", "predicted_target","target_correlation", "normalized_target_correlation", "target_rank_norm", "target_rank_notnorm", "threshold_idx", "fraction_incorrect_in_kabx2", "fraction_correct_in_kabx2", "query_pert_iname", "query_target_gene", "query_target_details", "query_target_protein_enzyme_or_complex", "query_target_process", "predicted_target_gene", "predicted_target_protein_enzyme_or_complex", "predicted_target_process", "same_target", "same_process")
    #   known_df_output_anno = select(known_df_output_anno, -fraction_incorrect_in_kabx2, -target_rank_norm, -target_rank_notnorm, -threshold_idx)
    #   known_df_output_anno = arrange(known_df_output_anno, desc(target_correlation))
    
    if(tani_path != ''){
      #Add max_tanimoto fdr/precision information from kabx2 set
      tani_rocinfo = readRDS(tani_roc_path)
      
      # test = unknown_df_output
      known_df_output_anno$threshold_idx = sapply(known_df_output_anno$max_tani, function(x){max(which(x > tani_rocinfo$threshold))})
      known_df_output_anno$fraction_correct_kabx2_tani = tani_rocinfo$precision[known_df_output_anno$threshold_idx]
      known_df_output_anno = select(known_df_output_anno, -threshold_idx)
    }else{
      known_df_output_anno$fraction_correct_kabx2_tani = NA
    }
    
    if(print_plot & tani_path != ''){
      p3 = ggplot(known_df_output_anno, aes(x = target_correlation, y = max_tani, color = same_target))+
        geom_point()+
        labs(x = "Average correlation", y = "Max Tanimoto", title = "Knowns")+
        xlim(0,1)+
        ylim(0,1)+
        theme_bw(base_size = 16)
      print(p3)
    }

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
  mech_color_df = data.frame(pert_mechanism = c("Negative control", "DNA synthesis", "Membrane Integrity", "Peptidoglycan biogenesis", "Protein synthesis"), mech_color = c("black", "red", "blue","green4","magenta3"))
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
    labs(x = "Target", y = "Pearson Correlation", title = paste0(pert_iname0, " (Closest drug: ", closest_compound, ")"), color = "Mechanism")+
    facet_grid(~predicted_pert_mechanism_factor, space="free_x", scales = "free_x")+
    ylim(-0.2,1)+
    scale_color_manual(values = c("Negative control" = "black", "DNA synthesis" = "red", "Membrane Integrity" = "blue", "Peptidoglycan biogenesis" = "green4", "Protein synthesis" = "magenta3"))+
    theme_bw(base_size = text_size)+
    theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text = element_blank())
  
  print(p0)
  
  
  #Color the points by Mechanism
  p1 = ggplot(max_cor_df_filt, aes(x = reorder(predicted_target, -target_correlation), y = target_correlation, color = predicted_pert_mechanism))+
    geom_point(size = 2)+
    #geom_text(label = paste0("Closest compound: \n", closest_compound), x = 18, y = 0.95, size = 4.5, font_face = "plain")+
    geom_hline(yintercept = 0, linetype = "dotted", color = "darkgray")+
    labs(x = "Target", y = "Pearson Correlation", title = pert_iname0, color = "Mechanism")+
    ylim(-0.2,1)+
    scale_color_manual(values = c("Negative control" = "black", "DNA synthesis" = "red", "Membrane Integrity" = "blue", "Peptidoglycan biogenesis" = "green4", "Protein synthesis" = "magenta3"))+
    theme_bw(base_size = text_size)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = color_vec))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(p1)
  
  plot_list = list(p0, p1)
  
  if(!is.na(save_file_path)){
    pdf(save_file_path, width = 7, height = 5)
    print(p0)
    dev.off()
  }
  
  return(plot_list)
}
