#Functions for calculating moderated z-scores
#jbagnall@broadinstitute.org
#Some functions adapted from cmapM package written in Matlab. https://github.com/cmap/cmapM/blob/master/sig_tools/math/modzs.m

library(cmapR)
library(DESeq2)

trim_matrix = function(data, remove_samples = NA, remove_genes = NA, savefilepath = NA){
  #Data is a matrix
  #remove_samples is a character vector of samples to be removed due to low quality--can be a filepath
  #remove_genes is a character vector of gene_ids to be removed due to low counts (half or more are zero)--can be a filepath
  #Returns matrix trimmed on rows and columns
  
  if(grepl(".txt", remove_samples)){
    if(file.size(remove_samples) != 0){
      remove_samples = read.table(remove_samples)$V1
    }else{
      remove_samples = NA
    }
    
  }
  if(grepl(".txt", remove_genes)){
    if(file.size(remove_genes) != 0){
      remove_genes = read.table(remove_genes)$V1
    }else{
      remove_genes = NA
    }
    
  }
  
  if(!(all(is.na(remove_genes)) | all(remove_genes == ""))){
    if(!all(remove_genes %in% rownames(data))){
      stop("Error: Not all genes in remove_genes are in the data set")
    }
  }
  if(!(all(is.na(remove_samples)) | all(remove_samples == ""))){
    if(!all(remove_samples %in% colnames(data))){
      stop("Error: Not all samples in remove_samples are in the data set")
    }
  }
  data_filt = data[!(rownames(data) %in% remove_genes), !(colnames(data) %in% remove_samples)]
  
  if(is.na(savefilepath)){
    return(data_filt)
  }else{
    saveRDS(data_filt, savefilepath)
  }
}

run_pca = function(data, col_meta = NA, color_col = "pert_mechanism", scale0 = T){
  #Input data should be genes on rows and samples on columns
  #Output returns dataframe with PC values, and prints plot of first 4 PCs
  #For general QC and visualization
  
  pca_scaled <- prcomp(t(data), scale = scale0) 
  pca_scaled.var <- pca_scaled$sdev^2
  pca_scaled.var.per <- round(pca_scaled.var/sum(pca_scaled.var)*100, 1)
  barplot(pca_scaled.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
  
  screendata_pcascaled <- data.frame(sample_id = rownames(pca_scaled$x), pc1 = pca_scaled$x[,1], pc2 = pca_scaled$x[,2], pc3 = pca_scaled$x[,3], pc4 = pca_scaled$x[,4], stringsAsFactors = FALSE)
  
  if(is.na(col_meta)){
    p0 = ggplot(data = screendata_pcascaled, aes(x = pc1, y = pc2)) +
      geom_point(alpha = 0.2)+
      xlab(paste("PC1 - ", pca_scaled.var.per[1], "%", sep="")) +
      ylab(paste("PC2 - ", pca_scaled.var.per[2], "%", sep="")) +
      theme_bw() 
    print(p0)
    
    p1 = ggplot(data = screendata_pcascaled, aes(x = pc3, y = pc4)) +
      geom_point(alpha = 0.2)+
      xlab(paste("PC3 - ", pca_scaled.var.per[3], "%", sep="")) +
      ylab(paste("PC4 - ", pca_scaled.var.per[4], "%", sep="")) +
      theme_bw() 
    print(p1)
  }else{
    #color by color_col
    screendata_pcascaled = left_join(screendata_pcascaled, col_meta, by = c("sample_id" = "id"))
    p0 = ggplot(data = screendata_pcascaled, aes(x = pc1, y = pc2)) +
      geom_point(alpha = 0.2, aes_string(color = color_col))+
      xlab(paste("PC1 - ", pca_scaled.var.per[1], "%", sep="")) +
      ylab(paste("PC2 - ", pca_scaled.var.per[2], "%", sep="")) +
      theme_bw() 
    print(p0)
    
    p1 = ggplot(data = screendata_pcascaled, aes(x = pc3, y = pc4)) +
      geom_point(alpha = 0.2, aes_string(color = color_col))+
      xlab(paste("PC3 - ", pca_scaled.var.per[3], "%", sep="")) +
      ylab(paste("PC4 - ", pca_scaled.var.per[4], "%", sep="")) +
      theme_bw() 
    print(p1)
  }

  
  
  return(screendata_pcascaled)
}

qnorm_median = function(values){
  #Adapted quantilenorm function from Matlab
  #Quantile normalization for each column of the input data
  #Input values is typically a matrix of moderated z-scores
  #Need to pick which factors to group the population before calling this function
  #For each column, the rows are sorted and ranked, and then the ranks are aligned to the median for each row.
  #Returns quantile normalized version of values, normalized per column
  
  #initialize
  tiedFlag = TRUE
  valSize = dim(values)
  ndx = matrix(data = 1, nrow = valSize[1], ncol = valSize[2])
  N = valSize[1] #num rows
  nanvals = is.na(values)
  numNans = colSums(nanvals)
  rankedVals = matrix(data = NaN, nrow = valSize[1], ncol = valSize[2])
  normalizedVals = values
  
  # create space for output
  if(tiedFlag){
    rr = list()
  }


  # for each column we want to ordered values and the ranks with ties
  for(col in 1:valSize[2]){
    sortedVals = sort(values[,col])
    ndx[,col] = order(values[,col])
    if(tiedFlag){
      rr[[col]] = sort(rank(values[which(!is.na(values[,col])),col], ties.method = "average"))
    }
    M = N-numNans[col];
    #interpolate over the non-NaN data to get ranked data for all points
    rankedVals[,col] = pracma::interp1(x = seq(from = 1, to = N, by = (N-1)/(M-1)), y = sortedVals[1:M], xi = seq(from = 1, to = N, by = 1), method = "linear");
  }
  
  #take the median of the ranked values for each row
  median_vals = rowMedians(rankedVals, na.rm = T)

  #Extract the values from the normalized distribution
  for(col in 1:valSize[2]){
    M = N-numNans[col]
    if(tiedFlag){
      normalizedVals[ndx[1:M,col],col] = pracma::interp1(seq(1, N, by = 1), median_vals, 1+((N-1)*(rr[[col]]-1)/(M-1)), method = "linear")
      # normalizedVals[ndx(1:M,col),col] = pracma::interp1(1:N,median_vals,1+((N-1)*(rr{col}-1)/(M-1)));
    }else{
      normalizedVals[ndx(1:M,col),col] = pracma::interp1(seq(1, N, by = 1),median_vals,1:(N-1)/(M-1):N);
    }
  }
  
  return(normalizedVals)
}

robust_zscore = function(data, dim0 = 1, min_mad = pracma::eps(x=1), estimate_prct = 0.01){
  #Adapted from cmapM::robust_zscore.m
  #Input: data is a matrix
  #Z-score along rows for dim0 = 1 or columns for dim0 = 2
  #min_mad = pracma::eps(x = 1), minimum mad is smallest floating point to avoid dividing by zero
  #Z-score is calculated for each column of data, where zscore = (data -median(data))/mad(data)
  #Note that in R, stats::mad already includes the 1.4826 scaling factor, so it is not explicitly included as it is in the matlab code
  #estimate_prct is the percentile of the data set for the lowest mad value, 0.01 is default to 1%
  #Outputs centered, scaled version of data (z-score), same size as data
  #If data is empty , return empty
  
  median_vec = apply(X = data, MARGIN = dim0, FUN = function(x){stats::median(x, na.rm = T)})
  sigma_vec = apply(X = data, MARGIN = dim0, FUN = function(x){stats::mad(x, na.rm = T)})
  
# Adjust for low variance by using the "estimate" method. 'estimate': The default method. Estimate the minimum
# variance from the data. EST_MAD = PRCTILE(SIGMA, est_prct), where SIGMA is a vector of MAD computed
# along dim for the entire dataset. The minimum variance is set to: MAX(EST_MAD, MIN_MAD).
  min_mad = max(quantile(as.vector(sigma_vec), estimate_prct, type = 5), min_mad, na.rm = T); #matlab percentile is similar to quantile type 5
  
  sigma_vec0 = base::pmax(sigma_vec, min_mad)
  
  data_zs = sweep(data, MARGIN = dim0, STATS = median_vec, FUN='-')
  data_zs = sweep(data_zs, MARGIN = dim0, STATS = sigma_vec0, FUN='/')
  
  return_list = list(data_zs, median_vec, sigma_vec0)
  return(return_list)
}

modzs_collapse_avg = function(data, ridx){
  #Adapted from cmapM::modzs.m (moderated z-score), where they take a weighted average of replicates
  #Input data & row indices
  #Performs spearman correlation across the columns (columns should be replicates)
  #Output collapsed values
  
  #defaults from matlab code
  low_thresh_wt = 0.01
  low_thresh_cc = 0
  
  if(is.vector(data)){
    print(paste("No replicate, ridx:", ridx))
    czs = data;
    norm_wt = 1;
    cc = 1;
  }else{
    dim_data = dim(data)
    nr = dim_data[1]
    nc = dim_data[2]
    
    if(nc > 1){
      this_zs = data[ridx,]
      num_nan = sum(is.na(this_zs)) #number of NA elements
      if(num_nan == 0){
        cc = cor(this_zs, method = "spearman")
      }else{
        warning(paste0("Missing values in data, using pairwise-mode in correlation: ", num_nan))
        cc = cor(this_zs, method = "spearman", use = "pairwise.complete.obs")
      }
      
      #Replace diagonal with zeros
      base::diag(cc) = 0
      wt = 0.5*rowSums(cc); #not sure why the half
      
      #clip low weights--this is not default for the matlab function
      #wt = base::pmax(wt, low_thresh_wt)
      
      #use wt_avg
      #Make sure the weights sum to 1, so normalize them
      sum_wt = sum(abs(wt)) 
      # norm_wt = pracma::mrdivide(wt,sum_wt) #matlab: wt/sum_wt, could be same as scalar divide
      norm_wt = wt / sum_wt #assuming sum_wt is scalar
      
      #replace NA with 0
      data = replace(data, which(is.na(data)), values = 0)
      czs = data %*% norm_wt #matrix multiply
    }else{
      czs = data;
      norm_wt = 1;
      cc = 1;
    }
  }
  #Returns collapsed replicates. For now doesn't return weights and correlation matrix
  return(czs) 
  
}