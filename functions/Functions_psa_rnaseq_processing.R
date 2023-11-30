#Functions for RNAseq data count pre-processing
#jbagnall@broadinstitute.org
#11/92/21

library(dplyr)
library(ggplot2)

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
