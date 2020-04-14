###################
# Functions       #
###################

# This file contains all functions that are needed behind to process RT-qPCR data.

# Include the class definitions
#' @include classes.R

#######

#' Read in qPCR data.
#'
#' \code{read_qpcr} reads in and formats the RT-qPCR expression data.
#'
#' @name read_qpcr
#' @title read_qpcr
#' @param qpcr An empty instance of the qpcr class.
#' @param files A vector of paths pointing to qPCR raw data files.
#' @param keep_undetected A logical indicating whether to keep undetected CT values and set them to 40, or to simply delete those entries. Defaults to TRUE.
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom readxl read_excel
#' @return An updated instance of the qpcr class.
read_qpcr<-function(qpcr,files,keep_undetected=TRUE){
  # Read in the data
  data<-lapply(files, readxl::read_excel, sheet="Results", skip=35)
  # Bind rows
  data<-lapply(data, FUN=function(x){
    x[]<-lapply(x,as.character)
  })
  data<-data.frame(dplyr::bind_rows(data),stringsAsFactors = F)
  # Format
  data<-data[,c("Sample.Name","Target.Name","CT")]
  data<-data[complete.cases(data),]
  # Keep or discard undetected values
  if(keep_undetected){
    data$CT[data$CT=="Undetermined"]<-40
  } else {
    data<-data[data$CT!="Undetermined",]
  }
  data$CT<-as.numeric(data$CT)
  names(data)<-c("Sample","Gene","CT")
  qpcr@ct<-data
  return(qpcr)
}

#' Add sample annotation
#'
#' \code{add_annot} allows one to add a sample annotation file to the qpcr object.
#'
#' @name add_annot
#' @title add_annot
#' @param qpcr An instance of the qpcr class.
#' @param file A path to a sample annotation file. This file should at least have one column named 'Sample'.
#' @return An qpcr object with added sample annotation.
#' @importFrom readxl read_excel
#' @export
add_annot<-function(qpcr, file){
  # Read in the annotation file
  annot<-readxl::read_excel(file, sheet = 1)
  if(!any(names(annot)=="Sample")){
    stop("Please provide an annotation file with at least one column (named 'Sample').")
  }
  # Add annotation to the object
  qpcr@annotation<-data.frame(annot,stringsAsFactors = F)
  return(qpcr)
}

# Calculation of geometric mean

geomean<-function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Mark undetected values

add_nd<-function(x){
  if(sum(x==40)>sum(x!=40)){
    res<-"ND"
  } else {
    res<-""
  }
  return(res)
}

# Resolve technical replicates and reference genes

resolve<-function(x, aggregation){
  if(sum(x==40)>sum(x!=40)){
    res<-40
  } else {
    res<-switch (aggregation,
      "geomean" = geomean(x[x!=40]),
      "mean" = mean(x[x!=40]),
      "median" = median(x[x!=40])
    )
  }
  return(res)
}


#' Analyse qpcr data
#'
#' \code{analyse_qpcr} allows one to perform the core analysis of qpcr data.
#'
#' @name analyse_qpcr
#' @title analyse_qpcr
#' @param qpcr An instance of the qpcr class. Expression data and sample annotations must have been loaded already.
#' @param reference_genes A vector listing the genes one wishes to use as housekeeping genes or reference genes.
#' @param levels A list with the set levels for each factor in the sample annotation data. Defaults to NULL.
#' @param aggregation Defines how CT values of technical replicates and reference genes should be aggregated; either by their arithmetic ('mean') or geometric mean ('geomean') or their median ('median').
#' @return An updated instance of the qpcr class with analysed expression data.
#' @import dplyr
#' @export
analyse_qpcr<-function(qpcr, reference_genes, levels=NULL, aggregation = c("geomean","mean","median")){
  # Check arguments
  if(!all(reference_genes%in%qpcr@ct$Gene)){
    stop("Not all reference genes are present in the provided expression data.")
  }
  # Resolve technical replicates
  qpcr@ct<-qpcr@ct%>%
    group_by(Sample,Gene)%>%
    mutate(tech_CT = CT,
           CT = resolve(CT, aggregation),
           Detected = add_nd(CT),
           sd_CT = sd(tech_CT, na.rm = T))
  # Separate out the housekeeping genes
  qpcr@references<-qpcr@ct[qpcr@ct$Gene%in%reference_genes,]
  qpcr@ct<-qpcr@ct[!qpcr@ct$Gene%in%reference_genes,]
  # TODO QC and filter reference genes
  # Aggregate reference genes
  qpcr@references<-qpcr@references%>%
    group_by(Sample)%>%
    mutate(Reference = resolve(CT, aggregation))
  # Merge CTs of genes of interest with the housekeeping genes
  qpcr@dct<-merge(qpcr@ct[,c("Sample","Gene","CT")], qpcr@references[,c("Sample","Reference")], by="Sample", all = F)
  qpcr@dct<-qpcr@dct[!duplicated(qpcr@dct),]
  # Calculate delta CT values
  qpcr@dct$DCT<-qpcr@dct$CT-qpcr@dct$Reference
  return(qpcr)
}








