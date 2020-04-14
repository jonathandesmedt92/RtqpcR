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
#' @param qpcr An empty instance of the qpcr class.
#' @param files A vector of paths pointing to qPCR raw data files.
#' @param keep_undetected A logical indicating whether to keep undetected CT values and set them to 40, or to simply delete those entries. Defaults to TRUE.
#' @export
#' @import dplyr
#' @import readxl
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
  data$CT[data$CT=="Undetermined"]<-40
  data$CT<-as.numeric(data$CT)
  names(data)<-c("Sample","Gene","CT")
  qpcr@ct<-data
  return(data)
}

