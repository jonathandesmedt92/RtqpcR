#################
# Classes       #
#################

# This file contains all classes that are needed behind to process RT-qPCR data.

# qpcr class

#' An S4 class to contain and process RT-qPCR data
#'
#' @export qpcr
#'
#' @slot ct A dataframe containing the raw CT values.
#' @slot references A dataframe containing the reference expression data.
#' @slot dct A dataframe containing the delta CT values
#' @slot annotation A dataframe containing the sample annotation data.
#' @slot melting A dataframe containing the melting curve data of each sample.

qpcr<-setClass("qpcr",
               slots = c(
                 ct = "data.frame",
                 references = "data.frame",
                 dct = "data.frame",
                 annotation = "data.frame",
                 melting = "data.frame"
               ))
