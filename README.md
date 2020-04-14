RtqpcR: Analysis of RT-qPCR expression data.
--------------------------------------------

This package facilitates the analysis of RT-qPCR expression data
(reverse transcriptase quantitative polymerase chain reaction). Here
below we provide a tutorial for a standard qPCR analysis.

Dependencies
------------

Before starting any qPCR analysis, please install the following
packages:

    BiocManager::install("readxl")
    BiocManager::install("dplyr")
    BiocManager::install("ggplot2")

Installation
------------

Install the RTqpcR package by running:

    devtools::install_github("jonathandesmedt92/RtqpcR")
    library(RtqpcR)

Tutorial
--------

Reading the raw expression data.
================================

To date, RtqpcR is only compatible with qPCR data from a Viaa7 platform
(384-well plates). The first step is to read in the expression data. At
this point, the user has to decide how to deal with undetected mRNA
expression values. We generally recommend setting these expression
values to 40 (i.e. maximal CT value). These values will be plotted in
graphs, however they will not be included in statistical analyses.

    # Initialise an RtqpcR object
    obj<-qpcr()

    # Read in the expression data
    obj<-read_qpcr(qpcr = obj, files = "data/test.xlsx")

Including Plots
---------------

You can also embed plots, for example:

![](README_files/figure-markdown_strict/pressure-1.png)

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
