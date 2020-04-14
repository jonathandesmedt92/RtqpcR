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
  qpcr@dct<-merge(qpcr@ct[,c("Sample","Gene","CT","Detected")], qpcr@references[,c("Sample","Reference")], by="Sample", all = F)
  qpcr@dct<-qpcr@dct[!duplicated(qpcr@dct),]
  # Calculate delta CT values
  qpcr@dct$DCT<-qpcr@dct$CT-qpcr@dct$Reference
  return(qpcr)
}

# t test for fold changes

fc_ttest<-function(x,y){
  x<--log2(x)
  y<--log2(y)
  return(t.test(x,y))
}



#' Create bar plot
#'
#' \code{bar_plot} makes a bar plot of the expression data
#'
#' @name bar_plot
#' @title bar_plot
#' @param qpcr An instance of the qpcr class
#' @param xvar Variable to be displayed on the x-axis. This needs to be a column n the annotation data.
#' @param baseline_samples Samples that should be used as baseline for foldchange calculation.
#' @param genes Vector of the genes that will be included in the plot.
#' @param comparisons List of comparisons that will be evaluated with the indicated statistical test. Each element of this list should contain a vector of length two, containing two values of the variable displayed on the x-axis.
#' @param xlabels Named vector containing aliases for the values of the x-axis variable. These aliases will be used for plotting.
#' @param linear Logical indicating whether the expression values should be in a linear range. If FALSE, log2 values will be plotted.
#' @param y_breaks Numerical vector with custom y-axis breaks.
#' @param map_signif_level Logical indicating whether to use asterisks to denote significance. If FALSE, p-values will be displayed.
#' @param legend Logical indicating whether a legend needs to be displayed.
#' @param step_increase Offset for each significance line.
#' @param tip_length Tip length of each significance line.
#' @param bar_fill Name of the variable, based on which the bars will be colored.
#' @param colors Colors to fill the bars with.
#' @param ND_y_nudge Nudge value for the display of labels for undetected values.
#' @import ggplot2
#' @import dplyr
#' @import ggsignif
#' @import ggrepel
#' @export
bar_plot<-function(qpcr, xvar="Condition", baseline_samples, genes, ND_y_nudge=0.5, comparisons=NULL, xlabels=NULL, linear=T, y_breaks=NULL, map_signif_level=F, legend=T, step_increase=0, tip_length=0, bar_fill, colors=c("#bfbfbfbf","#666666")){
  # Merge dct with annotation
  dat<-merge(qpcr@dct,qpcr@annotation, by="Sample", all=F)
  dat[,names(dat)!="DCT"][]<-lapply(dat[,names(dat)!="DCT"][], as.character)
  # Calculate DDCT, fold change, mean fold change
  dat<-dat%>%
    group_by(Gene, get(xvar))%>%
    mutate(meanDCT = mean(DCT))%>%
    group_by(Gene)%>%
    mutate(DDCT = DCT-unique(meanDCT[Sample%in%baseline_samples]),
           Foldchange = 2^(-DDCT),
           Meanfoldchange = 2^(unique(meanDCT[Sample%in%baseline_samples])-meanDCT))%>%
    group_by(Gene, get(xvar))%>%
    mutate(Nr_detected = sum(CT!=40))%>%
    group_by(Gene)%>%
    mutate(MFC_adj = if(any(unique(Nr_detected)>1)){Meanfoldchange}else{1})
  # Filter for the genes required
  plotdat<-dat[dat$Gene%in%genes,]
  if(bar_fill==xvar){
    nddat<-plotdat[,c("Gene",xvar,"Foldchange","Detected")]
  } else {
    nddat<-plotdat[,c("Gene",xvar,"Foldchange","Detected",bar_fill)]
  }
  if(any(nddat$Detected=="ND")){
    nddat<-nddat%>%
      group_by(Gene, get(xvar))%>%
      filter(any(Detected=="ND"))%>%
      mutate(Detected="ND")%>%
      ungroup(Gene, get(xvar))%>%
      group_by(Gene, get(xvar), Detected)%>%
      summarise(Ypos = max(Foldchange)+ND_y_nudge)
  } else {
    nddat$Ypos<-max(nddat$Foldchange)+ND_y_nudge
  }
  names(nddat)[2]<-xvar
  p<-ggplot2::ggplot(plotdat, aes(x = get(xvar)))+
    geom_bar(stat="identity", position = position_dodge(), aes(y = MFC_adj, fill=get(bar_fill)))+
    geom_point(data=plotdat,aes(y=Foldchange, group=get(bar_fill)))+
    scale_fill_manual(values=colors)+
    theme_classic()+
    xlab("")+
    ylab("Relative expression")+
    theme(axis.text.x = element_text(angle=70, hjust=1, face="bold",size=10),
          axis.text.y = element_text(face="bold",size=10),
          legend.title = element_blank())+
    geom_hline(yintercept = 1)+
    ggrepel::geom_text_repel(data=nddat,aes(x=get(xvar),y=Ypos,label=Detected), size=3, fontface="bold")+
    facet_grid(.~Gene, scales = "free")
  if(!linear){
    p<-p+
      scale_y_continuous(trans = "log2", breaks = y_breaks)+
      suppressWarnings(ggsignif::geom_signif(data=plotdat,
                  aes(x=get(xvar),y=Foldchange),comparisons = comparisons,
                  map_signif_level = map_signif_level,
                  test="t.test",
                  tip_length = tip_length,
                  step_increase = step_increase))
  } else {
    p<-p+
      suppressWarnings(ggsignif::geom_signif(data=plotdat,
                     aes(x=get(xvar),y=Foldchange),comparisons = comparisons,
                     map_signif_level = map_signif_level,
                     test="fc_ttest",
                     tip_length = tip_length,
                     step_increase = step_increase))
  }
  if(!legend){
    p<-p+theme(legend.position = "none")
  }
  if(!is.null(xlabels)){
    p<-p+ scale_x_discrete(labels = xlabels)
  }
  return(p)
}






