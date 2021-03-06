% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{bar_plot}
\alias{bar_plot}
\title{bar_plot}
\usage{
bar_plot(
  qpcr,
  xvar = "Condition",
  baseline_samples,
  genes,
  ND_y_nudge = 0.5,
  comparisons = NULL,
  xlabels = NULL,
  linear = T,
  y_breaks = NULL,
  map_signif_level = F,
  legend = T,
  step_increase = 0,
  tip_length = 0,
  bar_fill,
  colors = c("#bfbfbfbf", "#666666")
)
}
\arguments{
\item{qpcr}{An instance of the qpcr class}

\item{xvar}{Variable to be displayed on the x-axis. This needs to be a column n the annotation data.}

\item{baseline_samples}{Samples that should be used as baseline for foldchange calculation.}

\item{genes}{Vector of the genes that will be included in the plot.}

\item{ND_y_nudge}{Nudge value for the display of labels for undetected values.}

\item{comparisons}{List of comparisons that will be evaluated with the indicated statistical test. Each element of this list should contain a vector of length two, containing two values of the variable displayed on the x-axis.}

\item{xlabels}{Named vector containing aliases for the values of the x-axis variable. These aliases will be used for plotting.}

\item{linear}{Logical indicating whether the expression values should be in a linear range. If FALSE, log2 values will be plotted.}

\item{y_breaks}{Numerical vector with custom y-axis breaks.}

\item{map_signif_level}{Logical indicating whether to use asterisks to denote significance. If FALSE, p-values will be displayed.}

\item{legend}{Logical indicating whether a legend needs to be displayed.}

\item{step_increase}{Offset for each significance line.}

\item{tip_length}{Tip length of each significance line.}

\item{bar_fill}{Name of the variable, based on which the bars will be colored.}

\item{colors}{Colors to fill the bars with.}
}
\description{
Create bar plot
}
\details{
\code{bar_plot} makes a bar plot of the expression data
}
