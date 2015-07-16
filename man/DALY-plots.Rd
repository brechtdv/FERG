\name{DALY-plots}
\alias{DALY-plots}
\alias{line_plot}
\alias{region_barplot}
\alias{full_barplot}
\alias{scatter_plot}

\title{DALY Plots}

\description{Functions to generate various plots of DALY results.}

\usage{
## Uncertainty intervals
line_plot(agents, samples, asc = TRUE, names = NULL,
          ylim = NULL, ylab = "DALY (global)", ylog = TRUE,
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          col = NULL) 
        
## Barplots by (sub)region
region_barplot(x, names = NULL, xlab = NULL, ylab = NULL,
               col = NULL, nrow = 1)

## Barplots of YLD/YLL or age groups
full_barplot(..., labels = NULL, names = NULL,
             sort = c("none", "asc", "desc"), xlab = "Agent",
             ylab = "Proportion", col = NULL, nrow = 1,
             axis.text.x = element_text())

## Scatter plot of population versus patient burden
scatter_plot(x_val, y_val, fb = FALSE, names = NULL, col = NULL,
             xlog = TRUE, ylog = TRUE, labels = FALSE)
}

\arguments{
  \item{agents}{list of agents to be included in plot}
  \item{samples}{list of samples to be plotted}
  \item{asc}{logical flag; should ascending order be used?}
  \item{names}{vector of agent names}
  \item{ylim}{y-axis limits}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{xlog}{logical flag; should x-axis be shown on log scale?}
  \item{ylog}{logical flag; should y-axis be shown on log scale?}
  \item{axis.text.x}{x axis tick labels}
  \item{col}{vector of colours}
  
  \item{x}{matrix of point estimates by hazard and (sub)region}
  \item{nrow}{number of rows in legend}
  
  \item{\dots}{list of samples, one for each category}
  \item{labels}{category labels}
  \item{sort}{sorting option; should be any of \code{c("none", "asc", "desc")}}
  
  \item{x_val}{matrix of point estimate and uncertainty interval, per agent}
  \item{y_val}{matrix of point estimate and uncertainty interval, per agent}
  \item{fb}{logical flag; are foodborne results used?}

}

\value{All functions produce a \pkg{ggplot2} graphic.}

\seealso{\code{\link{DALY_map}}, \code{\link{DALY-summarize}}, \code{\link{FERG_report}}}

\keyword{report}
