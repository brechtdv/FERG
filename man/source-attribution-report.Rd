\name{source-attribution-report}
\alias{source-attribution-report}
\alias{summary_stats}
\alias{sa_plot}
\alias{sa_barplot}
\alias{sa_report}

\title{Source Attribution Reporting Functions}

\description{.}

\usage{
## Summarize samples
summary_stats(x, dig = 2, ...)

## Generate source attribution line plot
sa_plot(samples, what = "food", main = NULL)

## Generate source attribution barplot
sa_barplot(samples)

## Generate source attribution PDF report
sa_report(hazard, exp_samples, food_samples, exp_source, food_source)
}

\arguments{
  \item{x}{vector of samples}
  \item{dig}{integer indicating the number of decimal places to be used}
  \item{\dots}{other arguments to be passed to underlying functions}
  \item{samples}{list or matrix samples}
  \item{what}{item to be plotted}
  \item{main}{main plot title}
  \item{hazard}{name of the hazard}
  \item{exp_samples}{list of exposure route samples}
  \item{food_samples}{list of food item samples}
  \item{exp_source}{vector of exposure route names}
  \item{food_source}{vector of food item names}
}

\value{\code{summary_stats} returns a vector of the samples' mean, median, and 90\% and 95\% uncertainty interval. \code{sa_plot} and \code{sa_barplot} generate a \pkg{ggplot2} graph. \code{sa_report} generates a PDF report.}

\seealso{\code{\link{source-attribution-functions}}}

\keyword{SA}
