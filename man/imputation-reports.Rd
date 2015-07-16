\name{imputation-reports}
\alias{imputation-reports}
\alias{plot_missing}
\alias{summarize_imputation}
\alias{get_total_cases}
\alias{imputation_report}

\title{Functions to Report on Imputation Activity}

\description{Functions to plot a worldmap of available and missing data; to summarize imputation results; to calculate total (ie, raw+imputed) cases; and to generate a PDF report of the imputation activity.}

\usage{
plot_missing(inc, WHOsub, agent, agent_full, save = FALSE)

summarize_imputation(db, db_imputed)

get_total_cases(db_merged, pop_2010)

imputation_report(db, db_imputed, db_merged, pop, agent, agent_full,
                  tab1_options = NULL, tab2_options = NULL)
}

\arguments{
  \item{inc}{vector of mean raw incidences}
  \item{WHOsub}{vector of WHO subregions}
  \item{agent}{agent name, abbreviated}
  \item{agent_full}{agent name, full}
  \item{save}{logical flag; should image be saved?}

  \item{db_merged}{merged data frame}
  \item{pop_2010}{vector of population size by country (*1e5)}
    
  \item{db}{raw data frame}
  \item{db_imputed}{imputed data frame}
  \item{pop}{vector of population size by country (*1e5)}
  \item{tab1_options}{list of formatting options for tables of absolute numbers}
  \item{tab2_options}{list of formatting options for tables of rates}
}

\value{\code{get_total_cases} returns a data frame with WHO subregion; mean, median, 2.5th percentile, 97.5th percentile, and standard deviation of overall subregional number of cases; and subregional population size. \code{summarize_imputation} returns a list of regional fit; global geometric mean; within-region variance; and between-region variance.\code{plot_missing} creates a PNG file (if \code{save = FALSE}); \code{imputation_report} generates a PDF file.}

\seealso{\code{\link{imputation-functions}} for main imputation functions.}

\keyword{report}
