\name{imputation-functions}
\alias{imputation-functions}
\alias{impute}
\alias{merge}
\alias{save_merged_db}

\title{Main Imputation Functions}

\description{Functions to impute dataset; merge raw and imputed dataset; and save merged dataset.}

\usage{
impute(inc, WHOsub, model = "lognormal.txt", inits = NULL, 
    burnin = 5000, update = 5000)

merge(db, db_imputed)

save_merged_db(db_merged, file, sheet)
}

\arguments{
  \item{inc}{vector of mean raw incidences}
  \item{WHOsub}{vector of WHO subregions}
  \item{model}{file name of imputation model}
  \item{inits}{initial values for model parameters}
  \item{burnin}{number of burn-in samples}
  \item{update}{number of retained samples}
  \item{db}{raw data frame}
  \item{db_imputed}{imputed data frame}
  \item{db_merged}{merged data frame}
  \item{file}{file name}
  \item{sheet}{sheet name}
}

\value{\code{impute} returns the imputed dataset; \code{merge} returns the merged dataset; \code{save_merged_db} uses \pkg{XLConnect} to save the merged dataset to Excel.}

\seealso{\code{\link{imputation-reports}} for functions that report on imputation activity.}

\keyword{DataAdjust}
