\name{imputation-extract}
\alias{imputation-extract}
\alias{extractIncidence}
\alias{check_db}

\title{Extract Incidence Data for Imputation}

\description{Get incidence data from Excel database file and check structure.}

\usage{
extractIncidence(file, sheet)
check_db(db)}

\arguments{
  \item{file}{file name}
  \item{sheet}{sheet name}
  \item{db}{data frame returned by \code{extractIncidence}}
}

\value{\code{extractIncidence} returns a data frame with country name, WHO subregion, mean incidence, and incidence quantiles. \code{check_db} has no return value.}

\seealso{\code{\link{readDatabase}}}

\keyword{Excel2R}
