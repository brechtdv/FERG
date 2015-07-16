\name{getDALY-region}
\alias{getDALY-region}
\alias{getDALY_per_region}
\alias{getDALY_per_subregion}

\title{Aggregate DALYs per WHO Region and Subregion}

\description{Functions to aggregate DALYs across each of the 6 WHO regions and each of the 14 WHO subregions.}

\usage{
getDALY_per_region(agent, date_DALY, fb = FALSE)
getDALY_per_subregion(agent, date_DALY, fb = FALSE)
}

\arguments{
  \item{agent}{name of the agent}
  \item{date_DALY}{date of DALY calculation per country}
  \item{fb}{logical flag; should foodborne DALYs be used?}
}

\value{A list of \sQuote{\code{DALYrun}} objects per (sub)region and incidences per (sub)region.}

\seealso{\code{\link{getDALY-country}}}

\keyword{DALY}
