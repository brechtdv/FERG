\name{getDALY-country}
\alias{getDALY-country}
\alias{getDALY_per_country}

\title{Calculate DALYs per Country}

\description{
Function to calculate DALYs per country.
}

\usage{
getDALY_per_country(agent, DB, DM, today, LE = c("WHO", "GBD", "CD"),
                    n_samples = 1000)
}

\arguments{
  \item{agent}{name of agent}
  \item{DB}{database}
  \item{DM}{disease model}
  \item{today}{today's date, in the format \code{yyyymmdd}}
  \item{LE}{life expectancy table to be used; defaults to "WHO"}
  \item{n_samples}{number of samples for the Monte Carlo simulations; defaults to 1000}
}

\value{An object of class \sQuote{\code{DALYrun}}.}

\seealso{\code{\link{getDALY-region}}}


\keyword{DALY}
