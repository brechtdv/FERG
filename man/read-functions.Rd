\name{read-functions}
\alias{read-functions}
\alias{readDatabase}
\alias{check_XL}
\alias{read_XL}

\title{Functions to Read, Check and Transform Excel Database}

\description{Read Excel database using \pkg{XLConnect}; check if sheets are correctly ordered; and transform database to 'DB' list.}

\usage{
readDatabase(file = NULL)

check_XL(XL)

read_XL(XL_wb, DM, KEY)
}

\arguments{
  \item{file}{file name}
  \item{XL_wb}{XL workbook}
  \item{XL}{XL workbook}
  \item{DM}{Disease Model}
  \item{KEY}{Key sheet}
}

\value{\code{readDatabase} returns an \pkg{XLConnect} workbook object; \code{check_XL} has no return value; and \code{read_XL} returns a 'DB' list.}

\keyword{Excel2R}
