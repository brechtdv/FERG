\name{DALYdata-class}
\Rdversion{1.1}
\docType{class}
\alias{DALYdata-class}

\title{Class \code{"DALYdata"}}
\description{The \code{"DALYdata"} class represents data elements needed in DALY calculations.}
\section{Objects from the Class}{
Objects can be created by a call to setDALYdata().
}
\section{Slots}{
  \describe{
    \item{\code{dist}:}{Object of class \code{"character"} ~~ Name of the distribution.}
    \item{\code{strat}:}{Object of class \code{"character"} ~~ Name of the stratification.}
    \item{\code{param}:}{Object of class \code{"par"} ~~ Actual data parameters.}
  }
}

\seealso{
\code{\linkS4class{DALYmodel}}, \code{\linkS4class{DALYrun}}
}

\keyword{classes}
