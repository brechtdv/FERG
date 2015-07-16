\name{DALYmodel-class}
\Rdversion{1.1}
\docType{class}
\alias{DALYmodel-class}

\title{Class \code{"DALYmodel"}}
\description{The \code{"DALYmodel"} class represents parameterized DALY calculation models.}

\section{Objects from the Class}{
Objects can be created by a call to setDALYmodel().
}

\section{Slots}{
  \describe{
    \item{\code{disease}:}{Object of class \code{"character"} ~~ Disease name.}
    \item{\code{population}:}{Object of class \code{"matrix"} ~~ Population size.}
    \item{\code{LE_table}:}{Object of class \code{"matrix"} ~~ Life Expectancy table.}
    \item{\code{model}:}{Object of class \code{"list"} ~~ Disease model.}
    \item{\code{data}:}{Object of class \code{"list"} ~~ List of \code{"DALYdata"} objects.}
  }
}

\seealso{
\code{\linkS4class{DALYdata}}, \code{\linkS4class{DALYrun}}
}

\keyword{classes}
