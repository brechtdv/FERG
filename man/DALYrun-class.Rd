\name{DALYrun-class}
\Rdversion{1.1}
\docType{class}
\alias{DALYrun-class}

\title{Class \code{"DALYrun"}}
\description{The \code{"DALYmodel"} class represents parameterized DALY calculation models.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DALYrun", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{samples}:}{Object of class \code{"list"} ~~ List of sampled values, per node and data element.}
    \item{\code{K}:}{Object of class \code{"numeric"} ~~ Age weighting constant.}
    \item{\code{r}:}{Object of class \code{"numeric"} ~~ Discount rate.}
    \item{\code{disease}:}{Object of class \code{"character"} ~~ Disease name.}
    \item{\code{population}:}{Object of class \code{"matrix"} ~~ Population size.}
    \item{\code{LE_table}:}{Object of class \code{"matrix"} ~~ Life Expectancy table.}
    \item{\code{model}:}{Object of class \code{"list"} ~~ Disease model.}
    \item{\code{data}:}{Object of class \code{"list"} ~~ List of \code{"DALYdata"} objects.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{DALYmodel}"}, directly.
}

\seealso{
\code{\linkS4class{DALYdata}}, \code{\linkS4class{DALYmodel}}
}

\keyword{classes}
