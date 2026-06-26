\name{bdsmatrix-methods}
\docType{methods}
\title{Methods for bdsmatrix objects}
\alias{[,bdsmatrix,ANY,ANY,ANY-method}
\alias{[,gchol.bdsmatrix,ANY,ANY,ANY-method}
\alias{\%*\%,matrix,bdsmatrix-method}
\alias{\%*\%,numeric,bdsmatrix-method}
\alias{Math,bdsmatrix-method}
\alias{Math2,bdsmatrix-method}
\alias{Ops,bdsmatrix,bdsmatrix-method}
\alias{Ops,bdsmatrix,numeric-method}
\alias{Ops,numeric,bdsmatrix-method}
\alias{all,bdsmatrix-method}
\alias{any,bdsmatrix-method}
\alias{diag,bdsmatrix-method}
\alias{dim,bdsmatrix-method}
\alias{dimnames,bdsmatrix-method}
\alias{dimnames<-,bdsmatrix-method}
\alias{dimnames<-,bdsmatrix,ANY-method}
\alias{max,bdsmatrix-method}
\alias{min,bdsmatrix-method}
\alias{prod,bdsmatrix-method}
\alias{range,bdsmatrix-method}
\alias{sum,bdsmatrix-method}
\alias{show,bdsmatrix-method}
\alias{unique,bdsmatrix,missing-method}

\usage{
\S4method{[}{bdsmatrix,ANY,ANY,ANY}(x, i, j, ..., drop = TRUE)
\S4method{[}{gchol.bdsmatrix,ANY,ANY,ANY}(x, i, j, ..., drop = TRUE)
\S4method{dimnames}{bdsmatrix}(x) <- value
\S4method{\%*\%}{matrix,bdsmatrix}(x, y)
\S4method{\%*\%}{numeric,bdsmatrix}(x, y)
\S4method{Math}{bdsmatrix}(x)
\S4method{Math2}{bdsmatrix}(x, digits)
\S4method{Ops}{bdsmatrix,bdsmatrix}(e1, e2)
\S4method{Ops}{bdsmatrix,numeric}(e1, e2)
\S4method{Ops}{numeric,bdsmatrix}(e1, e2)
\S4method{all}{bdsmatrix}(x, ..., na.rm = FALSE)
\S4method{any}{bdsmatrix}(x, ..., na.rm = FALSE)
\S4method{diag}{bdsmatrix}(x, nrow = 1, ncol = 1)
\S4method{dim}{bdsmatrix}(x)
\S4method{dimnames}{bdsmatrix}(x)
\S4method{max}{bdsmatrix}(x, ..., na.rm = FALSE)
\S4method{min}{bdsmatrix}(x, ..., na.rm = FALSE)
\S4method{prod}{bdsmatrix}(x, ..., na.rm = FALSE)
\S4method{range}{bdsmatrix}(x, ..., na.rm = FALSE)
\S4method{sum}{bdsmatrix}(x, ..., na.rm = FALSE)
\S4method{show}{bdsmatrix}(object)
\S4method{unique}{bdsmatrix,missing}(x, incomparables = FALSE, ...)
}
 
\description{
S4 methods for bdsmatrix objects.
}

\arguments{
  \item{x}{bdsmatrix object}
  \item{i}{row indices}
  \item{j}{column indices}
  \item{...}{Additional arguments passed to or from other methods.}
  \item{drop}{logical, whether to drop dimensions}
  \item{y}{...}
  \item{digits}{...}
  \item{e1, e2}{operands}
  \item{na.rm}{logical}
  \item{nrow, ncol}{dimensions}
  \item{object}{S4 object}
  \item{incomparables}{a vector of values that cannot be compared. See base \code{unique}.}
  \item{value}{a valid list of dimnames (length 2) or NULL}
}

\keyword{methods}
\keyword{array}
