\name{bdsmatrix-class}
\docType{class}
\alias{bdsmatrix-class}

\title{Class "bdsmatrix"}

\description{
A sparse block-structured matrix class used in kinship calculations.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bdsmatrix", ...)}.
}

\section{Slots}{
  \describe{
    \item{blocksize}{integer vector of block sizes}
    \item{blocks}{numeric vector of diagonal block entries}
    \item{rmat}{matrix, dense border region}
    \item{offdiag}{numeric, off-diagonal constant values}
    \item{.Dim}{integer dimensions}
    \item{.Dimnames}{list or NULL dimension names}
  }
}

\section{Extends}{
Class \code{"matrix"}, directly.
}

\keyword{classes}
