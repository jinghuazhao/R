\name{gchol.bdsmatrix-class}
\docType{class}
\alias{gchol.bdsmatrix-class}
\alias{\%*\%,matrix,gchol.bdsmatrix-method}
\alias{\%*\%,numeric,gchol.bdsmatrix-method}
\alias{[,gchol.bdsmatrix-method}
\alias{coerce,gchol.bdsmatrix,matrix-method}
\alias{diag,gchol.bdsmatrix-method}
\alias{dim,gchol.bdsmatrix-method}
\alias{show,gchol.bdsmatrix-method}

\title{Class "gchol.bdsmatrix" }
\description{ an class generated from gchol(bdsmatrix object) }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gchol.bdsmatrix", ...)}.
    or gchol(bdsmatrix object)
}
\section{Slots}{
  \describe{
    \item{\code{blocksize}:}{Object of class \code{"integer"} ~~ }
    \item{\code{blocks}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{rmat}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{rank}:}{Object of class \code{"integer"} ~~ }
    \item{\code{.Dim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{.Dimnames}:}{Object of class \code{"list or NULL"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "matrix", y = "gchol.bdsmatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "numeric", y = "gchol.bdsmatrix")}: ... }
    \item{[}{\code{signature(x = "gchol.bdsmatrix")}: ... }
    \item{coerce}{\code{signature(from = "gchol.bdsmatrix", to = "matrix")}: ... }
    \item{diag}{\code{signature(x = "gchol.bdsmatrix")}: ... }
    \item{dim}{\code{signature(x = "gchol.bdsmatrix")}: ... }
    \item{show}{\code{signature(object = "gchol.bdsmatrix")}: ... }
  }
}

\keyword{classes}
