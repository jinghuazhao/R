\name{gchol-class}
\docType{class}
\alias{gchol-class}
\alias{coerce,gchol,matrix-method}
\alias{diag,gchol-method}
\alias{show,gchol-method}

\title{Class "gchol" }
\description{  a class of gchol }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gchol", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{.Dim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{.Dimnames}:}{Object of class \code{"list or NULL"} ~~ }
    \item{\code{rank}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "gchol", to = "matrix")}: ... }
    \item{diag}{\code{signature(x = "gchol")}: ... }
    \item{show}{\code{signature(object = "gchol")}: ... }
  }
}

\keyword{classes}
