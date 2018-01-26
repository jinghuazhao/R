\name{bdsmatrix-class}
\docType{class}
\alias{bdsmatrix-class}
\alias{\%*\%,matrix,bdsmatrix-method}
\alias{\%*\%,numeric,bdsmatrix-method}
\alias{Math2,bdsmatrix-method}
\alias{Math,bdsmatrix-method}
\alias{Ops,bdsmatrix,numeric-method}
\alias{Ops,bdsmatrix,bdsmatrix-method}
\alias{Ops,numeric,bdsmatrix-method}
\alias{[,bdsmatrix-method}
\alias{all,bdsmatrix-method}
\alias{any,bdsmatrix-method}
\alias{coerce,bdsmatrix,matrix-method}
\alias{coerce,bdsmatrix,vector-method}
\alias{diag,bdsmatrix-method}
\alias{diag<-,bdsmatrix-method}
\alias{dim,bdsmatrix-method}
\alias{dimnames,bdsmatrix-method}
\alias{dimnames<-,bdsmatrix-method}
\alias{max,bdsmatrix-method}
\alias{min,bdsmatrix-method}
\alias{prod,bdsmatrix-method}
\alias{range,bdsmatrix-method}
\alias{show,bdsmatrix-method}
\alias{sum,bdsmatrix-method}
\alias{unique,bdsmatrix,missing-method}

\title{Class "bdsmatrix" }
\description{ A bdsmarix class }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bdsmatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{blocksize}:}{Object of class \code{"integer"} vector of sizes for the matrices on the diagonal }
    \item{\code{blocks}:}{Object of class \code{"numeric"} contents of the diagonal blocks, strung out as a vector }
    \item{\code{rmat}:}{Object of class \code{"matrix"} the dense portion of the matrix, forming a right and lower border }
    \item{\code{offdiag}:}{Object of class \code{"numeric"} 0s }
    \item{\code{.Dim}:}{Object of class \code{"integer"} dimensions }
    \item{\code{.Dimnames}:}{Object of class \code{"list or NULL"} a list of dimension names for the matrix }
  }
}
\section{Extends}{
Class \code{"matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "matrix", y = "bdsmatrix")}:  }
    \item{\%*\%}{\code{signature(x = "numeric", y = "bdsmatrix")}:  }
    \item{Math2}{\code{signature(x = "bdsmatrix")}:  }
    \item{Math}{\code{signature(x = "bdsmatrix")}:  }
    \item{Ops}{\code{signature(e1 = "bdsmatrix", e2 = "numeric")}:  }
    \item{Ops}{\code{signature(e1 = "bdsmatrix", e2 = "bdsmatrix")}:  }
    \item{Ops}{\code{signature(e1 = "numeric", e2 = "bdsmatrix")}:  }
    \item{[}{\code{signature(x = "bdsmatrix")}:  }
    \item{all}{\code{signature(x = "bdsmatrix")}:  }
    \item{any}{\code{signature(x = "bdsmatrix")}:  }
    \item{coerce}{\code{signature(from = "bdsmatrix", to = "matrix")}:  }
    \item{coerce}{\code{signature(from = "bdsmatrix", to = "vector")}:  }
    \item{diag}{\code{signature(x = "bdsmatrix")}:  }
    \item{diag<-}{\code{signature(x = "bdsmatrix")}:  }
    \item{dim}{\code{signature(x = "bdsmatrix")}:  }
    \item{dimnames}{\code{signature(x = "bdsmatrix")}:  }
    \item{dimnames<-}{\code{signature(x = "bdsmatrix")}:  }
    \item{max}{\code{signature(x = "bdsmatrix")}:  }
    \item{min}{\code{signature(x = "bdsmatrix")}:  }
    \item{prod}{\code{signature(x = "bdsmatrix")}:  }
    \item{range}{\code{signature(x = "bdsmatrix")}:  }
    \item{show}{\code{signature(object = "bdsmatrix")}:  }
    \item{sum}{\code{signature(x = "bdsmatrix")}:  }
    \item{unique}{\code{signature(x = "bdsmatrix", incomparables = "missing")}:  }
  }
}

\keyword{classes}
