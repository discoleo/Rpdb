% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proj.R
\name{proj.line3d}
\alias{proj.line3d}
\alias{proj.line3d.numeric}
\alias{proj.line3d.matrix}
\title{Projections of Objects}
\usage{
proj.line3d(p, x, y, z, ...)

\method{proj}{line3d.numeric}(p, x, y = NULL, z = NULL, ...)

\method{proj}{line3d.matrix}(p, x, y = NULL, z = NULL, ...)
}
\arguments{
\item{p}{numeric array or matrix with the x,y,z coordinates of the point or collection of points.}

\item{x, y, z}{an R object containing the coordinates of the end-points of the line segment; can be specified either as individual values or as a matrix.}

\item{\dots}{currently not used.}
}
\value{
The \code{proj.line3d} function returns an object
of class \sQuote{proj} containing the 3-dimensional coordinates of
the projected point, as well as a scalar value representing the fraction
of the path on the line segment.
}
\description{
Computes 3D Projections
}
\details{
The purpose of the \sQuote{proj} helper functions is to compute
the spatial projections of various objects on other geometrical structures.
Currently, only projection of points on a line segment are implemented.

The methods for the \sQuote{proj.line3d} function compute the projections
of 1 point or a collection of points on a line segment defined by the 2
delimiting points.

The \sQuote{numeric} method projects a point, while the \sQuote{matrix} method
projects a collection of points on another geometric object.
}
\examples{
p = c(1,2,3)
line = matrix(c(0,5,2,3,1,4), nrow = 2)
proj.line3d(p, line)

}
\keyword{classes}
\keyword{manip}
