% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format.R
\name{format.pdb.title}
\alias{format.pdb.title}
\alias{format.pdb.title.character}
\title{Format "Title" Field of PDB File}
\usage{
format.pdb.title(x, ...)

\method{format.pdb.title}{character}(x, ...)
}
\arguments{
\item{x}{a vector of strings;}

\item{\dots}{currently not used;}
}
\value{
properly formatted text string.
}
\description{
Formats the title field of a Protein Data Bank (PDB) file.
}
\details{
This function is mostly used internally, e.g. by \code{\link{write.pdb}},
but may be useful in other instances as well.
}
\examples{
format.pdb.title(c("Molecule 1", "is just an example"))
}
\seealso{
\code{\link{write.pdb}}
}
