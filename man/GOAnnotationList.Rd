% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/GOAnnotationList.R
\name{GOAnnotationList}
\alias{GOAnnotationList}
\title{GOAnnotationList}
\usage{
GOAnnotationList(proteins, type, ...)
}
\arguments{
\item{proteins}{The list of proteins for which the annotation should be produced.}

\item{type}{The branch of the Gene Ontology to use, either "MF" for molecular function,
"BP" for biological process, or "CC" for cellular container}
}
\description{
\code{GOAnnotationList} retrieves and manages dense GO annotations
for a list of proteins. By dense it means that it doesn't only shows the tip/leaf
of the gene ontology to which a protein is associated, but includes all the is_a
and part_of ancestors for that term being assigned to the protein.
}
\details{
The object is a sub-class of AnnotationList.
}

